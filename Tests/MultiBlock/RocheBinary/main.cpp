// 2D contact-binary (Roche lobe + bridge) example using AMReX NonLocalBC on
// a multi-block curvilinear grid.
//
// Layout: two fixed "stars" (masked holes, not evolved) sit side by side on
// the x-axis, each wrapped in four curvilinear blocks covering one quadrant
// of the annulus a < r < r_ann1 (r_ann2) around it. A ninth curvilinear
// the gap between the two annuli across the L1 region:
//
//                 +-------+                 +-------+
//                 |  N1   |                 |  N2   |
//         +-------+ .---. +-------+-+-----+-+ .---. +-------+
//         |  W1   |(st 1)|  E1   | BRIDGE  |  W2  |(st 2)|  E2  |
//         +-------+ '---' +-------+-+-----+-+ '---' +-------+
//                 |  S1   |                 |  S2   |
//                 +-------+                 +-------+
//
// Polar blocks (E1,N1,W1,S1 around star 1 at x = x1 < 0, E2,...,S2 around
// star 2 at x = x2 > 0) use the StarDisk mapping translated to their star:
//
//   phi(i) = -pi/4 + (i/n_phi)*(pi/2),  rho(j) = a + (r_annK-a)*(j/n_r)^stretch
//   x = c + rho*(cos,sin)(theta_b + phi),   theta_b in {0, pi/2, pi, 3pi/2}
//
// with r_annK = r_ann1 (star 1) resp. r_ann2 (star 2); the radii may differ
// so each star's Roche lobe can be contained (mandatory for q != 1).
//
// with logical (i = tangential CCW, j = radial). The bridge block is a Coons
// (transfinite) patch with logical (i = along the arc, bottom -> top;
// j = corridor axis, star 1 -> star 2): its j-low edge is the outer arc of
// E1 (angles -45..+45 about star 1), its j-high edge the outer arc of W2
// (angles 135..225 about star 2), and its i-low/i-high edges are straight
// lines connecting the arc endpoints. All mappings are left-handed, so the
// same vertex-based metrics code (abs shoelace volumes, A = (-ey, ex) face
// areas) serves every block.
//
// Seams (20 one-sided NonLocalBC fills, no periodicity anywhere):
//   * 8 per star: tangential E<->N<->W<->S<->E, offset-only mappings
//     (identical to StarDisk);
//   * bridge j-low  <-> E1 j-high: offset-only (both i run bottom -> top);
//   * bridge j-high <-> W2 j-high: i-reversed mapping (W2's i runs
//     top -> bottom), dtos sign = (-1, 1) -- the non-trivial seam.
//
// Physics: isothermal Euler (rho, rho*u) with p = cs^2 rho, first-order
// finite volume with Rusanov fluxes on the curvilinear mesh, explicit Euler
// stepping. The simulation runs in the binary's corotating frame (both stars
// stationary) with Omega^2 = G (m1+m2) / sep^3 about the center of mass:
//   source = -rho grad(Phi) - 2 Omega z-hat x (rho u)
//   Phi    = -G m1/r1 - G m2/r2 - 0.5 Omega^2 (x^2 + y^2).
// The initial gas is the exact hydrostatic isothermal atmosphere of the
// Roche potential, rho = rho_l1 * exp(-(Phi - Phi_L1)/cs^2) (floored), at
// rest in the corotating frame: an overcontact (common-envelope) configura-
// tion with a gas bridge through L1. roche.perturb boosts the density on the
// star-1 side to drive mass transfer through the bridge. Setting
// roche.ic_mode = "semi_detached" instead leaves the accretor side
// (x >= x_l1) at the floor: the L1-plane discontinuity is then the initial
// condition, modelling a semi-detached binary in which only the donor fills
// its Roche lobe (the mass-transfer stream is born at L1).
// roche.counter_rotate=1 overrides the rest condition: the gas starts at
// inertial rest, which in the corotating frame is solid-body counter-rotation
// v = omega*(y,-x) (Mach ~3.5 at the outer edge). The unbalanced Coriolis
// force launches shocks and the circulation brakes back to corotation over
// ~10 orbits. This is not a discrete equilibrium, so run it with
// outer_bc=closed (the reservoir would pin the open edges to counter-rotation).
//
// Boundary conditions: reflecting wall at each star surface (polar j-low).
// Open edges (polar j-high, bridge i-edges) use a far-field reservoir: ghost
// cells hold the initial hydrostatic state, so the boundary cannot chase its
// own collapse (zero-gradient copies feed a runaway in the tenuous gas).
//
// Build: cd Tests/MultiBlock/RocheBinary && make -j
// Run:   mpirun -n 2 ./main2d.gnu.MPI.ex [roche.name=value ...]
// Output: RocheBinary/{e1,n1,w1,s1,e2,n2,w2,s2,br}/plt#### (rho,momx,momy,x,y).

#include "AMReX_NonLocalBC.H"

#include "AMReX.H"
#include "AMReX_AmrCore.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"
#include "AMReX_PlotFileUtil.H"
#include "AMReX_REAL.H"
#include "AMReX_Reduce.H"

#include <cmath>
#include <limits>
#include <string>
#include <fstream>
#include <iomanip>

static_assert(AMREX_SPACEDIM == 2, "RocheBinary is a 2D-only test");

using namespace amrex;

void MyMain();

int main(int argc, char** argv) {
    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, {},
                      std::cout, std::cerr,
                      [](const char* msg) { throw std::runtime_error(msg); });
    MyMain();
    amrex::Finalize();
}

namespace {

constexpr int ncomp = 3;           // rho, mom_x, mom_y
constexpr int URHO = 0;
constexpr int UMX  = 1;
constexpr int UMY  = 2;

constexpr int nghost = 1;          // first-order: one ghost layer

enum idirs { ix, iy };

// Initial-condition mode (host-parsed from roche.ic_mode into this int, since
// RocheParams is copied by value into device lambdas and cannot hold a string).
enum IcMode { IC_OVERCONTACT = 0, IC_SEMI_DETACHED = 1 };

struct RocheParams {
    Real m1        = 1.0;   // mass of star 1 (at x1 < 0)
    Real m2        = 1.0;   // mass of star 2 (at x2 > 0)
    Real sep       = 4.0;   // stellar separation D
    Real a         = 1.0;   // star radius (both stars)
    Real r_ann1    = 1.6;   // outer radius of star 1's annulus mesh
    Real r_ann2    = 1.6;   // outer radius of star 2's annulus mesh
    Real stretch   = 1.0;   // radial grading exponent (>1 clusters inward)
    Real cs        = 0.25;  // isothermal sound speed
    Real rho_l1    = 1.e-2; // IC density at L1 (sets the overcontact depth)
    Real rho_l1_1  = -1.0;  // <0 -> use rho_l1 (star-1 lobe normalization)
    Real rho_l1_2  = -1.0;  // <0 -> use rho_l1 (star-2 lobe normalization)
    Real neck_blend = 1.0;  // multiplier on the L1-neck blend width
    Real rho_floor = 1.e-8;
    Real soft      = 0.0;   // gravitational softening length
    Real cfl       = 0.5;
    Real vmax      = 0.0;   // >0 caps |v| (semi-detached vacuum-jet fix)
    Real perturb   = 0.0;   // relative density boost on the star-1 side
    int  ic_mode   = IC_OVERCONTACT;  // IC_OVERCONTACT | IC_SEMI_DETACHED
    int  outer_closed = 0;  // 1: reflecting wall at outer edges (closed system
                            //    for the differential-depth IC, else reservoir)
    int  counter_rotate = 0;// 1: IC gas at inertial rest, which in the corotating
                            //    frame is solid-body counter-rotation v = omega*
                            //    (y,-x); the unbalanced Coriolis force drives
                            //    shocks and the circulation brakes to corotation
                            //    over ~10 orbits. Not a discrete equilibrium:
                            //    run with outer_bc=closed (the reservoir would
                            //    pin the boundary to the counter-rotating state).
    int  n_phi     = 64;    // tangential cells per polar block (= bridge i)
    int  n_r       = 64;    // radial cells per polar block
    int  n_bridge  = 64;    // bridge cells along the corridor axis (j)
    // Derived quantities, filled by derive():
    Real omega     = 0.0;   // binary angular frequency, G = 1
    Real x1        = 0.0;   // star 1 position (COM at origin)
    Real x2        = 0.0;   // star 2 position
    Real x_l1      = 0.0;   // L1 location on the axis
    Real phi_l1    = 0.0;   // Roche potential at L1
    Real phi_xx_l1 = 0.0;   // d^2 Phi/dx^2 at L1 (neck-blend width scale)
    Real blend_width = 0.0; // L1-neck normalization-blend half-width
};

// Effective Roche potential in the corotating frame (COM at origin).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real potential (Real x, Real y, RocheParams const& p) {
    const Real dx1 = x - p.x1, dx2 = x - p.x2;
    const Real r1 = std::sqrt(dx1 * dx1 + y * y + p.soft * p.soft);
    const Real r2 = std::sqrt(dx2 * dx2 + y * y + p.soft * p.soft);
    return -p.m1 / r1 - p.m2 / r2
         - 0.5_rt * p.omega * p.omega * (x * x + y * y);
}

// dPhi/dx on the binary axis (y = 0). Host-only, used to locate L1.
Real potential_dx (Real x, RocheParams const& p) {
    const Real dx1 = x - p.x1, dx2 = x - p.x2;
    const Real s1 = dx1 * dx1 + p.soft * p.soft;
    const Real s2 = dx2 * dx2 + p.soft * p.soft;
    return p.m1 * dx1 / (s1 * std::sqrt(s1))
         + p.m2 * dx2 / (s2 * std::sqrt(s2))
         - p.omega * p.omega * x;
}
// Reflecting-wall ghost fill: mirror the interior state across a face with
// edge direction (ex, ey) on this left-handed grid, flipping the normal
// momentum component (free slip). Used for the star surface and, in
// closed-outer-boundary mode, the outer edges.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wall_ghost (Real ex, Real ey, Real rho_in, Real mx_in, Real my_in,
                 Real& rho_g, Real& mx_g, Real& my_g) {
    const Real anorm = std::sqrt(ey * ey + ex * ex);
    const Real nx = -ey / anorm;
    const Real ny =  ex / anorm;
    const Real mn = mx_in * nx + my_in * ny;
    rho_g = rho_in;
    mx_g  = mx_in - 2.0_rt * mn * nx;
    my_g  = my_in - 2.0_rt * mn * ny;
}

void derive (RocheParams& p) {
    p.omega = std::sqrt((p.m1 + p.m2) / (p.sep * p.sep * p.sep));  // G = 1
    p.x1 = -p.sep * p.m2 / (p.m1 + p.m2);
    p.x2 =  p.sep * p.m1 / (p.m1 + p.m2);
    // L1 is the maximum of Phi on the axis between the two stars; bracket it
    // between the star surfaces and bisect on dPhi/dx.
    Real lo = p.x1 + p.a, hi = p.x2 - p.a;
    Real flo = potential_dx(lo, p), fhi = potential_dx(hi, p);
    if (!(flo > 0.0_rt && fhi < 0.0_rt)) {
        amrex::Abort("RocheBinary: cannot bracket L1; check masses/sep/a/soft");
    }
    for (int it = 0; it < 200; ++it) {
        const Real mid = 0.5_rt * (lo + hi);
        if (potential_dx(mid, p) > 0.0_rt) lo = mid; else hi = mid;
    }
    p.x_l1 = 0.5_rt * (lo + hi);
    p.phi_l1 = potential(p.x_l1, 0.0_rt, p);
    // Resolve per-lobe normalizations (sentinel -> rho_l1; backward compat:
    // both unset => rho_l1_1 == rho_l1_2 == rho_l1, IC unchanged).
    if (p.rho_l1_1 < 0.0_rt) p.rho_l1_1 = p.rho_l1;
    if (p.rho_l1_2 < 0.0_rt) p.rho_l1_2 = p.rho_l1;
    // Curvature of Phi along the corridor axis at L1 (a saddle: max along x).
    // Sets the neck-blend width for the differential-depth IC to ~cs^2 in
    // potential: w = cs*sqrt(2/|Phi_xx|) => DeltaPhi(w) = cs^2.
    {
        const Real hb = 1e-5_rt * p.sep;
        const Real d2 = (potential_dx(p.x_l1 + hb, p)
                       - potential_dx(p.x_l1 - hb, p)) / (2.0_rt * hb);
        p.phi_xx_l1 = d2;
        p.blend_width = (std::abs(d2) > 0.0_rt)
            ? p.neck_blend * p.cs * std::sqrt(2.0_rt / std::abs(d2))
            : 1.0_rt;
    }
    // Geometry constraints: the annuli must not overlap, and the bridge
    // corridor (x1 + r_ann1, x2 - r_ann2) must contain L1.
    if (p.r_ann1 + p.r_ann2 >= p.sep) {
        amrex::Abort("RocheBinary: r_ann1 + r_ann2 >= sep (annuli overlap)");
    }
    if (!(p.x1 + p.r_ann1 < p.x_l1 && p.x_l1 < p.x2 - p.r_ann2)) {
        amrex::Abort("RocheBinary: L1 not inside the bridge corridor; "
                     "adjust r_ann1/r_ann2");
    }
}

enum class BlockKind { Polar, Bridge };

struct BlockGeom {
    BlockKind kind;
    Real cx      = 0.0;  // polar: x of star center
    Real cy      = 0.0;  // polar: y of star center
    Real theta_b = 0.0;  // polar: orientation angle
    Real r_ann   = 1.6;  // polar: outer radius of this block's annulus
    bool jhi_seam = false;  // polar: j-high edge is a seam, not outflow
};

// Physical position of vertex (i, j).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuArray<Real, 2> vertex_pos (int i, int j, BlockGeom const& g,
                              RocheParams const& p) {
    if (g.kind == BlockKind::Polar) {
        const Real dphi = (0.5_rt * M_PI) / p.n_phi;
        const Real phi  = -0.25_rt * M_PI + i * dphi;
        const Real s    = Real(j) / Real(p.n_r);
        const Real rho  = p.a + (g.r_ann - p.a) * std::pow(s, p.stretch);
        const Real ang  = g.theta_b + phi;
        return {g.cx + rho * std::cos(ang), g.cy + rho * std::sin(ang)};
    }
    // Bridge: Coons patch. s runs along the arcs (bottom -> top, same
    // direction as E1's i at the shared edge), t along the corridor.
    const Real s = Real(i) / Real(p.n_phi);
    const Real t = Real(j) / Real(p.n_bridge);
    const Real thl = -0.25_rt * M_PI + s * (0.5_rt * M_PI);  // -45..45 deg
    const Real thr =  1.25_rt * M_PI - s * (0.5_rt * M_PI);  // 225..135 deg
    const Real lx = p.x1 + p.r_ann1 * std::cos(thl);
    const Real ly =        p.r_ann1 * std::sin(thl);
    const Real rx = p.x2 + p.r_ann2 * std::cos(thr);
    const Real ry =        p.r_ann2 * std::sin(thr);
    // Corners: A = L(0), B = L(1), C = R(1), D = R(0).
    const Real cr1 = p.r_ann1 * 0.5_rt * std::sqrt(2.0_rt);  // r_ann1 * cos(45)
    const Real cr2 = p.r_ann2 * 0.5_rt * std::sqrt(2.0_rt);  // r_ann2 * cos(45)
    const Real ax = p.x1 + cr1, ay = -cr1;
    const Real bx = p.x1 + cr1, by =  cr1;
    const Real cx = p.x2 - cr2, cy =  cr2;
    const Real dx = p.x2 - cr2, dy = -cr2;
    // Straight bottom (A->D) and top (B->C) edges.
    const Real btm_x = (1.0_rt - t) * ax + t * dx;
    const Real btm_y = (1.0_rt - t) * ay + t * dy;
    const Real top_x = (1.0_rt - t) * bx + t * cx;
    const Real top_y = (1.0_rt - t) * by + t * cy;
    // Transfinite (Coons) blend; the bilinear term removes double-counting.
    const Real bl_x = (1.0_rt - s) * (1.0_rt - t) * ax + s * (1.0_rt - t) * bx
                    + (1.0_rt - s) * t * dx          + s * t * cx;
    const Real bl_y = (1.0_rt - s) * (1.0_rt - t) * ay + s * (1.0_rt - t) * by
                    + (1.0_rt - s) * t * dy          + s * t * cy;
    return {(1.0_rt - t) * lx + t * rx + (1.0_rt - s) * btm_x + s * top_x - bl_x,
            (1.0_rt - t) * ly + t * ry + (1.0_rt - s) * btm_y + s * top_y - bl_y};
}

// Rusanov (local Lax-Friedrichs) flux for isothermal Euler across a face with
// outward area vector (ax, ay). UL is the state inside, UR outside. The
// dissipative term acts on the deviation from the initial equilibrium state
// (U0L/U0R): in equilibrium the dissipation vanishes identically, so the
// steep hydrostatic atmosphere is not eroded by numerical diffusion on long
// times (well-balanced; plain-U dissipation evaporates the lobes in a few
// sound times).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rusanov_flux (Real const* UL, Real const* UR, Real const* U0L,
                   Real const* U0R, Real ax, Real ay, Real cs, Real* F) {
    const Real area = std::sqrt(ax * ax + ay * ay);
    const Real nx = ax / area;
    const Real ny = ay / area;
    const Real rhoL = UL[URHO], rhoR = UR[URHO];
    const Real unL = (UL[UMX] * nx + UL[UMY] * ny) / rhoL;
    const Real unR = (UR[UMX] * nx + UR[UMY] * ny) / rhoR;
    const Real pL = cs * cs * rhoL;
    const Real pR = cs * cs * rhoR;
    const Real smax = amrex::max(std::abs(unL), std::abs(unR)) + cs;
    F[URHO] = area * (0.5_rt * (rhoL * unL + rhoR * unR)
                      - 0.5_rt * smax * ((UR[URHO] - U0R[URHO])
                                       - (UL[URHO] - U0L[URHO])));
    F[UMX]  = area * (0.5_rt * (UL[UMX] * unL + UR[UMX] * unR + (pL + pR) * nx)
                      - 0.5_rt * smax * ((UR[UMX] - U0R[UMX])
                                       - (UL[UMX] - U0L[UMX])));
    F[UMY]  = area * (0.5_rt * (UL[UMY] * unL + UR[UMY] * unR + (pL + pR) * ny)
                      - 0.5_rt * smax * ((UR[UMY] - U0R[UMY])
                                       - (UL[UMY] - U0L[UMY])));
}

class RocheBlock : public AmrCore {
  public:
    RocheBlock(BlockGeom g_, Geometry const& level_0_geom, RocheParams const& p,
               AmrInfo const& amr_info = AmrInfo())
        : AmrCore(level_0_geom, amr_info), params{p}, g{g_} {
        AmrCore::InitFromScratch(0.0);
        InitData();
    }

    // Vertex coordinates, cell volumes, and the hydrostatic atmosphere.
    void InitData() {
        const RocheParams p = params;
        const BlockGeom bg = g;
        for (MFIter mfi(xyv, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Array4<Real> v = xyv.array(mfi);
            ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::ignore_unused(k);
                const auto pos = vertex_pos(i, j, bg, p);
                v(i, j, 0, 0) = pos[0];
                v(i, j, 0, 1) = pos[1];
            });
        }
        for (MFIter mfi(U, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Array4<Real> vol_a = vol.array(mfi);
            Array4<Real> u = U.array(mfi);
            Array4<Real const> v = xyv.const_array(mfi);
            ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::ignore_unused(k);
                const Real x00 = v(i, j, 0, 0),     y00 = v(i, j, 0, 1);
                const Real x10 = v(i + 1, j, 0, 0), y10 = v(i + 1, j, 0, 1);
                const Real x11 = v(i + 1, j + 1, 0, 0), y11 = v(i + 1, j + 1, 0, 1);
                const Real x01 = v(i, j + 1, 0, 0), y01 = v(i, j + 1, 0, 1);
                // Shoelace. All mappings here are left-handed, so the signed
                // area comes out negative.
                vol_a(i, j, k) = 0.5_rt * std::abs((x00 * y10 - x10 * y00)
                                                 + (x10 * y11 - x11 * y10)
                                                 + (x11 * y01 - x01 * y11)
                                                 + (x01 * y00 - x00 * y01));
                const Real xc = 0.25_rt * (x00 + x10 + x11 + x01);
                const Real yc = 0.25_rt * (y00 + y10 + y11 + y01);
                // Hydrostatic isothermal atmosphere in the Roche potential,
                // normalized to rho_l1 at L1; at rest in the corotating frame.
                const Real phi = potential(xc, yc, p);
                Real rho;
                if (p.ic_mode == IC_SEMI_DETACHED && xc >= p.x_l1) {
                    // Semi-detached IC: accretor side left at the floor; the
                    // L1-plane discontinuity is the stream's birth.
                    rho = p.rho_floor;
                } else {
                    // Per-lobe normalization: rho_l1_1 (star 1) smoothly
                    // blended to rho_l1_2 (star 2) across the L1 neck over
                    // width blend_width (~cs^2 in potential, so the bridge is
                    // shared with no hard jump). Equal values recover the
                    // single-normalization overcontact IC exactly.
                    const Real rl1 = (p.rho_l1_1 == p.rho_l1_2)
                        ? p.rho_l1
                        : (p.rho_l1_2 + 0.5_rt * (p.rho_l1_1 - p.rho_l1_2)
                                     * (1.0_rt - std::tanh((xc - p.x_l1)
                                                            / p.blend_width)));
                    rho = amrex::max(rl1 * std::exp(-(phi - p.phi_l1)
                                                    / (p.cs * p.cs)),
                                      p.rho_floor);
                    if (p.perturb != 0.0_rt && xc < p.x_l1) {
                        rho *= 1.0_rt + p.perturb;
                    }
                }
                u(i, j, k, URHO) = rho;
                if (p.counter_rotate) {
                    // Inertial-rest gas seen in the corotating frame: solid-
                    // body counter-rotation v = -Omega z x r = omega*(y, -x).
                    // Mach ~ omega*r/cs, so ~3.5 at the outer edge initially.
                    u(i, j, k, UMX) = rho * p.omega * yc;
                    u(i, j, k, UMY) = -rho * p.omega * xc;
                } else {
                    u(i, j, k, UMX) = 0.0_rt;
                    u(i, j, k, UMY) = 0.0_rt;
                }
            });
        }
        // Interior of the reference state (its ghost cells are snapshot from
        // U's pipeline-filled ghosts in MyMain -- see SnapshotInitGhosts).
        MultiFab::Copy(Uinit, U, 0, 0, ncomp, 0);
    }

    // Physical boundary conditions on the ghost cells of U. Polar blocks:
    // reflecting wall at the star (j-low); at j-high a far-field reservoir
    // (initial state), or a reflecting wall when outer_closed (closed system
    // for the differential-depth IC -- the reservoir would pin the boundary
    // to the unequal initial profile and sustain the imbalance). Bridge
    // block: same choice on both i edges (its j edges are seams). Ghost cells
    // between grids of a block belong to FillGhosts().
    void FillPhysicalBCs() {
        const Box& domain = Geom(0).Domain();
        const bool closed = params.outer_closed;
        for (MFIter mfi(U); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            Array4<Real> u = U.array(mfi);
            Array4<Real const> u0 = Uinit.const_array(mfi);
            Array4<Real const> v = xyv.const_array(mfi);
            if (g.kind == BlockKind::Polar) {
                const Box glo = amrex::adjCellLo(vbx, iy, nghost)
                    & amrex::adjCellLo(domain, iy, nghost);
                const Box ghi = amrex::adjCellHi(vbx, iy, nghost)
                    & amrex::adjCellHi(domain, iy, nghost);
                // j-low: wall. Ghost (i, j0-1) mirrors interior (i, j0) with
                // the normal momentum component flipped (free slip).
                if (glo.ok()) {
                ParallelFor(glo,
                            [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    amrex::ignore_unused(k);
                    // Outward area vector of the j-low face of cell (i, j+1),
                    // left-handed grid: A = (-ey, ex) for edge
                    // (i,j+1)->(i+1,j+1).
                    const Real ex = v(i + 1, j + 1, 0, 0) - v(i, j + 1, 0, 0);
                    const Real ey = v(i + 1, j + 1, 0, 1) - v(i, j + 1, 0, 1);
                    const Real anorm = std::sqrt(ey * ey + ex * ex);
                    const Real nx = -ey / anorm;
                    const Real ny =  ex / anorm;
                    const Real mx = u(i, j + 1, k, UMX);
                    const Real my = u(i, j + 1, k, UMY);
                    const Real mn = mx * nx + my * ny;
                    u(i, j, k, URHO) = u(i, j + 1, k, URHO);
                    u(i, j, k, UMX)  = mx - 2.0_rt * mn * nx;
                    u(i, j, k, UMY)  = my - 2.0_rt * mn * ny;
                });
                }
                // j-high (unless a seam): reservoir (initial state), or a
                // reflecting wall when closed. Skipped on seams either way:
                // overwriting seam-filled ghosts breaks flux cancellation.
                if (!g.jhi_seam && ghi.ok()) {
                ParallelFor(ghi,
                            [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    amrex::ignore_unused(k);
                    if (closed) {
                        const Real ex = v(i + 1, j, 0, 0) - v(i, j, 0, 0);
                        const Real ey = v(i + 1, j, 0, 1) - v(i, j, 0, 1);
                        Real rg, mxg, myg;
                        wall_ghost(ex, ey, u(i, j - 1, k, URHO),
                                   u(i, j - 1, k, UMX), u(i, j - 1, k, UMY),
                                   rg, mxg, myg);
                        u(i, j, k, URHO) = rg;
                        u(i, j, k, UMX)  = mxg;
                        u(i, j, k, UMY)  = myg;
                    } else {
                        for (int n = 0; n < ncomp; ++n) {
                            u(i, j, k, n) = u0(i, j - 1, k, n);
                        }
                    }
                });
                }
            } else {
                // Bridge: reservoir (or reflecting wall when closed) on i edges.
                const Box glo = amrex::adjCellLo(vbx, ix, nghost)
                    & amrex::adjCellLo(domain, ix, nghost);
                const Box ghi = amrex::adjCellHi(vbx, ix, nghost)
                    & amrex::adjCellHi(domain, ix, nghost);
                if (glo.ok()) {
                ParallelFor(glo,
                            [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    amrex::ignore_unused(k);
                    if (closed) {
                        const Real ex = v(i + 1, j + 1, 0, 0) - v(i + 1, j, 0, 0);
                        const Real ey = v(i + 1, j + 1, 0, 1) - v(i + 1, j, 0, 1);
                        Real rg, mxg, myg;
                        wall_ghost(ex, ey, u(i + 1, j, k, URHO),
                                   u(i + 1, j, k, UMX), u(i + 1, j, k, UMY),
                                   rg, mxg, myg);
                        u(i, j, k, URHO) = rg;
                        u(i, j, k, UMX)  = mxg;
                        u(i, j, k, UMY)  = myg;
                    } else {
                        for (int n = 0; n < ncomp; ++n) {
                            u(i, j, k, n) = u0(i + 1, j, k, n);
                        }
                    }
                });
                }
                if (ghi.ok()) {
                ParallelFor(ghi,
                            [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    amrex::ignore_unused(k);
                    if (closed) {
                        const Real ex = v(i, j + 1, 0, 0) - v(i, j, 0, 0);
                        const Real ey = v(i, j + 1, 0, 1) - v(i, j, 0, 1);
                        Real rg, mxg, myg;
                        wall_ghost(ex, ey, u(i - 1, j, k, URHO),
                                   u(i - 1, j, k, UMX), u(i - 1, j, k, UMY),
                                   rg, mxg, myg);
                        u(i, j, k, URHO) = rg;
                        u(i, j, k, UMX)  = mxg;
                        u(i, j, k, UMY)  = myg;
                    } else {
                        for (int n = 0; n < ncomp; ++n) {
                            u(i, j, k, n) = u0(i - 1, j, k, n);
                        }
                    }
                });
                }
            }
        }
    }

    // One explicit Euler step: Rusanov fluxes through the four curvilinear
    // faces plus the Roche-potential and Coriolis sources.
    void Advance(Real dt) {
        const RocheParams p = params;
        for (MFIter mfi(U, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Array4<Real const> u = U.const_array(mfi);
            Array4<Real const> u0 = Uinit.const_array(mfi);
            Array4<Real> un = Unew.array(mfi);
            Array4<Real const> v = xyv.const_array(mfi);
            Array4<Real const> vol_a = vol.const_array(mfi);
            ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::ignore_unused(k);
                const Real x00 = v(i, j, 0, 0),     y00 = v(i, j, 0, 1);
                const Real x10 = v(i + 1, j, 0, 0), y10 = v(i + 1, j, 0, 1);
                const Real x11 = v(i + 1, j + 1, 0, 0), y11 = v(i + 1, j + 1, 0, 1);
                const Real x01 = v(i, j + 1, 0, 0), y01 = v(i, j + 1, 0, 1);
                // Outward area vectors. The (i, j) -> (x, y) map is
                // left-handed, so with edges taken in v00->v10->v11->v01
                // order the outward normal is A = (-ey, ex).
                const Real Ailo[2] = {y01 - y00, x00 - x01};
                const Real Ajlo[2] = {y00 - y10, x10 - x00};
                const Real Aihi[2] = {y10 - y11, x11 - x10};
                const Real Ajhi[2] = {y11 - y01, x01 - x11};

                Real UL[ncomp];
                Real UR[ncomp];
                Real ZL[ncomp];
                Real ZR[ncomp];
                Real F[ncomp];
                Real fsum[ncomp] = {0.0_rt, 0.0_rt, 0.0_rt};
                for (int n = 0; n < ncomp; ++n) { UL[n] = u(i, j, k, n); }
                for (int n = 0; n < ncomp; ++n) { ZL[n] = u0(i, j, k, n); }

                for (int n = 0; n < ncomp; ++n) { UR[n] = u(i - 1, j, k, n); }
                for (int n = 0; n < ncomp; ++n) { ZR[n] = u0(i - 1, j, k, n); }
                rusanov_flux(UL, UR, ZL, ZR, Ailo[0], Ailo[1], p.cs, F);
                for (int n = 0; n < ncomp; ++n) { fsum[n] += F[n]; }

                for (int n = 0; n < ncomp; ++n) { UR[n] = u(i + 1, j, k, n); }
                for (int n = 0; n < ncomp; ++n) { ZR[n] = u0(i + 1, j, k, n); }
                rusanov_flux(UL, UR, ZL, ZR, Aihi[0], Aihi[1], p.cs, F);
                for (int n = 0; n < ncomp; ++n) { fsum[n] += F[n]; }

                for (int n = 0; n < ncomp; ++n) { UR[n] = u(i, j - 1, k, n); }
                for (int n = 0; n < ncomp; ++n) { ZR[n] = u0(i, j - 1, k, n); }
                rusanov_flux(UL, UR, ZL, ZR, Ajlo[0], Ajlo[1], p.cs, F);
                for (int n = 0; n < ncomp; ++n) { fsum[n] += F[n]; }

                for (int n = 0; n < ncomp; ++n) { UR[n] = u(i, j + 1, k, n); }
                for (int n = 0; n < ncomp; ++n) { ZR[n] = u0(i, j + 1, k, n); }
                rusanov_flux(UL, UR, ZL, ZR, Ajhi[0], Ajhi[1], p.cs, F);
                for (int n = 0; n < ncomp; ++n) { fsum[n] += F[n]; }

                const Real dtv = dt / vol_a(i, j, k);
                Real rho_new = UL[URHO] - dtv * fsum[URHO];
                Real mx_new  = UL[UMX]  - dtv * fsum[UMX];
                Real my_new  = UL[UMY]  - dtv * fsum[UMY];
                // Density floor: keep velocity if possible, else stop the gas.
                if (rho_new < p.rho_floor) {
                    const Real rescale = (rho_new > 0.0_rt)
                        ? p.rho_floor / rho_new : 0.0_rt;
                    rho_new = p.rho_floor;
                    mx_new *= rescale;
                    my_new *= rescale;
                }
                // Binary gravity + centrifugal at the cell center.
                const Real xc = 0.25_rt * (x00 + x10 + x11 + x01);
                const Real yc = 0.25_rt * (y00 + y10 + y11 + y01);
                const Real dx1 = xc - p.x1, dx2 = xc - p.x2;
                const Real r1sq = dx1 * dx1 + yc * yc + p.soft * p.soft;
                const Real r2sq = dx2 * dx2 + yc * yc + p.soft * p.soft;
                const Real gx = -p.m1 * dx1 / (r1sq * std::sqrt(r1sq))
                              - p.m2 * dx2 / (r2sq * std::sqrt(r2sq))
                              + p.omega * p.omega * xc;
                const Real gy = -p.m1 * yc / (r1sq * std::sqrt(r1sq))
                              - p.m2 * yc / (r2sq * std::sqrt(r2sq))
                              + p.omega * p.omega * yc;
                mx_new += dt * rho_new * gx;
                my_new += dt * rho_new * gy;
                // Coriolis: -2 Omega z-hat x v.
                const Real mx_c = mx_new, my_c = my_new;
                mx_new += dt * 2.0_rt * p.omega * my_c;
                my_new -= dt * 2.0_rt * p.omega * mx_c;

                // Velocity cap (roche.vmax > 0): limit |v| after all sources.
                // Expansion into near-vacuum (e.g. the semi-detached IC's L1
                // discontinuity: donor at rho_l1 next to the floor) gives
                // floor-density cells a pressure kick whose v = mom/rho blows
                // up since rho ~ 0; that collapses dt and stalls the run. The
                // cap clips those jets while preserving direction. Off by
                // default, so the overcontact mode is bitwise unchanged.
                if (p.vmax > 0.0_rt) {
                    const Real vmag = std::sqrt(mx_new * mx_new
                                                + my_new * my_new) / rho_new;
                    if (vmag > p.vmax) {
                        const Real fac = p.vmax / vmag;
                        mx_new *= fac;
                        my_new *= fac;
                    }
                }

                un(i, j, k, URHO) = rho_new;
                un(i, j, k, UMX)  = mx_new;
                un(i, j, k, UMY)  = my_new;
            });
        }
        std::swap(U, Unew);
    }

    // Ghost cells shared between grids of this block (not needed at the
    // block's outer edges -- those are covered by seams and physical BCs).
    void FillGhosts() { U.FillBoundary(Geom(0).periodicity()); }

    // Ghost cells of the reference state must be filled with the SAME rules
    // as U's ghosts (seam copies, wall mirrors, reservoir copies) -- anything
    // else (e.g. the analytic profile inside the star hole) makes the
    // deviation U - U0 nonzero at boundary faces even in equilibrium, and the
    // well-balanced dissipation turns into a mass pump. Call after the fill
    // pipeline has run once on U at t = 0.
    void SnapshotInitGhosts() {
        MultiFab::Copy(Uinit, U, 0, 0, ncomp, nghost);
        if (params.counter_rotate) {
            // The well-balanced reference is the hydrostatic atmosphere AT
            // REST, not the counter-rotating IC. Zero the reference momentum
            // so the Rusanov dissipation acts on the full momentum field:
            // the counter-rotating shear v = omega*r is strong, and a
            // dissipation that vanishes on it is a central scheme that goes
            // unstable (velocity runaway to 1e9 in <500 steps). The density
            // reference is unchanged, so the steep atmosphere is still
            // protected. Standard-Rusanov momentum diffusion then brakes the
            // counter-rotation over ~8 orbits (R^2/(smax*dx) ~ 290).
            Uinit.setVal(0.0_rt, UMX, 2, nghost);
        }
    }

    // CFL-limited timestep for this block: V / sum_faces (|u.n| + cs) |A|.
    Real ComputeDt() const {
        const RocheParams p = params;
        ReduceOps<ReduceOpMin> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(U); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            Array4<Real const> u = U.const_array(mfi);
            Array4<Real const> v = xyv.const_array(mfi);
            Array4<Real const> vol_a = vol.const_array(mfi);
            reduce_op.eval(vbx, reduce_data,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                amrex::ignore_unused(k);
                const Real x00 = v(i, j, 0, 0),     y00 = v(i, j, 0, 1);
                const Real x10 = v(i + 1, j, 0, 0), y10 = v(i + 1, j, 0, 1);
                const Real x11 = v(i + 1, j + 1, 0, 0), y11 = v(i + 1, j + 1, 0, 1);
                const Real x01 = v(i, j + 1, 0, 0), y01 = v(i, j + 1, 0, 1);
                const Real Ailo[2] = {y01 - y00, x00 - x01};
                const Real Ajlo[2] = {y00 - y10, x10 - x00};
                const Real Aihi[2] = {y10 - y11, x11 - x10};
                const Real Ajhi[2] = {y11 - y01, x01 - x11};
                const Real rho = u(i, j, k, URHO);
                const Real ux = u(i, j, k, UMX) / rho;
                const Real uy = u(i, j, k, UMY) / rho;
                Real lam = 0.0_rt;
                auto acc = [&](Real ax, Real ay) {
                    lam += std::abs(ux * ax + uy * ay)
                        + p.cs * std::sqrt(ax * ax + ay * ay);
                };
                acc(Ailo[0], Ailo[1]);
                acc(Aihi[0], Aihi[1]);
                acc(Ajlo[0], Ajlo[1]);
                acc(Ajhi[0], Ajhi[1]);
                return {vol_a(i, j, k) / lam};
            });
        }
        Real dt = amrex::get<0>(reduce_data.value());
        ParallelDescriptor::ReduceRealMin(dt);
        return p.cfl * dt;
    }

    // Min and max density in this block (NaN detector).
    std::pair<Real, Real> RhoMinMax() const {
        ReduceOps<ReduceOpMin, ReduceOpMax> reduce_op;
        ReduceData<Real, Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(U); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            Array4<Real const> u = U.const_array(mfi);
            reduce_op.eval(vbx, reduce_data,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                return {u(i, j, k, URHO), u(i, j, k, URHO)};
            });
        }
        auto hv = reduce_data.value();
        ParallelDescriptor::ReduceRealMin(amrex::get<0>(hv));
        ParallelDescriptor::ReduceRealMax(amrex::get<1>(hv));
        return {amrex::get<0>(hv), amrex::get<1>(hv)};
    }

    // Largest gas speed in this block (hydrostatic-hold metric).
    Real MaxVel() const {
        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(U); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            Array4<Real const> u = U.const_array(mfi);
            reduce_op.eval(vbx, reduce_data,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                const Real rho = u(i, j, k, URHO);
                const Real mx = u(i, j, k, UMX);
                const Real my = u(i, j, k, UMY);
                return {std::sqrt(mx * mx + my * my) / rho};
            });
        }
        Real v = amrex::get<0>(reduce_data.value());
        ParallelDescriptor::ReduceRealMax(v);
        return v;
    }
    // Mass-weighted speed moments for this block: (sum rho*|u|*vol,
    // sum rho*vol). The ratio, reduced over all blocks, is the bulk flow
    // speed -- the honest "is the gas corotating?" metric. max|v| is
    // dominated by a few tenuous jet cells; the bulk speed tracks the mean.
    std::pair<Real, Real> BulkSpeedMoments() const {
        ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
        ReduceData<Real, Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(U); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            Array4<Real const> u = U.const_array(mfi);
            Array4<Real const> vol_a = vol.const_array(mfi);
            reduce_op.eval(vbx, reduce_data,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                const Real rho = u(i, j, k, URHO);
                const Real mx = u(i, j, k, UMX);
                const Real my = u(i, j, k, UMY);
                const Real vmag = std::sqrt(mx * mx + my * my) / rho;
                return {rho * vmag * vol_a(i, j, k), rho * vol_a(i, j, k)};
            });
        }
        const auto r = reduce_data.value();
        Real num = amrex::get<0>(r);
        Real den = amrex::get<1>(r);
        ParallelDescriptor::ReduceRealSum(num);
        ParallelDescriptor::ReduceRealSum(den);
        return {num, den};
    }

    // Number of NaN density cells in this block (MPI_SUM-safe, unlike
    // min/max reductions which can mask NaNs).
    Long NaNCount() const {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(U); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            Array4<Real const> u = U.const_array(mfi);
            reduce_op.eval(vbx, reduce_data,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                return {std::isnan(u(i, j, k, URHO)) ? 1 : 0};
            });
        }
        Long count = amrex::get<0>(reduce_data.value());
        ParallelDescriptor::ReduceLongSum(count);
        return count;
    }

    // Total mass in this block: sum of rho * volume.
    Real TotalMass() const {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(U); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            Array4<Real const> u = U.const_array(mfi);
            Array4<Real const> vol_a = vol.const_array(mfi);
            reduce_op.eval(vbx, reduce_data,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                return {u(i, j, k, URHO) * vol_a(i, j, k)};
            });
        }
        Real mass = amrex::get<0>(reduce_data.value());
        ParallelDescriptor::ReduceRealSum(mass);
        return mass;
    }

    // Total cell volume in this block (geometry check).
    Real TotalVolume() const {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(U); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            Array4<Real const> vol_a = vol.const_array(mfi);
            reduce_op.eval(vbx, reduce_data,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                return {vol_a(i, j, k)};
            });
        }
        Real v = amrex::get<0>(reduce_data.value());
        ParallelDescriptor::ReduceRealSum(v);
        return v;
    }

    // Vertex bounding box of this block (geometry check).
    void VertexBBox(Real* lo, Real* hi) const {
        lo[0] = lo[1] = std::numeric_limits<Real>::max();
        hi[0] = hi[1] = std::numeric_limits<Real>::lowest();
        for (MFIter mfi(xyv); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();  // xyv is nodal: this IS the vertex box
            Array4<Real const> v = xyv.const_array(mfi);
            for (int j = bx.smallEnd(iy); j <= bx.bigEnd(iy); ++j) {
                for (int i = bx.smallEnd(ix); i <= bx.bigEnd(ix); ++i) {
                    lo[0] = amrex::min(lo[0], v(i, j, 0, 0));
                    lo[1] = amrex::min(lo[1], v(i, j, 0, 1));
                    hi[0] = amrex::max(hi[0], v(i, j, 0, 0));
                    hi[1] = amrex::max(hi[1], v(i, j, 0, 1));
                }
            }
        }
        ParallelDescriptor::ReduceRealMin(lo[0]);
        ParallelDescriptor::ReduceRealMin(lo[1]);
        ParallelDescriptor::ReduceRealMax(hi[0]);
        ParallelDescriptor::ReduceRealMax(hi[1]);
    }

    // Mass flux through the bridge midplane (the j = n_bridge/2 plane of the
    // bridge block), positive toward star 2. Only meaningful on the bridge.
    Real MidplaneFlux() const {
        const int jm = params.n_bridge / 2 - 1;  // cell whose j-hi face is mid
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(U); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            if (jm < vbx.smallEnd(iy) || jm > vbx.bigEnd(iy)) continue;
            const Box row = vbx & Box(IntVect{AMREX_D_DECL(0, jm, 0)},
                                      IntVect{AMREX_D_DECL(params.n_phi - 1, jm, 0)});
            Array4<Real const> u = U.const_array(mfi);
            Array4<Real const> v = xyv.const_array(mfi);
            reduce_op.eval(row, reduce_data,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                amrex::ignore_unused(k, j);
                // Outward j-hi face area of cell (i, jm).
                const Real ax = v(i + 1, j + 1, 0, 1) - v(i, j + 1, 0, 1);
                const Real ay = v(i, j + 1, 0, 0) - v(i + 1, j + 1, 0, 0);
                const Real mx = 0.5_rt * (u(i, j, k, UMX) + u(i, j + 1, k, UMX));
                const Real my = 0.5_rt * (u(i, j, k, UMY) + u(i, j + 1, k, UMY));
                return {mx * ax + my * ay};
            });
        }
        Real f = amrex::get<0>(reduce_data.value());
        ParallelDescriptor::ReduceRealSum(f);
        return f;
    }

    // Poison test: set every ghost cell of U to NaN. After the fill pipeline
    // (FillGhosts + seams + physical BCs), UnfilledGhostCount() must be 0.
    void PoisonGhosts() {
        for (MFIter mfi(U); mfi.isValid(); ++mfi) {
            const Box vbx = mfi.validbox();
            Array4<Real> u = U.array(mfi);
            ParallelFor(amrex::grow(vbx, nghost),
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (!vbx.contains(i, j, k)) {
                    u(i, j, k, URHO) = std::numeric_limits<Real>::quiet_NaN();
                    u(i, j, k, UMX)  = std::numeric_limits<Real>::quiet_NaN();
                    u(i, j, k, UMY)  = std::numeric_limits<Real>::quiet_NaN();
                }
            });
        }
    }

    // Number of NaN cells in the four face-adjacent ghost slabs of U.
    // (Corner ghost cells are intentionally excluded: nothing fills them and
    // the first-order stencil never reads them.)
    Long UnfilledGhostCount() const {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(U); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            Array4<Real const> u = U.const_array(mfi);
            const Box slabs[4] = {amrex::adjCellLo(vbx, ix, nghost),
                                  amrex::adjCellHi(vbx, ix, nghost),
                                  amrex::adjCellLo(vbx, iy, nghost),
                                  amrex::adjCellHi(vbx, iy, nghost)};
            for (const Box& slab : slabs) {
                reduce_op.eval(slab, reduce_data,
                               [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
                    const bool bad = std::isnan(u(i, j, k, URHO))
                                  || std::isnan(u(i, j, k, UMX))
                                  || std::isnan(u(i, j, k, UMY));
                    return {bad ? 1 : 0};
                });
            }
        }
        Long count = amrex::get<0>(reduce_data.value());
        ParallelDescriptor::ReduceLongSum(count);
        return count;
    }

    MultiFab U{};      // (rho, mom_x, mom_y), cell-centered, one ghost layer
    MultiFab Unew{};
    MultiFab Uinit{};  // initial state: the open-boundary reservoir
    MultiFab xyv{};    // physical vertex coordinates, nodal, 2 components
    MultiFab vol{};    // cell volumes
    RocheParams params{};
    BlockGeom g{};

  private:
    void ErrorEst(int, ::amrex::TagBoxArray&, Real, int) override {
        throw std::runtime_error("single-level only");
    }
    void MakeNewLevelFromScratch(int level, Real, const ::amrex::BoxArray& ba,
                                 const ::amrex::DistributionMapping& dm) override {
        if (level > 0) throw std::runtime_error("single-level only");
        U.define(ba, dm, ncomp, nghost);
        Unew.define(ba, dm, ncomp, nghost);
        Uinit.define(ba, dm, ncomp, nghost);
        xyv.define(amrex::convert(ba, IntVect{AMREX_D_DECL(1, 1, 1)}), dm, 2, 0);
        vol.define(ba, dm, 1, 0);
    }
    void MakeNewLevelFromCoarse(int, Real, const ::amrex::BoxArray&,
                                const ::amrex::DistributionMapping&) override {
        throw std::runtime_error("single-level only");
    }
    void RemakeLevel(int, Real, const ::amrex::BoxArray&,
                     const ::amrex::DistributionMapping&) override {
        throw std::runtime_error("single-level only");
    }
    void ClearLevel(int level) override {
        if (level > 0) throw std::runtime_error("single-level only");
        U.clear();
        Unew.clear();
        Uinit.clear();
        xyv.clear();
        vol.clear();
    }
};

using namespace NonLocalBC;

// One-sided fill of U from src -> dest across a seam. The mapping may be
// offset-only (aligned edges) or include a sign flip (reversed edges).
struct SeamBoundaryFn {
    RocheBlock* dest;
    const RocheBlock* src;
    MultiBlockIndexMapping dtos;
    Box boundary_to_fill;
    std::unique_ptr<MultiBlockCommMetaData> cmd{};
    FabArrayBase::BDKey cached_dest_bd_key{};
    FabArrayBase::BDKey cached_src_bd_key{};
    ApplyDtosAndProjectionOnReciever<MultiBlockIndexMapping, MapComponents<Identity>>
        packing{PackComponents{0, 0, ncomp}, dtos};

    AMREX_NODISCARD CommHandler FillBoundary_nowait() {
        if (!cmd || cached_dest_bd_key != dest->U.getBDKey()
            || cached_src_bd_key != src->U.getBDKey()) {
            cmd = std::make_unique<MultiBlockCommMetaData>(
                dest->U, boundary_to_fill, src->U, dest->U.nGrowVect(), dtos);
            cached_dest_bd_key = dest->U.getBDKey();
            cached_src_bd_key = src->U.getBDKey();
        }
        return ParallelCopy_nowait(no_local_copy, dest->U, src->U, *cmd, packing);
    }
    void FillBoundary_do_local_copy() const {
        if (cmd->m_LocTags && !cmd->m_LocTags->empty()) {
            LocalCopy(packing, dest->U, src->U, *cmd->m_LocTags);
        }
    }
    void FillBoundary_finish(CommHandler handler) const {
        ParallelCopy_finish(dest->U, std::move(handler), *cmd, packing); // NOLINT
    }
};

struct FillBoundaryFn {
    std::vector<SeamBoundaryFn> boundaries;
    void operator()() {
        std::vector<CommHandler> comms;
        comms.reserve(boundaries.size());
        for (auto& b : boundaries) {
            comms.emplace_back(b.FillBoundary_nowait());
        }
        for (auto& b : boundaries) b.FillBoundary_do_local_copy();
        for (std::size_t i = 0; i < boundaries.size(); ++i) {
            boundaries[i].FillBoundary_finish(std::move(comms[i])); // NOLINT
        }
    }
};

void WritePlotfile(const RocheBlock& block, const char* dir, Real t, int step) {
    static const Vector<std::string> varnames{"rho", "momx", "momy", "x", "y"};
    char buf[256];
    snprintf(buf, sizeof(buf), "RocheBinary/%s/plt%04d", dir, step);
    MultiFab plotmf(block.U.boxArray(), block.U.DistributionMap(), 5, 0);
    MultiFab::Copy(plotmf, block.U, 0, 0, ncomp, 0);
    for (MFIter mfi(plotmf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real> pf = plotmf.array(mfi);
        Array4<Real const> v = block.xyv.const_array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            amrex::ignore_unused(k);
            pf(i, j, k, 3) = 0.25_rt * (v(i, j, 0, 0) + v(i + 1, j, 0, 0)
                                      + v(i + 1, j + 1, 0, 0) + v(i, j + 1, 0, 0));
            pf(i, j, k, 4) = 0.25_rt * (v(i, j, 0, 1) + v(i + 1, j, 0, 1)
                                      + v(i + 1, j + 1, 0, 1) + v(i, j + 1, 0, 1));
        });
    }
    Vector<const MultiFab*> mf{&plotmf};
    Vector<Geometry> geoms{block.Geom(0)};
    Vector<int> level_steps{step};
    Vector<IntVect> ref_ratio{};
    WriteMultiLevelPlotfile(buf, 1, mf, varnames, geoms, t, level_steps, ref_ratio);
}

// Tangential (i-edge) seam ghost slabs: the i-low (i = -1) and i-high
// (i = n_phi) ghost columns of a block's domain.
Box seam_ghost_lo(const Box& domain) {
    return Box{IntVect{AMREX_D_DECL(domain.smallEnd(ix) - 1, domain.smallEnd(iy), 0)},
               IntVect{AMREX_D_DECL(domain.smallEnd(ix) - 1, domain.bigEnd(iy), 0)}};
}
Box seam_ghost_hi(const Box& domain) {
    return Box{IntVect{AMREX_D_DECL(domain.bigEnd(ix) + 1, domain.smallEnd(iy), 0)},
               IntVect{AMREX_D_DECL(domain.bigEnd(ix) + 1, domain.bigEnd(iy), 0)}};
}

// Radial/corridor (j-edge) seam ghost slabs: the j-low (j = -1) and j-high
// (j = n_j) ghost rows of a block's domain.
Box seam_ghost_jlo(const Box& domain) {
    return Box{IntVect{AMREX_D_DECL(domain.smallEnd(ix), domain.smallEnd(iy) - 1, 0)},
               IntVect{AMREX_D_DECL(domain.bigEnd(ix), domain.smallEnd(iy) - 1, 0)}};
}
Box seam_ghost_jhi(const Box& domain) {
    return Box{IntVect{AMREX_D_DECL(domain.smallEnd(ix), domain.bigEnd(iy) + 1, 0)},
               IntVect{AMREX_D_DECL(domain.bigEnd(ix), domain.bigEnd(iy) + 1, 0)}};
}

} // namespace

void MyMain() {
    RocheParams p;
    Real stop_time = 1.0_rt;  // orbital period is 2 pi/omega ~ 35.5
    int max_steps = 1000000;
    int plot_int = 100;
    int print_int = 10;
    int poison_test = 0;
    {
        ParmParse pp("roche");
        pp.query("m1", p.m1);
        pp.query("m2", p.m2);
        pp.query("sep", p.sep);
        pp.query("a", p.a);
        pp.query("r_ann1", p.r_ann1);
        pp.query("r_ann2", p.r_ann2);
        pp.query("stretch", p.stretch);
        pp.query("cs", p.cs);
        pp.query("rho_l1", p.rho_l1);
        pp.query("rho_l1_1", p.rho_l1_1);
        pp.query("rho_l1_2", p.rho_l1_2);
        pp.query("neck_blend", p.neck_blend);
        pp.query("rho_floor", p.rho_floor);
        pp.query("soft", p.soft);
        pp.query("cfl", p.cfl);
        pp.query("vmax", p.vmax);
        pp.query("perturb", p.perturb);
        pp.query("n_phi", p.n_phi);
        pp.query("n_r", p.n_r);
        pp.query("n_bridge", p.n_bridge);
        pp.query("stop_time", stop_time);
        pp.query("max_steps", max_steps);
        pp.query("plot_int", plot_int);
        pp.query("print_int", print_int);
        pp.query("poison_test", poison_test);
        std::string ic_mode_str = "overcontact";
        pp.query("ic_mode", ic_mode_str);
        if (ic_mode_str == "semi_detached" || ic_mode_str == "semi-detached") {
            p.ic_mode = IC_SEMI_DETACHED;
        } else if (ic_mode_str == "overcontact") {
            p.ic_mode = IC_OVERCONTACT;
        } else {
            amrex::Abort("roche.ic_mode must be 'overcontact' or 'semi_detached'");
        }
        std::string outer_bc_str = "reservoir";
        pp.query("outer_bc", outer_bc_str);
        if (outer_bc_str == "closed") {
            p.outer_closed = 1;
        } else if (outer_bc_str == "reservoir") {
            p.outer_closed = 0;
        } else {
            amrex::Abort("roche.outer_bc must be 'reservoir' or 'closed'");
        }
        pp.query("counter_rotate", p.counter_rotate);
    }
    derive(p);
    amrex::Print().SetPrecision(10)
        << "RocheBinary: omega = " << p.omega
        << ", period = " << 2.0_rt * M_PI / p.omega
        << ", stars at " << p.x1 << ", " << p.x2
        << ", L1 at " << p.x_l1 << ", Phi_L1 = " << p.phi_l1 << '\n';
    amrex::Print() << "IC mode: "
                   << (p.ic_mode == IC_SEMI_DETACHED ? "semi_detached"
                                                     : "overcontact")
                   << '\n';
    if (p.rho_l1_1 != p.rho_l1_2) {
        amrex::Print().SetPrecision(10)
            << "Differential contact depth: rho_l1_1 = " << p.rho_l1_1
            << ", rho_l1_2 = " << p.rho_l1_2
            << ", neck blend width = " << p.blend_width
            << " (|Phi_xx(L1)| = " << std::abs(p.phi_xx_l1) << ")\n";
    }
    if (p.outer_closed) {
        amrex::Print() << "Outer boundary: closed (reflecting walls)\n";
    }
    if (p.counter_rotate) {
        amrex::Print() << "IC: gas at inertial rest -- counter-rotating at "
                          "omega*(y,-x); use outer_bc=closed\n";
    }

    // Polar blocks share one logical domain box; the bridge has its own.
    Box domain_polar(IntVect{}, IntVect{AMREX_D_DECL(p.n_phi - 1, p.n_r - 1, 0)});
    Box domain_bridge(IntVect{}, IntVect{AMREX_D_DECL(p.n_phi - 1, p.n_bridge - 1, 0)});
    const int N = p.n_phi;
    Array<int, AMREX_SPACEDIM> per_none{AMREX_D_DECL(0, 0, 0)};

    auto make_geom = [&](const BlockGeom& bg, const Box& domain) {
        std::array<Real, AMREX_SPACEDIM> lo{AMREX_D_DECL(std::numeric_limits<Real>::max(),
                                                         std::numeric_limits<Real>::max(),
                                                         std::numeric_limits<Real>::max())};
        std::array<Real, AMREX_SPACEDIM> hi{AMREX_D_DECL(std::numeric_limits<Real>::lowest(),
                                                         std::numeric_limits<Real>::lowest(),
                                                         std::numeric_limits<Real>::lowest())};
        for (int j = 0; j <= domain.bigEnd(iy) + 1; ++j) {
            for (int i = 0; i <= domain.bigEnd(ix) + 1; ++i) {
                const auto pos = vertex_pos(i, j, bg, p);
                lo[0] = amrex::min(lo[0], pos[0]);
                lo[1] = amrex::min(lo[1], pos[1]);
                hi[0] = amrex::max(hi[0], pos[0]);
                hi[1] = amrex::max(hi[1], pos[1]);
            }
        }
        return Geometry{domain, RealBox{lo, hi}, CoordSys::cartesian, per_none};
    };

    AmrInfo amr_info{};

    // Star 1 (at x1) and star 2 (at x2), each with E, N, W, S blocks.
    const BlockGeom g_e1{BlockKind::Polar, p.x1, 0.0_rt, 0.0_rt, p.r_ann1, /*jhi_seam=*/true};
    const BlockGeom g_n1{BlockKind::Polar, p.x1, 0.0_rt, 0.5_rt * M_PI, p.r_ann1};
    const BlockGeom g_w1{BlockKind::Polar, p.x1, 0.0_rt, 1.0_rt * M_PI, p.r_ann1};
    const BlockGeom g_s1{BlockKind::Polar, p.x1, 0.0_rt, 1.5_rt * M_PI, p.r_ann1};
    const BlockGeom g_e2{BlockKind::Polar, p.x2, 0.0_rt, 0.0_rt, p.r_ann2};
    const BlockGeom g_n2{BlockKind::Polar, p.x2, 0.0_rt, 0.5_rt * M_PI, p.r_ann2};
    const BlockGeom g_w2{BlockKind::Polar, p.x2, 0.0_rt, 1.0_rt * M_PI, p.r_ann2, /*jhi_seam=*/true};
    const BlockGeom g_s2{BlockKind::Polar, p.x2, 0.0_rt, 1.5_rt * M_PI, p.r_ann2};
    const BlockGeom g_br{BlockKind::Bridge, 0.0_rt, 0.0_rt, 0.0_rt};

    RocheBlock e1(g_e1, make_geom(g_e1, domain_polar), p, amr_info);
    RocheBlock n1(g_n1, make_geom(g_n1, domain_polar), p, amr_info);
    RocheBlock w1(g_w1, make_geom(g_w1, domain_polar), p, amr_info);
    RocheBlock s1(g_s1, make_geom(g_s1, domain_polar), p, amr_info);
    RocheBlock e2(g_e2, make_geom(g_e2, domain_polar), p, amr_info);
    RocheBlock n2(g_n2, make_geom(g_n2, domain_polar), p, amr_info);
    RocheBlock w2(g_w2, make_geom(g_w2, domain_polar), p, amr_info);
    RocheBlock s2(g_s2, make_geom(g_s2, domain_polar), p, amr_info);
    RocheBlock br(g_br, make_geom(g_br, domain_bridge), p, amr_info);

    // Geometry report: per-block bounding box and volume. Each polar block
    // is a quarter annulus of area pi (r_annK^2 - a^2) / 4.
    {
        const Real vol_an1 = M_PI * (p.r_ann1 * p.r_ann1 - p.a * p.a) / 4.0_rt;
        const Real vol_an2 = M_PI * (p.r_ann2 * p.r_ann2 - p.a * p.a) / 4.0_rt;
        amrex::Print().SetPrecision(10)
            << "Polar block volume (analytic): " << vol_an1
            << " (star 1), " << vol_an2 << " (star 2)\n";
        RocheBlock* all[9] = {&e1, &n1, &w1, &s1, &e2, &n2, &w2, &s2, &br};
        const char* names[9] = {"e1", "n1", "w1", "s1", "e2", "n2", "w2", "s2", "br"};
        for (int b = 0; b < 9; ++b) {
            Real lo[2], hi[2];
            all[b]->VertexBBox(lo, hi);
            amrex::Print().SetPrecision(10)
                << "  " << names[b] << ": bbox x in [" << lo[0] << ", " << hi[0]
                << "], y in [" << lo[1] << ", " << hi[1]
                << "], volume = " << all[b]->TotalVolume() << '\n';
        }
    }

    // Aligned-block index mapping: src_idx[d] = dest_idx[d] - offset[d].
    auto make_dtos = [](IntVect off) {
        MultiBlockIndexMapping d{};
        d.permutation = IntVect{AMREX_D_DECL(0, 1, 2)};
        d.sign        = IntVect{AMREX_D_DECL(1, 1, 1)};
        d.offset      = off;
        return d;
    };
    // i-reversed mapping: src[0] = -(dest[0] - off[0]), src[1] = dest[1] - off[1].
    auto make_dtos_flip = [](IntVect off) {
        MultiBlockIndexMapping d{};
        d.permutation = IntVect{AMREX_D_DECL(0, 1, 2)};
        d.sign        = IntVect{AMREX_D_DECL(-1, 1, 1)};
        d.offset      = off;
        return d;
    };

    std::vector<SeamBoundaryFn> bnd;
    // Tangential seams close each annulus: E -> N -> W -> S -> E (offset-only,
    // identical to StarDisk). Ghost column i = N maps to interior i = 0 of
    // the next block (offset = (N, 0)); ghost i = -1 to interior i = N-1 of
    // the previous (offset = (-N, 0)).
    auto add_annulus_seams = [&](RocheBlock* be, RocheBlock* bn, RocheBlock* bw,
                                 RocheBlock* bs) {
        bnd.push_back({be, bn, make_dtos(N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                       seam_ghost_hi(domain_polar)});
        bnd.push_back({bn, be, make_dtos(-N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                       seam_ghost_lo(domain_polar)});
        bnd.push_back({bn, bw, make_dtos(N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                       seam_ghost_hi(domain_polar)});
        bnd.push_back({bw, bn, make_dtos(-N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                       seam_ghost_lo(domain_polar)});
        bnd.push_back({bw, bs, make_dtos(N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                       seam_ghost_hi(domain_polar)});
        bnd.push_back({bs, bw, make_dtos(-N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                       seam_ghost_lo(domain_polar)});
        bnd.push_back({bs, be, make_dtos(N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                       seam_ghost_hi(domain_polar)});
        bnd.push_back({be, bs, make_dtos(-N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                       seam_ghost_lo(domain_polar)});
    };
    add_annulus_seams(&e1, &n1, &w1, &s1);
    add_annulus_seams(&e2, &n2, &w2, &s2);

    // Bridge seams.
    //   bridge j-low <-> E1 j-high: offset-only, both i run bottom -> top.
    //     br ghost (i, -1)  <- e1 interior (i, n_r-1):  offset (0, -n_r)
    //     e1 ghost (i, n_r) <- br interior (i, 0):      offset (0,  n_r)
    bnd.push_back({&br, &e1, make_dtos(IntVect{AMREX_D_DECL(0, -p.n_r, 0)}),
                   seam_ghost_jlo(domain_bridge)});
    bnd.push_back({&e1, &br, make_dtos(IntVect{AMREX_D_DECL(0, p.n_r, 0)}),
                   seam_ghost_jhi(domain_polar)});
    //   bridge j-high <-> W2 j-high: W2's i runs top -> bottom, so the
    //   mapping reverses i: src = (n_phi-1-i, ...).
    //     br ghost (i, n_bridge)  <- w2 interior (n_phi-1-i, n_r-1)
    //     w2 ghost (i, n_r)       <- br interior (n_phi-1-i, n_bridge-1)
    bnd.push_back({&br, &w2,
                   make_dtos_flip(IntVect{AMREX_D_DECL(p.n_phi - 1,
                                                       p.n_bridge - p.n_r + 1, 0)}),
                   seam_ghost_jhi(domain_bridge)});
    bnd.push_back({&w2, &br,
                   make_dtos_flip(IntVect{AMREX_D_DECL(p.n_phi - 1,
                                                       p.n_r - p.n_bridge + 1, 0)}),
                   seam_ghost_jhi(domain_polar)});

    FillBoundaryFn FillSeams{std::move(bnd)};

    RocheBlock* blocks[9] = {&e1, &n1, &w1, &s1, &e2, &n2, &w2, &s2, &br};
    const char* names[9] = {"e1", "n1", "w1", "s1", "e2", "n2", "w2", "s2", "br"};

    int step = 0;
    Real t = 0.0_rt;
    auto region_mass = [&] {
        Real m1l = 0.0_rt, m2l = 0.0_rt;
        for (int b = 0; b < 4; ++b) m1l += blocks[b]->TotalMass();
        for (int b = 4; b < 8; ++b) m2l += blocks[b]->TotalMass();
        return std::array<Real, 3>{m1l, m2l, br.TotalMass()};
    };
    const Real mass0 = [&] {
        Real m = 0.0_rt;
        for (auto* b : blocks) m += b->TotalMass();
        return m;
    }();
    const auto rm0 = region_mass();
    amrex::Print().SetPrecision(10)
        << "Initial total mass: " << mass0
        << " (lobe1 " << rm0[0] << ", lobe2 " << rm0[1]
        << ", bridge " << rm0[2] << ")\n";

    for (int b = 0; b < 9; ++b) WritePlotfile(*blocks[b], names[b], t, step);

    // Establish U's t=0 ghost state with the production fill pipeline, then
    // snapshot interior + ghosts as the static reference U0. With U0's
    // ghosts filled by the same rules as U's, the deviation U - U0 is zero
    // at every face in equilibrium and the well-balanced dissipation
    // vanishes identically.
    for (auto* b : blocks) b->FillGhosts();
    FillSeams();
    for (auto* b : blocks) b->FillPhysicalBCs();
    for (auto* b : blocks) b->SnapshotInitGhosts();

    // Optional poison test: NaN-poison all ghost cells, run the full fill
    // pipeline once, and require every face-adjacent ghost cell to be filled.
    if (poison_test) {
        for (auto* b : blocks) b->PoisonGhosts();
        for (auto* b : blocks) b->FillGhosts();
        FillSeams();
        for (auto* b : blocks) b->FillPhysicalBCs();
        Long unfilled = 0;
        for (auto* b : blocks) unfilled += b->UnfilledGhostCount();
        amrex::Print() << "Poison test: unfilled face-ghost cells = "
                       << unfilled << '\n';
        if (unfilled > 0) {
            amrex::Abort("RocheBinary: poison test failed");
        }
    }
    // Time-series diagnostics for the L1-flux / lobe-equilibration figures
    // (parsed stdout is fragile; a structured file is reusable for setups
    // B/C/D). Written by the IO rank only; the reductions feeding it are
    // collective, so every rank holds the same values.
    std::ofstream diag;
    if (ParallelDescriptor::IOProcessor()) {
        diag.open("RocheBinary/diagnostics.dat", std::ios::out);
        diag << "# step t dt mass_drift m_lobe1 m_bridge m_lobe2 l1_flux max_v bulk_v\n";
    }
    auto write_diag = [&](int s, Real tt, Real dtt, Real mdrift, Real ml1,
                          Real mbr, Real ml2, Real flux, Real mxv, Real bvv) {
        if (ParallelDescriptor::IOProcessor()) {
            diag << s << ' ' << std::scientific << std::setprecision(12)
                 << tt << ' ' << dtt << ' ' << mdrift << ' ' << ml1 << ' '
                 << mbr << ' ' << ml2 << ' ' << flux << ' ' << mxv << ' '
                 << bvv << '\n';
            diag.flush();
        }
    };
    write_diag(0, 0.0_rt, 0.0_rt, 0.0_rt, rm0[0], rm0[2], rm0[1], 0.0_rt, 0.0_rt, 0.0_rt);

    while (t < stop_time && step < max_steps) {
        Real dt = blocks[0]->ComputeDt();
        for (int b = 1; b < 9; ++b) dt = amrex::min(dt, blocks[b]->ComputeDt());
        dt = amrex::min(dt, stop_time - t);

        // Intra-block fills first (upstream MultiBlock idiom): ghost cells
        // between fabs of the same block. Without this, a block split across
        // ranks reads stale inter-fab ghosts, and the NonLocalBC seam
        // machinery's ngrow-grown sliver overlaps read invalid src ghosts.
        for (auto* b : blocks) b->FillGhosts();
        FillSeams();
        for (auto* b : blocks) b->FillPhysicalBCs();
        for (auto* b : blocks) b->Advance(dt);
        t += dt;
        ++step;

        if (step % print_int == 0) {
            Real mass = 0.0_rt, maxv = 0.0_rt, bnum = 0.0_rt, bden = 0.0_rt;
            for (auto* b : blocks) {
                mass += b->TotalMass();
                maxv = amrex::max(maxv, b->MaxVel());
                const auto [n, d] = b->BulkSpeedMoments();
                bnum += n; bden += d;
            }
            const Real bulk_v = (bden > 0.0_rt) ? bnum / bden : 0.0_rt;
            const auto rm = region_mass();
            const Real l1_flux = br.MidplaneFlux();
            amrex::Print().SetPrecision(10)
                << "Step #" << step << ", t = " << t << ", dt = " << dt
                << ", mass drift = " << (mass - mass0) / mass0
                << ", m_lobe1 = " << rm[0] << ", m_bridge = " << rm[2]
                << ", m_lobe2 = " << rm[1]
                << ", L1 flux = " << l1_flux
                << ", max|v| = " << maxv
                << ", bulk|v| = " << bulk_v << '\n';
            write_diag(step, t, dt, (mass - mass0) / mass0,
                       rm[0], rm[2], rm[1], l1_flux, maxv, bulk_v);
            for (int b = 0; b < 9; ++b) {
                const auto [rlo, rhi] = blocks[b]->RhoMinMax();
                amrex::Print().SetPrecision(10)
                    << "    " << names[b] << ": mass = " << blocks[b]->TotalMass()
                    << ", rho in [" << rlo << ", " << rhi << "]"
                    << ", NaN cells = " << blocks[b]->NaNCount() << '\n';
            }
        }
        if (step % plot_int == 0 || t >= stop_time) {
            for (int b = 0; b < 9; ++b) WritePlotfile(*blocks[b], names[b], t, step);
        }
    }

    // A test should fail loudly: NaN anywhere means a broken boundary fill.
    Long nan_total = 0;
    for (auto* b : blocks) nan_total += b->NaNCount();
    if (nan_total > 0) {
        amrex::Abort("RocheBinary: NaN detected in final state");
    }
}
