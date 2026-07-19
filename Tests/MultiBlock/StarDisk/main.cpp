// 2D star + accretion disk example using AMReX NonLocalBC on a multi-block
// curvilinear grid.
//
// Layout: a fixed "star" occupies the central disk r < a (masked, not
// evolved). Four identical blocks surround it, each covering one quadrant
// of the circular annulus a < r < r_out:
//
//                 +-------+
//                 |   N   |
//         +-------+ .---. +-------+
//         |   W   |( star )   E   |
//         +-------+ '---' +-------+
//                 |   S   |
//                 +-------+
//
// Each block has logical coordinates (i = tangential, counterclockwise;
// j = radial, outward) and is mapped to physical space by
//
//   phi(i)    = -pi/4 + (i/n_phi)*(pi/2)         (angle relative to block)
//   rho(j)    = a + (r_out - a) * (j/n_r)^stretch
//   x = rho * cos(theta_b + phi),  y = rho * sin(theta_b + phi)
//
// with block orientation theta_b in {0, pi/2, pi, 3 pi/2} for E, N, W, S.
// Each block is a polar sector; the inner edge of the annulus is the circle
// r = a (the star surface).
//
// The four tangential seams (E<->N, N<->W, W<->S, S<->E) are glued with
// aligned (offset-only) MultiBlockIndexMapping fills, so the disk annulus is
// seamless all the way around -- no periodic boundary anywhere. The inner
// edge of each block is a reflecting wall (the star surface); the outer edge
// is zero-gradient outflow.
//
// Physics: isothermal Euler (rho, rho*u) with p = cs^2 * rho, advanced by
// first-order finite volume with Rusanov fluxes on the curvilinear mesh.
// Geometric factors (face area vectors, cell volumes) are computed from
// vertex coordinates, so no analytic Jacobians are needed. Gravity is a
// point mass GM at the origin. The initial disk is in radial equilibrium:
// rho ~ r^-q, v_phi = sqrt(GM/r - q*cs^2).
//
// Build: cd Tests/MultiBlock/StarDisk && make -j
// Run:   mpirun -n 2 ./main2d.gnu.MPI.ex [disk.name=value ...]
// Output: StarDisk/{e,n,w,s}/plt#### plotfiles (rho, momx, momy, x, y).

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

static_assert(AMREX_SPACEDIM == 2, "StarDisk is a 2D-only test");

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

struct DiskParams {
    Real a         = 1.0;   // star radius (inner disk boundary r = a)
    Real r_out     = 4.0;   // outer disk radius
    Real stretch   = 1.0;   // radial grading exponent (>1 clusters cells inward)
    Real gm        = 1.0;   // point mass of the star
    Real cs        = 0.1;   // isothermal sound speed
    Real rho0      = 1.0;   // density normalization at r = a
    Real q         = 1.0;   // density power-law index, rho ~ r^-q
    Real rho_floor = 1.e-8; // density floor
    Real soft      = 0.0;   // gravitational softening length
    Real cfl       = 0.5;
    int  n_phi     = 64;    // tangential cells per block
    int  n_r       = 64;    // radial cells per block
};

// Physical position of vertex (i, j) for a block with orientation theta_b.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuArray<Real, 2> vertex_pos (int i, int j, Real theta_b, DiskParams const& p) {
    const Real dphi = (0.5_rt * M_PI) / p.n_phi;
    const Real phi  = -0.25_rt * M_PI + i * dphi;
    const Real s    = Real(j) / Real(p.n_r);
    const Real rho  = p.a + (p.r_out - p.a) * std::pow(s, p.stretch);
    const Real ang  = theta_b + phi;
    return {rho * std::cos(ang), rho * std::sin(ang)};
}

// Rusanov (local Lax-Friedrichs) flux for isothermal Euler across a face with
// outward area vector (ax, ay). UL is the state inside, UR outside.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rusanov_flux (Real const* UL, Real const* UR, Real ax, Real ay, Real cs,
                   Real* F) {
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
                      - 0.5_rt * smax * (rhoR - rhoL));
    F[UMX]  = area * (0.5_rt * (UL[UMX] * unL + UR[UMX] * unR + (pL + pR) * nx)
                      - 0.5_rt * smax * (UR[UMX] - UL[UMX]));
    F[UMY]  = area * (0.5_rt * (UL[UMY] * unL + UR[UMY] * unR + (pL + pR) * ny)
                      - 0.5_rt * smax * (UR[UMY] - UL[UMY]));
}

class DiskBlock : public AmrCore {
  public:
    DiskBlock(Real theta_b_, Geometry const& level_0_geom, DiskParams const& p,
              AmrInfo const& amr_info = AmrInfo())
        : AmrCore(level_0_geom, amr_info), params{p}, theta_b{theta_b_} {
        AmrCore::InitFromScratch(0.0);
        InitData();
    }

    // Vertex coordinates, cell volumes, and the initial equilibrium disk.
    void InitData() {
        const DiskParams p = params;
        const Real tb = theta_b;
        for (MFIter mfi(xyv, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Array4<Real> v = xyv.array(mfi);
            ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::ignore_unused(k);
                const auto pos = vertex_pos(i, j, tb, p);
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
                // Shoelace. The (i=tangential-CCW, j=radial) parameterization
                // is left-handed, so the signed area comes out negative.
                vol_a(i, j, k) = 0.5_rt * std::abs((x00 * y10 - x10 * y00)
                                                 + (x10 * y11 - x11 * y10)
                                                 + (x11 * y01 - x01 * y11)
                                                 + (x01 * y00 - x00 * y01));
                const Real xc = 0.25_rt * (x00 + x10 + x11 + x01);
                const Real yc = 0.25_rt * (y00 + y10 + y11 + y01);
                const Real r  = std::sqrt(xc * xc + yc * yc);
                const Real rho = amrex::max(p.rho0 * std::pow(r / p.a, -p.q),
                                            p.rho_floor);
                // Radial equilibrium: v_phi^2 / r = GM/r^2 + (1/rho) dp/dr.
                const Real vphi = std::sqrt(amrex::max(p.gm / r - p.q * p.cs * p.cs,
                                                       0.0_rt));
                u(i, j, k, URHO) = rho;
                u(i, j, k, UMX)  = rho * vphi * (-yc / r);
                u(i, j, k, UMY)  = rho * vphi * ( xc / r);
            });
        }
    }

    // Physical boundary conditions on the radial (j) edges, applied to the
    // ghost cells of U: reflecting wall at the star (j-low), zero-gradient
    // outflow at the outer edge (j-high).
    void FillPhysicalBCs() {
        const Box& domain = Geom(0).Domain();
        for (MFIter mfi(U); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            Array4<Real> u = U.array(mfi);
            Array4<Real const> v = xyv.const_array(mfi);
            // Only ghost cells outside the physical domain are ours; ghost
            // cells between grids of this block belong to FillGhosts().
            const Box glo = amrex::adjCellLo(vbx, iy, nghost)
                & amrex::adjCellLo(domain, iy, nghost);
            const Box ghi = amrex::adjCellHi(vbx, iy, nghost)
                & amrex::adjCellHi(domain, iy, nghost);
            // j-low: wall. Ghost (i, j0-1) mirrors interior (i, j0) with the
            // normal momentum component flipped (free slip).
            if (glo.ok()) {
            ParallelFor(glo,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::ignore_unused(k);
                // Outward area vector of the j-low face of cell (i, j+1),
                // left-handed grid: A = (-ey, ex) for edge (i,j+1)->(i+1,j+1).
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
            // j-high: zero-gradient outflow.
            if (ghi.ok()) {
            ParallelFor(ghi,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::ignore_unused(k);
                for (int n = 0; n < ncomp; ++n) {
                    u(i, j, k, n) = u(i, j - 1, k, n);
                }
            });
            }
        }
    }

    // One explicit Euler step: Rusanov fluxes through the four curvilinear
    // faces plus the point-mass gravity source.
    void Advance(Real dt) {
        const DiskParams p = params;
        for (MFIter mfi(U, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Array4<Real const> u = U.const_array(mfi);
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
                Real F[ncomp];
                Real fsum[ncomp] = {0.0_rt, 0.0_rt, 0.0_rt};
                for (int n = 0; n < ncomp; ++n) { UL[n] = u(i, j, k, n); }

                for (int n = 0; n < ncomp; ++n) { UR[n] = u(i - 1, j, k, n); }
                rusanov_flux(UL, UR, Ailo[0], Ailo[1], p.cs, F);
                for (int n = 0; n < ncomp; ++n) { fsum[n] += F[n]; }

                for (int n = 0; n < ncomp; ++n) { UR[n] = u(i + 1, j, k, n); }
                rusanov_flux(UL, UR, Aihi[0], Aihi[1], p.cs, F);
                for (int n = 0; n < ncomp; ++n) { fsum[n] += F[n]; }

                for (int n = 0; n < ncomp; ++n) { UR[n] = u(i, j - 1, k, n); }
                rusanov_flux(UL, UR, Ajlo[0], Ajlo[1], p.cs, F);
                for (int n = 0; n < ncomp; ++n) { fsum[n] += F[n]; }

                for (int n = 0; n < ncomp; ++n) { UR[n] = u(i, j + 1, k, n); }
                rusanov_flux(UL, UR, Ajhi[0], Ajhi[1], p.cs, F);
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
                // Point-mass gravity at the cell center.
                const Real xc = 0.25_rt * (x00 + x10 + x11 + x01);
                const Real yc = 0.25_rt * (y00 + y10 + y11 + y01);
                const Real r2 = xc * xc + yc * yc + p.soft * p.soft;
                const Real g_over_r = p.gm / (r2 * std::sqrt(r2));
                mx_new -= dt * rho_new * g_over_r * xc;
                my_new -= dt * rho_new * g_over_r * yc;

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

    // CFL-limited timestep for this block: V / sum_faces (|u.n| + cs) |A|.
    Real ComputeDt() const {
        const DiskParams p = params;
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

    MultiFab U{};      // (rho, mom_x, mom_y), cell-centered, one ghost layer
    MultiFab Unew{};
    MultiFab xyv{};    // physical vertex coordinates, nodal, 2 components
    MultiFab vol{};    // cell volumes
    DiskParams params{};
    Real theta_b{};    // block orientation angle

  private:
    void ErrorEst(int, ::amrex::TagBoxArray&, Real, int) override {
        throw std::runtime_error("single-level only");
    }
    void MakeNewLevelFromScratch(int level, Real, const ::amrex::BoxArray& ba,
                                 const ::amrex::DistributionMapping& dm) override {
        if (level > 0) throw std::runtime_error("single-level only");
        U.define(ba, dm, ncomp, nghost);
        Unew.define(ba, dm, ncomp, nghost);
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
        xyv.clear();
        vol.clear();
    }
};

using namespace NonLocalBC;

// Aligned (no rotation, no flip, no component swap) one-sided fill of U from
// src -> dest across a tangential seam.
struct AlignedBoundaryFn {
    DiskBlock* dest;
    const DiskBlock* src;
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
    std::vector<AlignedBoundaryFn> boundaries;
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

void WritePlotfile(const DiskBlock& block, const char* dir, Real t, int step) {
    static const Vector<std::string> varnames{"rho", "momx", "momy", "x", "y"};
    char buf[256];
    snprintf(buf, sizeof(buf), "StarDisk/%s/plt%04d", dir, step);
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

// Tangential seam ghost slabs: the i-low (i = -1) and i-high (i = n_phi)
// ghost columns of the block's domain.
Box seam_ghost_lo(const Box& domain) {
    return Box{IntVect{AMREX_D_DECL(domain.smallEnd(ix) - 1, domain.smallEnd(iy), 0)},
               IntVect{AMREX_D_DECL(domain.smallEnd(ix) - 1, domain.bigEnd(iy), 0)}};
}
Box seam_ghost_hi(const Box& domain) {
    return Box{IntVect{AMREX_D_DECL(domain.bigEnd(ix) + 1, domain.smallEnd(iy), 0)},
               IntVect{AMREX_D_DECL(domain.bigEnd(ix) + 1, domain.bigEnd(iy), 0)}};
}

} // namespace

void MyMain() {
    DiskParams p;
    Real stop_time = 2.0_rt * M_PI;  // one orbit at the inner edge
    int max_steps = 1000000;
    int plot_int = 200;
    int print_int = 10;
    {
        ParmParse pp("disk");
        pp.query("a", p.a);
        pp.query("r_out", p.r_out);
        pp.query("stretch", p.stretch);
        pp.query("gm", p.gm);
        pp.query("cs", p.cs);
        pp.query("rho0", p.rho0);
        pp.query("q", p.q);
        pp.query("rho_floor", p.rho_floor);
        pp.query("soft", p.soft);
        pp.query("cfl", p.cfl);
        pp.query("n_phi", p.n_phi);
        pp.query("n_r", p.n_r);
        pp.query("stop_time", stop_time);
        pp.query("max_steps", max_steps);
        pp.query("plot_int", plot_int);
        pp.query("print_int", print_int);
    }

    // All four blocks share the same logical domain box; the physical
    // RealBox is the curved block's bounding box (used for plot axes only).
    Box domain(IntVect{}, IntVect{AMREX_D_DECL(p.n_phi - 1, p.n_r - 1, 0)});
    const int N = p.n_phi;
    Array<int, AMREX_SPACEDIM> per_none{AMREX_D_DECL(0, 0, 0)};

    auto make_geom = [&](Real theta_b) {
        std::array<Real, AMREX_SPACEDIM> lo{AMREX_D_DECL(std::numeric_limits<Real>::max(),
                                                         std::numeric_limits<Real>::max(),
                                                         std::numeric_limits<Real>::max())};
        std::array<Real, AMREX_SPACEDIM> hi{AMREX_D_DECL(std::numeric_limits<Real>::lowest(),
                                                         std::numeric_limits<Real>::lowest(),
                                                         std::numeric_limits<Real>::lowest())};
        for (int j = 0; j <= p.n_r; ++j) {
            for (int i = 0; i <= p.n_phi; ++i) {
                const auto pos = vertex_pos(i, j, theta_b, p);
                lo[0] = amrex::min(lo[0], pos[0]);
                lo[1] = amrex::min(lo[1], pos[1]);
                hi[0] = amrex::max(hi[0], pos[0]);
                hi[1] = amrex::max(hi[1], pos[1]);
            }
        }
        return Geometry{domain, RealBox{lo, hi}, CoordSys::cartesian, per_none};
    };

    AmrInfo amr_info{};

    DiskBlock be(0.0_rt,           make_geom(0.0_rt),           p, amr_info);
    DiskBlock bn(0.5_rt * M_PI,    make_geom(0.5_rt * M_PI),    p, amr_info);
    DiskBlock bw(1.0_rt * M_PI,    make_geom(1.0_rt * M_PI),    p, amr_info);
    DiskBlock bs(1.5_rt * M_PI,    make_geom(1.5_rt * M_PI),    p, amr_info);

    // Aligned-block index mapping: src_idx[d] = dest_idx[d] - offset[d].
    auto make_dtos = [](IntVect off) {
        MultiBlockIndexMapping d{};
        d.permutation = IntVect{AMREX_D_DECL(0, 1, 2)};
        d.sign        = IntVect{AMREX_D_DECL(1, 1, 1)};
        d.offset      = off;
        return d;
    };

    // Tangential seams close the annulus: E -> N -> W -> S -> E. Ghost
    // column i = N of a block maps to interior column i = 0 of the next
    // (offset = (N, 0)); ghost column i = -1 maps to interior i = N-1 of the
    // previous (offset = (-N, 0)). The radial index is unchanged.
    std::vector<AlignedBoundaryFn> bnd;
    bnd.push_back({&be, &bn, make_dtos(N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                   seam_ghost_hi(domain)});   // E.right <- N.left
    bnd.push_back({&bn, &be, make_dtos(-N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                   seam_ghost_lo(domain)});   // N.left  <- E.right
    bnd.push_back({&bn, &bw, make_dtos(N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                   seam_ghost_hi(domain)});   // N.right <- W.left
    bnd.push_back({&bw, &bn, make_dtos(-N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                   seam_ghost_lo(domain)});   // W.left  <- N.right
    bnd.push_back({&bw, &bs, make_dtos(N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                   seam_ghost_hi(domain)});   // W.right <- S.left
    bnd.push_back({&bs, &bw, make_dtos(-N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                   seam_ghost_lo(domain)});   // S.left  <- W.right
    bnd.push_back({&bs, &be, make_dtos(N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                   seam_ghost_hi(domain)});   // S.right <- E.left
    bnd.push_back({&be, &bs, make_dtos(-N * IntVect{AMREX_D_DECL(1, 0, 0)}),
                   seam_ghost_lo(domain)});   // E.left  <- S.right

    FillBoundaryFn FillSeams{std::move(bnd)};

    DiskBlock* blocks[4] = {&be, &bn, &bw, &bs};
    const char* names[4] = {"e", "n", "w", "s"};

    int step = 0;
    Real t = 0.0_rt;
    const Real mass0 = [&] {
        Real m = 0.0_rt;
        for (auto* b : blocks) m += b->TotalMass();
        return m;
    }();
    amrex::Print().SetPrecision(10)
        << "Initial total mass: " << mass0 << '\n';

    for (int b = 0; b < 4; ++b) WritePlotfile(*blocks[b], names[b], t, step);

    while (t < stop_time && step < max_steps) {
        Real dt = blocks[0]->ComputeDt();
        for (int b = 1; b < 4; ++b) dt = amrex::min(dt, blocks[b]->ComputeDt());
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
            Real mass = 0.0_rt;
            for (auto* b : blocks) mass += b->TotalMass();
            amrex::Print().SetPrecision(10)
                << "Step #" << step << ", t = " << t << ", dt = " << dt
                << ", mass drift = " << (mass - mass0) / mass0;
            for (int b = 0; b < 4; ++b) {
                const auto [rlo, rhi] = blocks[b]->RhoMinMax();
                amrex::Print() << ", " << names[b] << "_rho in [" << rlo
                               << ", " << rhi << "]";
            }
            amrex::Print() << '\n';
            for (int b = 0; b < 4; ++b) {
                amrex::Print().SetPrecision(10)
                    << "    " << names[b] << ": mass = " << blocks[b]->TotalMass()
                    << ", NaN cells = " << blocks[b]->NaNCount() << '\n';
            }
        }
        if (step % plot_int == 0 || t >= stop_time) {
            for (int b = 0; b < 4; ++b) WritePlotfile(*blocks[b], names[b], t, step);
        }
    }

    // A test should fail loudly: NaN anywhere means a broken boundary fill.
    Long nan_total = 0;
    for (auto* b : blocks) nan_total += b->NaNCount();
    if (nan_total > 0) {
        amrex::Abort("StarDisk: NaN detected in final state");
    }
}
