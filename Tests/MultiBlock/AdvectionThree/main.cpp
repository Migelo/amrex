// 2D three-block advection example using AMReX NonLocalBC.
//
// Layout (each block is a unit-square 64x64 in its own physical frame, drawn
// here just by adjacency, not by physical x-coordinate):
//
//          +-------+
//          |  B3   |     B1.top    <-> B3.bottom
//          +-------+     B1.bottom <-> B3.top
//          +---+---+---+
//          | B1| B2| B3|  B1.right <-> B2.left
//          +---+---+---+  B2.right <-> B3.left
//          +-------+     B2 is y-self-periodic
//          |  B3   |     B1.left, B3.right are open (no fill -> zero ghosts)
//          +-------+
//
// Each block carries one passive scalar `Mass`. Block 1 advects with
// velocity (+1, 0), block 2 with (+1, 0), block 3 with (0, +1).
// First-order upwind + Godunov dimensional splitting.
//
// This is the same machinery as Tests/MultiBlock/Advection (two blocks),
// extended to three and to 2D, exercising aligned (non-rotated) block
// stitching plus open and self-periodic boundaries.

#include "AMReX_NonLocalBC.H"

#include "AMReX.H"
#include "AMReX_AmrCore.H"
#include "AMReX_MultiFab.H"
#include "AMReX_PlotFileUtil.H"
#include "AMReX_REAL.H"

using namespace amrex;

void MyMain();

int main(int argc, char** argv) {
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#else
    amrex::ignore_unused(argc, argv);
#endif
    amrex::Initialize(MPI_COMM_WORLD, std::cout, std::cerr,
                      [](const char* msg) { throw std::runtime_error(msg); });
    MyMain();
    amrex::Finalize();
#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
}

enum idirs { ix, iy };
enum num_components { one_component = 1 };

static constexpr IntVect e_x = IntVect::TheDimensionVector(ix);
static constexpr IntVect e_y = IntVect::TheDimensionVector(iy);

class AdvectionAmrCore : public AmrCore {
  public:
    AdvectionAmrCore(Direction vel, Geometry const& level_0_geom,
                     AmrInfo const& amr_info = AmrInfo())
        : AmrCore(level_0_geom, amr_info), velocity{vel} {
        AmrCore::InitFromScratch(0.0);
        InitData();
    }

    // Initial blob in the lower-left interior of each block.
    void InitData() {
        const auto problo = Geom(0).ProbLoArray();
        const auto dx = Geom(0).CellSizeArray();
        const Real x0 = problo[0] + 0.2_rt;
        const Real y0 = problo[1] + 0.2_rt;
        for (MFIter mfi(mass, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Array4<Real> m = mass.array(mfi);
            ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Real x = problo[0] + (0.5_rt + i) * dx[0];
                Real y = problo[1] + (0.5_rt + j) * dx[1];
                Real r2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
                m(i, j, k) = r2 < 0.01_rt ? 1_rt : 0_rt;
            });
        }
    }

    void AdvanceInTime(double dt) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            DoOperatorSplitStep(dt, static_cast<Direction>(d));
        }
    }

    void DoOperatorSplitStep(double dt, Direction dir) {
        const double dx = Geom(0).CellSize(0);
        const double a_dt_over_dx = dt / dx * (velocity == dir);
        if (dir == Direction::x) {
            for (MFIter mfi(mass); mfi.isValid(); ++mfi) {
                Array4<Real> m = mass.array(mfi);
                Array4<Real> next = mass_next.array(mfi);
                ParallelFor(mfi.growntilebox(e_y), int(one_component),
                            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                                next(i, j, k, n) = m(i, j, k, n)
                                    - a_dt_over_dx * (m(i, j, k, n) - m(i - 1, j, k, n));
                            });
            }
            std::swap(mass, mass_next);
        } else if (dir == Direction::y) {
            for (MFIter mfi(mass); mfi.isValid(); ++mfi) {
                Array4<Real> m = mass.array(mfi);
                Array4<Real> next = mass_next.array(mfi);
                ParallelFor(mfi.growntilebox(e_x), int(one_component),
                            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) {
                                next(i, j, k, n) = m(i, j, k, n)
                                    - a_dt_over_dx * (m(i, j, k, n) - m(i, j - 1, k, n));
                            });
            }
            std::swap(mass, mass_next);
        }
    }

    MultiFab mass{};
    MultiFab mass_next{};
    Direction velocity{};

  private:
    void ErrorEst(int, ::amrex::TagBoxArray&, Real, int) override {
        throw std::runtime_error("single-level only");
    }
    void MakeNewLevelFromScratch(int level, Real, const ::amrex::BoxArray& ba,
                                 const ::amrex::DistributionMapping& dm) override {
        if (level > 0) throw std::runtime_error("single-level only");
        const IntVect ngrow{AMREX_D_DECL(1, 1, 0)};
        mass.define(ba, dm, one_component, ngrow);
        mass_next.define(ba, dm, one_component, ngrow);
        // Zero out ghosts so "open" boundaries (never filled) stay well-defined.
        mass.setVal(0.0);
        mass_next.setVal(0.0);
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
        mass.clear();
    }
};

using namespace NonLocalBC;

// Aligned (no rotation, no flip, no component swap) one-sided fill from src->dest.
struct AlignedBoundaryFn {
    AdvectionAmrCore* dest;
    const AdvectionAmrCore* src;
    MultiBlockIndexMapping dtos;
    Box boundary_to_fill;
    std::unique_ptr<MultiBlockCommMetaData> cmd{};
    FabArrayBase::BDKey cached_dest_bd_key{};
    FabArrayBase::BDKey cached_src_bd_key{};
    ApplyDtosAndProjectionOnReciever<MultiBlockIndexMapping, MapComponents<Identity>>
        packing{PackComponents{0, 0, one_component}, dtos};

    AMREX_NODISCARD CommHandler FillBoundary_nowait() {
        if (!cmd || cached_dest_bd_key != dest->mass.getBDKey()
            || cached_src_bd_key != src->mass.getBDKey()) {
            cmd = std::make_unique<MultiBlockCommMetaData>(
                dest->mass, boundary_to_fill, src->mass, dest->mass.nGrowVect(), dtos);
            cached_dest_bd_key = dest->mass.getBDKey();
            cached_src_bd_key = src->mass.getBDKey();
        }
        return ParallelCopy_nowait(no_local_copy, dest->mass, src->mass, *cmd, packing);
    }
    void FillBoundary_do_local_copy() const {
        if (cmd->m_LocTags && !cmd->m_LocTags->empty()) {
            LocalCopy(packing, dest->mass, src->mass, *cmd->m_LocTags);
        }
    }
    void FillBoundary_finish(CommHandler handler) const {
        ParallelCopy_finish(dest->mass, std::move(handler), *cmd, packing); // NOLINT
    }
};

struct FillBoundaryFn {
    std::vector<AlignedBoundaryFn> boundaries;
    void operator()(std::initializer_list<AdvectionAmrCore*> cores) {
        // Geometry-level fills first (handles intra-block periodicity, e.g. B2 in y).
        for (auto* c : cores) {
            c->mass.FillBoundary(c->Geom(0).periodicity());
        }
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

void WritePlotfile(const AdvectionAmrCore& core, const char* dir, Real t, int step) {
    static const Vector<std::string> varnames{"Mass"};
    char buf[256];
    snprintf(buf, sizeof(buf), "MultiBlock3/%s/plt%04d", dir, step);
    Vector<const MultiFab*> mf{&core.mass};
    Vector<Geometry> geoms{core.Geom(0)};
    Vector<int> level_steps{step};
    Vector<IntVect> ref_ratio{};
    WriteMultiLevelPlotfile(buf, 1, mf, varnames, geoms, t, level_steps, ref_ratio);
}

// Helpers to construct the boundary patch box (ghost slab in the dest block).
static Box face_ghost_x_high(const Box& domain) { // i = N+1, j = -1..N+1
    return grow(shift(Box{domain.bigEnd(ix) * e_x, domain.bigEnd()}, e_x), e_y);
}
static Box face_ghost_x_low(const Box& domain) { // i = -1
    return grow(shift(Box{domain.smallEnd(),
                          domain.bigEnd() - domain.bigEnd(ix) * e_x}, -e_x), e_y);
}
static Box face_ghost_y_high(const Box& domain) { // j = N+1
    return grow(shift(Box{domain.bigEnd(iy) * e_y, domain.bigEnd()}, e_y), e_x);
}
static Box face_ghost_y_low(const Box& domain) { // j = -1
    return grow(shift(Box{domain.smallEnd(),
                          domain.bigEnd() - domain.bigEnd(iy) * e_y}, -e_y), e_x);
}

void MyMain() {
    Box domain(IntVect{}, IntVect{AMREX_D_DECL(63, 63, 0)});
    const int N = domain.bigEnd(ix); // 63

    // Three blocks share the same logical domain box; physical RealBoxes are
    // chosen for nicer plotfile geometry (no other meaning).
    RealBox rb1{{AMREX_D_DECL(0.0_rt, 0.0_rt, 0_rt)},
                {AMREX_D_DECL(1.0_rt, 1.0_rt, 1_rt)}};
    RealBox rb2{{AMREX_D_DECL(1.2_rt, 0.0_rt, 0_rt)},
                {AMREX_D_DECL(2.2_rt, 1.0_rt, 1_rt)}};
    RealBox rb3{{AMREX_D_DECL(2.4_rt, 0.0_rt, 0_rt)},
                {AMREX_D_DECL(3.4_rt, 1.0_rt, 1_rt)}};

    Array<int, AMREX_SPACEDIM> per_none{AMREX_D_DECL(0, 0, 0)};
    Array<int, AMREX_SPACEDIM> per_y{AMREX_D_DECL(0, 1, 0)}; // B2: y-self-periodic
    Geometry g1{domain, rb1, CoordSys::cartesian, per_none};
    Geometry g2{domain, rb2, CoordSys::cartesian, per_y};
    Geometry g3{domain, rb3, CoordSys::cartesian, per_none};

    AmrInfo amr_info{};
#if AMREX_SPACEDIM > 2
    amr_info.blocking_factor[0][2] = 1;
#endif

    AdvectionAmrCore b1(Direction::x, g1, amr_info);
    AdvectionAmrCore b2(Direction::x, g2, amr_info);
    AdvectionAmrCore b3(Direction::y, g3, amr_info);

    // Aligned-block index mapping: src_idx[d] = dest_idx[d] - offset[d].
    auto make_dtos = [](IntVect off) {
        MultiBlockIndexMapping d{};
        d.permutation = IntVect{AMREX_D_DECL(0, 1, 2)};
        d.sign        = IntVect{AMREX_D_DECL(1, 1, 1)};
        d.offset      = off;
        return d;
    };

    std::vector<AlignedBoundaryFn> bnd;

    // ---- Horizontal chain: B1 <-> B2 <-> B3 ----
    // B1.right ghost (i = N+1) <- B2 interior at i = 0  =>  offset = (N+1, 0)
    bnd.push_back({&b1, &b2, make_dtos((N + 1) * e_x), face_ghost_x_high(domain)});
    // B2.left ghost (i = -1) <- B1 interior at i = N    =>  offset = (-(N+1), 0)
    bnd.push_back({&b2, &b1, make_dtos(-(N + 1) * e_x), face_ghost_x_low(domain)});
    // B2.right ghost (i = N+1) <- B3 interior at i = 0  =>  offset = (N+1, 0)
    bnd.push_back({&b2, &b3, make_dtos((N + 1) * e_x), face_ghost_x_high(domain)});
    // B3.left ghost (i = -1) <- B2 interior at i = N    =>  offset = (-(N+1), 0)
    bnd.push_back({&b3, &b2, make_dtos(-(N + 1) * e_x), face_ghost_x_low(domain)});

    // ---- Vertical loop: B1 <-> B3 (top/bottom both ways) ----
    // B1.top    (j = N+1) <- B3 j = 0   =>  offset = (0, N+1)
    bnd.push_back({&b1, &b3, make_dtos((N + 1) * e_y), face_ghost_y_high(domain)});
    // B3.bottom (j = -1)  <- B1 j = N   =>  offset = (0, -(N+1))
    bnd.push_back({&b3, &b1, make_dtos(-(N + 1) * e_y), face_ghost_y_low(domain)});
    // B1.bottom (j = -1)  <- B3 j = N   =>  offset = (0, -(N+1))
    bnd.push_back({&b1, &b3, make_dtos(-(N + 1) * e_y), face_ghost_y_low(domain)});
    // B3.top    (j = N+1) <- B1 j = 0   =>  offset = (0, N+1)
    bnd.push_back({&b3, &b1, make_dtos((N + 1) * e_y), face_ghost_y_high(domain)});

    // B1.left and B3.right are "open" -- no entry. Ghosts stay zero.

    FillBoundaryFn FillBoundary{std::move(bnd)};

    int step = 0;
    const Real dx = std::min({g1.CellSize(0), g2.CellSize(0), g3.CellSize(0),
                              g1.CellSize(1), g2.CellSize(1), g3.CellSize(1)});
    const Real cfl = 1.0_rt;
    const Real dt = cfl * dx;
    const Real final_time = 4.0_rt;
    Real t = 0.0_rt;

    WritePlotfile(b1, "b1", t, step);
    WritePlotfile(b2, "b2", t, step);
    WritePlotfile(b3, "b3", t, step);

    while (t < final_time) {
        FillBoundary({&b1, &b2, &b3});
        b1.AdvanceInTime(dt);
        b2.AdvanceInTime(dt);
        b3.AdvanceInTime(dt);
        t += dt;
        ++step;
        amrex::Print() << "Step #" << step << ", t = " << t << '\n';
        WritePlotfile(b1, "b1", t, step);
        WritePlotfile(b2, "b2", t, step);
        WritePlotfile(b3, "b3", t, step);
    }
}
