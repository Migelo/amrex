# AMReX Agents Guide

Use this guide whenever you orchestrate explorers/workers inside the AMReX repository. It keeps PR reviews, bug hunts, new features, docs, and user support consistent for contributors and AI copilots alike.

## Purpose & Personas

- **AMReX developers** – follow a predictable workflow so fixes/features land quickly and safely.
- **AMReX users** – receive accurate build/test/docs guidance sourced from the repo.

## Repository Map

- `Src/` – Core C++/Fortran implementation. `Amr*` manage hierarchy/regridding, `Base` hosts runtime utilities, and folders such as `Boundary`, `EB`, `LinearSolvers`, `Particle`, `FFT` contain subsystem code (check local `CMakeLists.txt` for toggles).
- `Tests/` – Organized by topic; when `AMReX_ENABLE_TESTS=ON`, they expose `ctest` targets. Many tutorials/tests also ship a `GNUmakefile` for standalone runs.
- `Docs/` – Authoritative Sphinx sources under `Docs/sphinx_documentation/`; Doxygen lives under `Docs/Doxygen/`.
- `Tools/` – Shared CMake modules, GNU make helpers, scripts.

## Build & Test Reference

- **Default CMake flow**:
  ```bash
  cmake -S . -B build \
    -DAMReX_ENABLE_TESTS=ON \
    -DAMReX_TEST_TYPE=Small
  cmake --build build -j
  ctest --test-dir build --output-on-failure
  ```
  Toggle options from `Docs/sphinx_documentation/source/BuildingAMReX.rst` (“Customization options”) and `Tools/CMake/AMReXOptions.cmake` for GPUs, dimensions, debug flags, etc.
- **Targeted builds/tests**: `cmake --build build -j --target <name>` for individual executables; `ctest --test-dir build -R <regex>` (or `ctest -R <regex>`) to rerun impacted cases only.
- **GNUmakefile workflows**: When a directory ships a `GNUmakefile`, `cd` there and run `make -j` with required variables (e.g., `DIM`, `USE_MPI`, `USE_CUDA`, `COMP`) as documented in `Docs/sphinx_documentation/source/BuildingAMReX.rst` and `Tools/GNUMake/README.md`.
- **Log everything**: Capture exact commands plus pass/fail output in PR descriptions or linked issues so reviewers can reproduce without guessing.

## Hard Rules & Defaults

- Work on short-lived topic branches; never commit directly to `development`. Rebase before opening or updating a PR.
- Plan before you edit. Define scope, owners, and validation steps, and use `rg` (or your editor) for fast searches to avoid wandering.
- Choose a build workflow from the “Build & Test Reference” and stick to it for the task; don’t restate commands elsewhere—link back and record the actual invocations you ran.
- Update user-facing docs under `Docs/sphinx_documentation/` whenever behavior changes, and mention the doc edits in the same PR.
- When adding Doxygen comments or Sphinx docs, write for AMReX users: explain behavior and usage, omit unnecessary implementation details, and keep the tone instructional instead of developer-oriented.
- Never log work inside `CHANGES.md`; release notes are curated separately, so leave that file untouched unless maintainers ask for a release update.
- Hand-off context inside PRs/issues (commands run, results, TODOs). Personal scratchpads are optional, local, and should be pruned frequently.
- AI delegation: assign disjoint write scopes, keep subtasks bounded, and serialize conflicting edits when necessary. Review every AI-generated diff like a teammate’s patch.
- Safety: never paste secrets, private datasets, or PII into prompts; scrub logs before sharing; cite sources and licensing when in doubt by escalating to maintainers.
- Testing accountability: every substantive change must state which tests ran (and command lines) in the PR/issue.
- Historical hotspots: mirrored kernels and dimension-specific paths often drift—compare siblings whenever you touch them.

## GPU Lambda Safety

- `AMREX_GPU_DEVICE`/`ParallelFor` lambdas run on the GPU, so never capture host-only pointers (e.g., `Geometry::CellSize()`, `Geometry::ProbLo()`). Take the device-safe views first—e.g., `auto const dx = geom.CellSizeArray();`, `auto const problo = geom.ProbLoArray();`—and pass those by value into the lambda.

## Task Playbooks

### PR & Bug Review
1. **Sync & scope** – Update the branch, read the PR/issue description, and note ownership expectations.
2. **Reproduce** – Follow the reporter’s steps or minimal `ctest`/`make` invocations from the Build & Test Reference.
3. **Inspect diffs** – Check style, compare mirrored kernels/dimensional variants, and look for missing constant updates.
4. **Verify** – Run the focused tests (GPU flags, `-DAMReX_TEST_TYPE=Small`, GNU make targets, etc.) and capture output.
5. **Report** – List blockers first, cite files/lines, and record the commands/tests in the PR.
6. **Follow-ups** – Document remaining work in the PR/issue; use a scratchpad only for personal reminders.

### Feature or Fix Implementation
1. **Confirm requirements** – Capture physics context, success metrics, and acceptance tests from the originating ticket.
2. **Configure builds** – Apply the Build & Test Reference workflow that matches the directory; use targeted targets/tests for fast iteration.
3. **Implement carefully** – Touch only owned files, note tricky code with short comments, and maintain coding-style guidelines.
4. **Update docs** – When behavior or UX changes, add/adjust content in the relevant `Docs/sphinx_documentation/` section (User’s Guide, Technical Reference, tutorials, etc.).
5. **Validate** – Run and record the necessary `ctest`/`make` commands plus any benchmarking or profiling relevant to the change.
6. **Hand off** – Summarize status, remaining risks, and next steps in the PR description; attach log excerpts only if scrubbed.

### Documentation & Tutorial Updates
- Mirror the hierarchy noted in `README.md` (User’s Guide, Example Codes, Guided Tutorials, Technical Reference) so published docs stay synchronized.
- Surface any new build knobs/workflows in `Docs/sphinx_documentation/source/BuildingAMReX.rst`.
- Reference runnable examples in `Tutorials/README.md` or `https://github.com/AMReX-Codes/amrex-tutorials` when guiding users.

## Guidance for AMReX Users Working with Agents

- **Orientation** – Summarize capabilities using the “Overview,” “Features,” and “Documentation” sections of `README.md`, then link users to the best resource (User’s Guide, Example Codes, Guided Tutorials, Technical Reference).
- **Build help** – Walk through the Build & Test Reference commands, clarifying how to toggle features via `-DVAR=value` or `make` variables the tutorial expects.
- **Learning resources** – Point to the tutorials repo plus any slides/videos referenced near the Documentation section of `README.md`.
- **Authoritative sources** – Read directly from `Docs/sphinx_documentation/` when answering doc questions so citations match the published site.
- **Support channels** – Encourage GitHub Discussions/issues for unresolved questions and remind users that contributions follow `CONTRIBUTING.md`.

## Optional Personal Notes

A local scratchpad can help you capture quick reminders, but keep it untracked, short-lived, and private. Anything others need to know belongs in PR descriptions, issues, or review comments.

## Quick Checklist

1. On a topic branch rebased on `development`? If not, fix it.
2. Do you have a written plan that names owners, validation, and Build & Test commands? If not, write one before editing.
3. Are tests/docs updated and the exact commands/results logged in the PR/issue? If not, add them.
4. Are delegation, safety, and hand-off notes captured in canonical threads (not just scratchpads)? If not, update them now.

## Current Situation (2026-07-19, branch `multiblock-three-block-advection`)

Working state for the StarDisk multi-block test. This section is a local scratchpad; prune before merging upstream.

### Goal

Adapt the `Tests/MultiBlock/AdvectionThree` setup (commit `31602db9c`) into a 2D star + accretion disk test: a fixed circular star `r < a` (masked, not evolved; round since 2026-07-19, previously the square `[-a,a]^2` via `r_in(phi)=a/cos(phi)`) surrounded by 4 curvilinear annulus blocks (E, N, W, S). Design decisions (user): **isothermal Euler hydrodynamics**, **fixed star region**.

### Environment

- No system compiler on this NixOS workstation. `shell.nix` + `.envrc` at repo root (tracked since `5ca8de56`; only `.direnv/` is untracked) provide gcc 15.2, gfortran, openmpi 5.0.10, cmake, python3 via direnv.
- **Weak VPS: build with max 2 cores (`make -j2`), run with max 2 MPI ranks** (user directive 2026-07-19).
- Build: `cd Tests/MultiBlock/StarDisk && direnv exec . make -j2` (~2 min from clean; ~4.5 min for `DEBUG=TRUE`).
- Executables: `main2d.gnu.MPI.ex` (opt) and `main2d.gnu.DEBUG.MPI.ex` (built with `DEBUG=TRUE`).
- Run: `direnv exec . mpirun -n 2 ./main2d.gnu.MPI.ex disk.max_steps=1 disk.print_int=1 disk.plot_int=1000000` (ParmParse overrides work; see fixed bug #4).
- ctest: one-time configure `direnv exec . cmake -S . -B build -DAMReX_SPACEDIM=2 -DAMReX_ENABLE_TESTS=ON -DAMReX_TEST_TYPE=All` (2D required for StarDisk/RocheBinary), then `direnv exec . ctest --test-dir build -R MultiBlock --output-on-failure`. New test dirs under Tests/MultiBlock are auto-registered via GLOB_RECURSE, but only after re-running cmake configure.
- Video pipeline: `cd Tests/MultiBlock/RocheBinary && ./run_roche_video.sh OUTNAME roche.key=value ...` runs the sim in an isolated `video_runs/OUTNAME/` workdir, renders one rho frame per plotfile step (render_frames.py, star masks auto-derived from m1/m2/sep/a), encodes h264 via system ffmpeg (env knobs FPS, CRF), and publishes to `~/stardisk-site/OUTNAME.mp4`. Verified end-to-end with q=0.5 (`roche.m2=0.5`): masks land on the asymmetric star positions (x=-1.333, +2.667).

### The test: `Tests/MultiBlock/StarDisk/` (committed)

`main.cpp`, `GNUmakefile` (DIM=2, MPI), `Make.package`, `CMakeLists.txt` (2D-only target; registered and passing in ctest).

- Blocks: 4 x `DiskBlock : AmrCore`, orientations theta_b = 0, pi/2, pi, 3pi/2. Logical (i=tangential CCW, j=radial). Mapping: `rho(j) = a + (R_out-a)*(j/n_r)^stretch`, `x = rho*(cos,sin)(theta_b+phi)`, i.e. each block is a polar sector of the circular annulus `a < r < R_out`. NOTE: (i=theta, j=r) ordering is **left-handed** -> shoelace volumes need `abs()`, outward face area vectors are `A=(-ey,ex)` (already fixed in code).
- Metrics from vertex coordinates (face area vectors, cell volumes); no analytic Jacobians. Isothermal Euler (rho, mom_x, mom_y), p = cs^2 rho, Rusanov flux per face, point-mass gravity, explicit Euler. IC: near-equilibrium Keplerian disk `rho ~ r^-q`, `v_phi = sqrt(GM/r - q*cs^2)`.
- BCs: inner j-edge reflecting wall (star surface), outer j-edge zero-gradient outflow, restricted to domain-edge ghost cells via `adjCellLo/Hi(vbx) & adjCellLo/Hi(domain)`.
- Seams: 4 tangential seams (E<->N<->W<->S<->E) = 8 one-sided `AlignedBoundaryFn` fills (offset-only `MultiBlockIndexMapping`, `PackComponents{0,0,3}`), closing the annulus with no periodicity anywhere.
- Kept diagnostics: per-print_int one-line summary (step, t, dt, mass drift, per-block rho ranges), per-block mass + NaN count, and a final `amrex::Abort` on any NaN (ctest pass/fail).

### Bugs fixed during bring-up

1. Left-handed parameterization (see above).
2. `seam_ghost_lo` built a 1-cell box instead of the full ghost column.
3. Missing intra-block `U.FillBoundary` for ghosts between grids of a block (`FillGhosts()`).
4. **ParmParse was dead**: `amrex::Initialize(MPI_COMM_WORLD, cout, cerr, handler)` overload never parses argv; switched to the argc/argv overload. (`AdvectionThree` has the same pattern; harmless there, fully hardcoded.)
5. `FillPhysicalBCs` applied wall/outflow to internal fab boundaries (fixed with domain intersection).
6. `RealBox` has no `RealVect` ctor -> `std::array`.

### Verified (1 MPI rank)

Runs to t=2pi (1116 steps). **Mass drift 2.5e-15 after step 1** -> scheme + seams exactly conservative. 4-fold block symmetry holds. Initial total mass 18.84908281 = 6*pi exactly (analytic for rho=r^-1 on the circular annulus 1<r<4). Slow +0.6% drift/orbit = physical flux through the open outer boundary during transient adjustment (not a conservation bug).

### RESOLVED (2026-07-18): seam fills broken with 2 MPI ranks

**Root cause (self-inflicted):** an instrumentation edit accidentally deleted the `FillGhosts()` call from the time loop. `AmrMesh::MakeBaseGrids` chops each 64x64 block into `NProcs` grids, so with 2 ranks every block has 2 fabs (split in j) and the inter-fab ghost ROWS (j=31/32, 128 cells/block) went unfilled; the flux kernel read stale garbage there. On 1 rank there is 1 fab/block (no inter-fab ghosts), so it looked perfect. The "128 seam cells" were the inter-fab rows, not seam columns; the NonLocalBC seam machinery was correct all along (probe-verified bitwise). Fix: call `FillGhosts()` (intra-block `U.FillBoundary`) **before** `FillSeams()` — the upstream MultiBlock idiom (`Tests/MultiBlock/Advection/main.cpp:187`).

**Slivers (secondary, benign):** `MultiBlockCommMetaData::define` grows both dst and src boxes by ngrow (Src/Base/AMReX_NonLocalBCImpl.H:362,369 via BoxArray::intersections, AMReX_BoxArray.cpp:1271), creating 1-cell cross-fab sliver overlaps that read the neighboring src fab's ghost ring and are applied (MPI finish) after the local main fill. Under the FillBoundary-first idiom those ghosts hold current valid data (bitwise identical to the main fill), so slivers are harmless. Not an upstream bug; a latent fragility worth knowing when writing new multi-block tests: **always FillBoundary intra-block ghosts before NonLocalBC seam fills**.

**AdvectionThree:** NOT broken on 2 ranks (earlier byte-wise plotfile diff was a fab-count artifact: 1 rank writes 1 fab/plotfile, 2 ranks write 2; cell-wise parse = bitwise identical).

**Verified:** 2-rank poison test 0 unfilled; 1- and 2-rank runs 0 NaN, bitwise-identical trajectories (mass drift ~2.5e-15 at step 1; +0.61% over one orbit = physical outflow); 4-fold symmetry; `ctest -R StarDisk` (cmake, AMReX_TEST_TYPE=All) passes on `mpiexec -n 2`. Re-verified 2026-07-19 for the circular star (1+2 ranks, full orbit, `stardisk_rho.png` shows a clean round hole, no seam artifacts).

### Status

Committed on `multiblock-three-block-advection` as `Tests/MultiBlock/StarDisk/` (main.cpp, GNUmakefile, Make.package, CMakeLists.txt). Debug instrumentation stripped; kept per-step mass/rho diagnostics + final NaN abort (ctest pass/fail). Untracked, uncommitted by design: `.direnv/` and this AGENTS.md section (`shell.nix` + `.envrc` are tracked since `5ca8de56`). Prune this section when stale.

### RocheBinary: two stars + bridge (IMPLEMENTED 2026-07-19)

**Handoff: `Tests/MultiBlock/RocheBinary/PLAN.md` is the full instructions doc for the next agent** (mission, hard rules, code map, do-not-break list, verification toolkit, publishing workflow, queued long-settling setups A–D). Read it before continuing this work.

Contact-binary extension of StarDisk: two fixed stars at (±2, 0), each wrapped in 4 polar-sector blocks (annulus 1<r<1.6), plus a 9th curvilinear **bridge block** (Coons patch) connecting the annuli through L1. Isothermal Euler in the corotating frame (Omega^2=G(m1+m2)/D^3): source = -rho grad(Phi_roche) - 2 Omega x (rho u). IC = exact hydrostatic isothermal atmosphere rho = rho_l1 exp(-(Phi-Phi_L1)/cs^2) at rest, L1 found by bisection on dPhi/dx (host). New test dir `Tests/MultiBlock/RocheBinary/` (main.cpp, GNUmakefile, Make.package, CMakeLists.txt, plot_roche.py, compare_ranks.py); ctest auto-registers via GLOB_RECURSE (test `MultiBlock_RocheBinary_2d`, needs AMReX_SPACEDIM=2 configure). StarDisk untouched. Build: `cd Tests/MultiBlock/RocheBinary && direnv exec . make -j2` (this VPS: max 2 cores, max 2 MPI ranks). Render: `python3 plot_roche.py 0 <step>` needs numpy+matplotlib — shell.nix python lacks them; use nix-shell with `(python3.withPackages (ps: [ps.numpy ps.matplotlib]))` over the same fetchTarball pin.

Key design points: 20 one-sided seam fills = 16 tangential (offset-only, as StarDisk) + bridge.jlo<->e1.jhi (offset-only) + **bridge.jhi<->w2.jhi with dtos sign=(-1,1)** (w2's i runs top->bottom, bridge's i bottom->top; offset=(n_phi-1, ±(n_bridge-n_r+1)); sign-flip precedent Tests/MultiBlock/Advection/main.cpp:256). Bridge mapping (i~+y, j~+x) is left-handed like the polar blocks, so the same vertex metrics code (abs shoelace, A=(-ey,ex)) serves all 9 blocks. Defaults: m1=m2=1, D=4, a=1, r_ann1=r_ann2=1.6 (Eggleton lobe 0.379*sep=1.52 contained; 2*r_ann<sep no overlap), cs=0.25, rho_l1=1e-2, floor 1e-8, stop_time=1.0 (~150 steps). Params under `roche.*`, incl. `perturb` (lobe-1 density boost drives L1 transfer), `poison_test` (NaN-poison ghosts + count unfilled face ghosts, must be 0).

**Per-star radii (q != 1)**: r_ann1/r_ann2 set each star's mesh radius independently; the primary's Roche lobe must fit (Eggleton: q=2 -> 0.44*D ~ 1.76+far side ~1.87, use r_ann1=1.85), and L1 must lie inside the bridge corridor (x1+r_ann1 < x_L1 < x2-r_ann2) — both asserted in derive().

**cs tuning rule**: rho_surf/rho_l1 = exp(Delta-Phi/cs^2) with Delta-Phi = Phi_L1 - Phi(star surface toward L1); keep it ~e^6 (~400). q=1/cs=0.25: Delta-Phi=0.365. q=2/cs=0.38: Delta-Phi1=0.856 -> contrast 376, Delta-Phi2=0.252 -> 5.7. Carrying cs=0.25 to q=2 gives contrast ~1e6 = broken run.

**Well-balanced flux (long-time survival)**: plain Rusanov dissipation erodes the steep hydrostatic profile on a diffusion time H^2/D ~ a few time units (v doubles every ~10 t, lobes evaporate by t~100; 20-orbit runs impossible). Fix: the Rusanov dissipative term acts on the deviation from the initial state, U-U0 (rusanov_flux takes U0L/U0R; centered pressure flux unchanged). In equilibrium the dissipation vanishes identically -> step-1 mass drift exactly 0, atmospheres survive 100+ orbits at truncation level. CRITICAL: U0's ghost cells must be filled with the SAME rules as U's ghosts (seam copies, wall mirrors, reservoir copies) — filling them analytically (the equilibrium profile continues inside the star hole) made the deviation nonzero at every physical BC face and turned the dissipation into a boundary mass pump (-82% in one step, v~2000). Solution: snapshot U0 interior+ghosts from U after the production fill pipeline ran once at t=0 (SnapshotInitGhosts, called in MyMain before the loop).

**Verified (1+2 ranks unless noted):** t=0 bboxes/volumes match analytic (polar 1.22510 vs pi(R^2-a^2)/4=1.22522 chord error; bridge 2.4701; lobe masses equal). Poison test 0 unfilled face ghosts (corners never filled, never read by 4-face stencil — excluded). **Step-1 mass drift -2.0e-16/+8.2e-16 (1/2 ranks)** -> all 20 seams + BCs exactly conservative. Static hold 100 steps: drift -6.0e-5, max|v| 0.024 stabilizing, lobe symmetry exact, L1 flux ~1e-19. Transfer run (perturb=0.05, t=4): drift -0.08%, max|v| 0.059 bounded, lobe1 -1.3e-3 net, bridge +4e-4, L1 flux ~+1e-4 toward star 2 (lobe1->bridge->lobe2 pathway, Coriolis-deflected circulation visible in drift panel). **1v2 ranks bitwise-identical state** (compare_ranks.py cell-wise plotfile diff, worst = 0; printed mass sums differ at 1e-16, reduction order). DEBUG build + 5-step run clean.

**Bugs fixed during bring-up (lessons for future multi-block tests):**
1. Nodal off-by-one in a diagnostics loop: xyv is nodal, mfi.validbox() IS the vertex box; adding +1 read garbage (trivial, caught by bbox report).
2. **FillPhysicalBCs overwrote seam-filled ghosts**: loop order FillGhosts -> FillSeams -> FillPhysicalBCs (StarDisk idiom) is only safe if BC edges and seam edges are disjoint. e1/w2 j-hi are seams (bridge), but the inherited polar j-hi zero-gradient BC overwrote the seam-filled ghosts -> one-sided flux across the seam -> +2.95e-6 mass/step created in the bridge. Fix: per-edge seam flag in BlockGeom, BCs skip seam edges. **Lesson: when adding seams to an existing block layout, audit every BC edge for overlap.**
3. **Zero-gradient outflow runaway in a tenuous hydrostatic atmosphere**: ghost copies of the interior let the discretely-imbalanced boundary shell collapse inward and the BC chases it: +25% mass and Mach>4 tenuous gas by t=4, dt collapsing. Fix: far-field reservoir BC (ghost = initial equilibrium state stored in Uinit) -> exchange driven only by genuine interior relaxation; drift -0.08%/t=4, bounded velocities. **Lesson: zero-gradient outflow is unsafe for steep stratified atmospheres; pin the boundary to equilibrium or use a reservoir.**

**Queued setup A (semi-detached) — DONE 2026-07-19:** new keys `roche.ic_mode` ("overcontact"|"semi_detached") and `roche.vmax` (velocity cap, 0=off). Semi-detached IC cuts `xc >= x_l1` to the floor (donor filled, accretor empty; L1-plane discontinuity). The donor-vs-vacuum jump drives an expansion-into-vacuum Riemann problem whose floor-cells get v=mom/rho ~ 5000, collapsing dt (~6e-7) and stalling; fix = velocity cap `roche.vmax=1.0` (4*cs) after the hydro update, preserving direction, transonic stream untouched. **vmax defaults 0 so overcontact is bitwise unchanged** (step-1 drift 0 even with vmax=1.0 — cap only touches momentum, never triggers in equilibrium). Verified q=1 cs=0.25: poison 0 (1+2 ranks), step-1 drift 0, 1v2 bitwise (worst=0), 4 MultiBlock ctests pass, DEBUG 5-step clean, 0 NaN to t=70 (2 orbits). lobe2 monotonic 4.9e-8 -> 7.6e-3 (0.7% of lobe1), L1 flux sustained ~5e-4, drift -2.2% at t=70 (open-edge outflow, ~1%/orbit). Run cut to 2 orbits per user ("until settled"); full 1-5% fill + L1-flux decay need ~5-10 more orbits (scheme stable). Published `~/stardisk-site/semidetached.mp4` + `semidetached_rho.png` + site section.

**Queued setup B (differential contact depth) — DONE 2026-07-19:** new keys `roche.rho_l1_1` / `roche.rho_l1_2` (sentinel <0 -> rho_l1; equal values take the original single-normalization IC path bitwise unchanged), `roche.neck_blend` (multiplier), `roche.outer_bc` ("reservoir" default | "closed"). Per-lobe normalization blends across L1 via tanh in x, half-width neck_blend*cs*sqrt(2/|Phi_xx(L1)|) (DeltaPhi=cs^2). Added `RocheBinary/diagnostics.dat` time-series (step,t,dt,drift,m_lobe1,m_bridge,m_lobe2,l1_flux,max_v) + `plot_l1_flux.py`. **CRITICAL: `outer_bc=closed` is required for this setup.** The inherited reservoir BC pins open-edge ghosts to Uinit = the *differential* profile, sustaining the imbalance -> steady-state flux (~0.003, no decay) + -2.3%/orbit mass leak. `closed` switches polar j-hi (non-seam) + bridge i-edges to reflecting walls (new `wall_ghost` host-device helper) -> exactly mass-conservative, flux decays. **Lesson: do-not-break #3 (reservoir for the tenuous overcontact atmosphere) does NOT extend to a non-equilibrium IC — use closed walls whenever the IC isn't a discrete equilibrium.** Verified q=1 cs=0.25 contrast 1.5 closed: poison 0 (1+2 ranks), step-1 drift 0, 1v2 bitwise (worst=0), 4 MultiBlock ctests pass, DEBUG 5-step clean, mass drift 1.3e-15 over 56500 steps (7 orbits). Physics: L1 flux peak 1.9e-3 (t~27) -> 3e-5 (t~249), ~60x decay (self-limiting); lobe diff 0.533 -> overshoot 0.41 (t~47, Coriolis slosh) -> quasi-steady 0.46. Full equilibration does NOT occur: the corotating-frame L1 stream + Coriolis set up a steady circulation (lobe1->bridge->lobe2->outer return) locking in a residual — the physics of real mass-transfer binaries. Contrast 4:1 is too violent (Mach~2 sloshing); 1.5 is the clean subsonic demonstrator. Published `~/stardisk-site/differential_{rho,flux,lobes,drift}.png` + site section.
**Queued setup C (counter-rotating start) — DONE 2026-07-19:** new key `roche.counter_rotate` (int, default 0). When set, InitData gives the gas solid-body counter-rotation `v = omega*(y,-x)` (inertial rest seen in the corotating frame; Mach ~3.5 at the outer edge); the else branch is the old mom=0 path, so the default IC is bitwise unchanged. Run with `outer_bc=closed` (counter-rotation is NOT a discrete equilibrium; the reservoir would pin open edges to the counter-rotating state — the setup-B lesson) and `vmax=2.0` (the transient dynamically rarefies the overcontact atmosphere toward the floor; uncapped v=mom/rho runs to 1e9 in <500 steps and collapses dt — same vacuum-jet class as semi-detached). Added mass-weighted bulk speed diagnostic `BulkSpeedMoments()` -> diagnostics.dat column `bulk_v` (printed as `bulk|v|`) + `plot_counter_braking.py`; max|v| alone is worst-cell-dominated and hides the bulk braking. **CRITICAL (well-balanced reference):** SnapshotInitGhosts copies U (counter-rotating, u!=0) into U0 by default; with U0_mom = counter-rotating the Rusanov dissipation acts on (U-U0) and *vanishes* on the strong counter-rotating shear v=omega*r -> central-scheme instability (v->1e9 in <500 steps). Fix: when counter_rotate is set, SnapshotInitGhosts zeroes U0's momentum (`Uinit.setVal(0,UMX,2,nghost)`), so the reference is the hydrostatic atmosphere AT REST — density dissipation still protects the steep profile (acts on rho-rho_eq), momentum dissipation acts on the full field (standard Rusanov on momentum) which stabilizes the shear and brakes counter-rotation on the diffusion time R^2/(smax*dx) ~ 290 ~ 8 orbits. **Lesson: the well-balanced reference is the EQUILIBRIUM atmosphere (u=0); a flowing non-equilibrium U0 turns off the dissipation on its own shear and is unstable.** Verified q=1 cs=0.25 vmax=2.0 closed: poison 0 unfilled (2 ranks), step-1 drift -8.2e-16 (exact), 4 MultiBlock ctests pass (default bitwise unchanged), DEBUG 5-step clean, 0 NaN over 100,500 steps to t=350 (~10 orbits). bulk|v| decays 0.368 -> 0.031 (~12x, to sloshing ~cs/8) by ~10 orbits; max|v| capped at 2.0 during the t<30 transient then releases, settling near cs=0.25; mass drift plateaus 0.25% (closed, mass-conservative cap); 180deg point symmetry throughout. Published `~/stardisk-site/counter_rot.mp4` (67 frames) + `counter_rot_rho.png` + `counter_braking.png` + site section. NOTE: max|v| never reaches sloshing (worst tenuous cell stays ~cs); the BULK flow is the honest corotation metric.

**Queued setup D (driven donor) — DONE 2026-07-19:** new keys `roche.drive_eps` (Real, default 0 = off) and `roche.drive_tau` (Real, default 0 -> one orbital period). When drive_eps>0, Advance adds a Newtonian-relaxation source on every donor-side cell (`xc < x_l1`): `rho += dt*((1+eps(t))*rho_eq - rho)/tau`, where `eps(t) = drive_eps*min(1,t/tau)` (linear ramp 0->drive_eps over tau) and `rho_eq = rho_l1*exp(-(Phi-Phi_L1)/cs^2)` is the original hydrostatic atmosphere. Injected gas is at rest in the corotating frame (no momentum source), so the donor envelope inflates and the pressure imbalance drives L1 mass transfer. drive_eps=0 skips the source entirely -> default bitwise unchanged (the `if (p.drive_eps > 0)` guard + Advance signature gaining a `Real t` param are the only changes). New script `plot_driven_donor.py` (10-column diagnostics.dat reader: L1 flux + eps(t) overlay, lobe masses, drift). **Physics (closed system):** with `outer_bc=closed` the accretor cannot shed mass, so the system fills until L1 pressure-equilibrates: the L1 flux ramps with eps(t) (peak ~1.7e-4 at t~80 ~ 2 orbits, lag ~ one orbit), then DECAYS to a small residual (~5e-6) sustained by the ongoing drive — a quasi-steady transfer state, NOT the sustained-flux picture (for that use `outer_bc=reservoir` so the accretor sheds through its open edge and the imbalance persists). Donor saturates at +4.6% (= eps minus L1 leak); accretor fills slowly (+1.7% over 20 orbits); total mass drift plateaus at +3.6% (drive injects only enough to replace L1 loss once donor saturated). Verified q=1 cs=0.25 drive_eps=0.05 drive_tau=35 closed: poison 0 unfilled (2 ranks), step-1 mass drift **exactly 0** (default drive_eps=0, well-balanced dissipation intact), DEBUG 5-step clean (drive-on, -ftrapv/-ffpe-trap, 0 NaN), 4 MultiBlock ctests pass (default bitwise unchanged), 0 NaN over 174,736 steps to t=710 (~20 orbits), max|v| plateaus ~0.47 (Mach ~1.9 in the L1 stream cell; bulk|v| ~0.04 subsonic), dt stable. Published `~/stardisk-site/driven_donor.mp4` (875 drift frames, 48.6s) + `driven_donor_rho.png` (3-panel t=0/final/drift) + `driven_{flux,lobes,drift}.png` + site section. **Lesson:** do-not-break #8 (truncation-level drift ~ %/orbit) is about the UNDRIVEN case; the driven setup intentionally injects mass and the drift budget is set by drive_eps * m_donor (~+5% here), NOT a bug.

### Web server (rho panel viewer, 2026-07-19)

Static site at `~/stardisk-site/` (outside the repo): `index.html` (dark page, one section per test) + `stardisk_rho.png` (StarDisk t=0/t=2pi/drift) + `roche_rho.png` (RocheBinary t=0/t=4 transfer run/drift) + `roche_evolution.mp4` (64-frame video of the perturb=0.05 run to t=4, embedded in the RocheBinary section) + `m2_05.mp4` (q=0.5 pipeline demo) + `bug_showcase.png` (bring-up failure modes vs fixed scheme, q=2 20-orbit work) + `q2_20orb.mp4` (20-orbit q=2 video, 568 frames). Served by `python3 -m http.server 8000` (nix-shell python via direnv), detached with `setsid nohup`, log `~/stardisk-site/server.log`. Reachable at:

- tailnet: http://100.67.152.108:8000 (machine `blu`)
- LAN: http://178.254.33.110:8000
- local: http://localhost:8000

Restart: `setsid nohup direnv exec /home/cernetic/amrex python3 -m http.server 8000 --directory /home/cernetic/stardisk-site --bind 0.0.0.0 </dev/null >/home/cernetic/stardisk-site/server.log 2>&1 &`. Stop: `pkill -f "http.server 8000"`.

Updating the site (no restart needed — the http.server serves new/changed files immediately):
1. Produce artifacts: rho panels via the nix-shell python recipe (`python3 plot_roche.py 0 <step>` in the test dir, writes `roche_rho.png`); showcase via `make_showcase.py`; videos via `run_roche_video.sh OUTNAME roche.key=value ...` (auto-publishes `~/stardisk-site/OUTNAME.mp4` itself).
2. `cp` other PNGs to `~/stardisk-site/`.
3. Edit `~/stardisk-site/index.html`: one `<h2>` + `<p class="sub">` + `<figure>` section per artifact (`<img>` for panels, `<video controls loop muted playsinline preload="metadata" poster=...><source ...></video>` for mp4), keep the existing dark style.
4. Verify: `curl -s -o /dev/null -w "%{http_code} %{size_download}\n" http://localhost:8000/FILE` and `http://localhost:8000/` — expect 200s with plausible sizes.

NOTE: `tailscale serve` (proper HTTPS on `blu.<tailnet>.ts.net`) is blocked — serve-config writes need root/operator and sudo is broken in non-interactive shells on this box (`/run/wrappers/bin` has no sudo). Fix once from a real terminal: `sudo tailscale set --operator=$USER`, then `tailscale serve --bg /home/cernetic/stardisk-site`.

### Vision oracle (Ollama Cloud `gemma4:31b-cloud`, 2026-07-19)

The main agent runs on `zai/glm-5.2` (text-only — `read` on an image returns `[image omitted: model does not support vision]`). To inspect any plot/frame/screenshot, query the vision oracle below instead of `read`.

**Setup (already done):** the authenticated system `ollama.service` (`/bin/ollama serve`, root) proxies cloud models to ollama.com via the device keypair in `~/.ollama`. No API key, no local weights. Model: `gemma4:31b-cloud` (vision+tools+thinking; ~10–45 s/query).

**Usage** — from the JS eval kernel:

```
import { ollamaVision } from "/home/cernetic/.config/ollama/vision.mjs";
const r = await ollamaVision("/home/cernetic/stardisk-site/roche_rho.png",
                             "Describe this figure; flag any seam artifacts.",
                             { temperature: 0.2 });
if (r.ok) console.log(r.content); else console.log(r);   // {ok,via,model,dt,tokens,content|error}
```

`ollamaVision(imagePath, prompt, opts?)` is daemon-first (`http://localhost:11434/api/chat`); it falls back to direct `https://ollama.com/api/chat` with a Bearer key from `~/.config/ollama/cloud_key` only if the daemon is down. Verified accurate on `stardisk_rho.png` (3 panels, circular star hole, inferno log ρ, RdBu drift) and `roche_rho.png` (two stars at ±2, L1 bridge mass transfer, no seams).

**CRITICAL — never pull local model tags on this 7.8 GB box.** `ollama pull qwen3.5` (6.6 GB local) OOM-killed the system `ollama.service` (`failed (oom-kill)`). Only ever use `<model>:<size>-cloud` tags (manifest-only, 342 B). If `systemctl is-active ollama` shows `failed`, restart from a real terminal (sudo is broken non-interactively): `sudo systemctl reset-failed ollama && sudo systemctl start ollama` — or drop an API key in `~/.config/ollama/cloud_key` to switch to the direct-cloud fallback.
