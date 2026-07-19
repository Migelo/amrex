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

Contact-binary extension of StarDisk: two fixed stars at (±2, 0), each wrapped in 4 polar-sector blocks (annulus 1<r<1.6), plus a 9th curvilinear **bridge block** (Coons patch) connecting the annuli through L1. Isothermal Euler in the corotating frame (Omega^2=G(m1+m2)/D^3): source = -rho grad(Phi_roche) - 2 Omega x (rho u). IC = exact hydrostatic isothermal atmosphere rho = rho_l1 exp(-(Phi-Phi_L1)/cs^2) at rest, L1 found by bisection on dPhi/dx (host). New test dir `Tests/MultiBlock/RocheBinary/` (main.cpp, GNUmakefile, Make.package, CMakeLists.txt, plot_roche.py, compare_ranks.py); ctest auto-registers via GLOB_RECURSE (test `MultiBlock_RocheBinary_2d`, needs AMReX_SPACEDIM=2 configure). StarDisk untouched. Build: `cd Tests/MultiBlock/RocheBinary && direnv exec . make -j2` (this VPS: max 2 cores, max 2 MPI ranks). Render: `python3 plot_roche.py 0 <step>` needs numpy+matplotlib — shell.nix python lacks them; use nix-shell with `(python3.withPackages (ps: [ps.numpy ps.matplotlib]))` over the same fetchTarball pin.

Key design points: 20 one-sided seam fills = 16 tangential (offset-only, as StarDisk) + bridge.jlo<->e1.jhi (offset-only) + **bridge.jhi<->w2.jhi with dtos sign=(-1,1)** (w2's i runs top->bottom, bridge's i bottom->top; offset=(n_phi-1, ±(n_bridge-n_r+1)); sign-flip precedent Tests/MultiBlock/Advection/main.cpp:256). Bridge mapping (i~+y, j~+x) is left-handed like the polar blocks, so the same vertex metrics code (abs shoelace, A=(-ey,ex)) serves all 9 blocks. Defaults: m1=m2=1, sep=4, a=1, r_ann=1.6 (Eggleton lobe 0.379*sep=1.52 contained; 2*r_ann<sep no overlap), cs=0.25, rho_l1=1e-2, floor 1e-8, stop_time=1.0 (~150 steps). Params under `roche.*`, incl. `perturb` (lobe-1 density boost drives L1 transfer), `poison_test` (NaN-poison ghosts + count unfilled face ghosts, must be 0).

**Verified (1+2 ranks unless noted):** t=0 bboxes/volumes match analytic (polar 1.22510 vs pi(R^2-a^2)/4=1.22522 chord error; bridge 2.4701; lobe masses equal). Poison test 0 unfilled face ghosts (corners never filled, never read by 4-face stencil — excluded). **Step-1 mass drift -2.0e-16/+8.2e-16 (1/2 ranks)** -> all 20 seams + BCs exactly conservative. Static hold 100 steps: drift -6.0e-5, max|v| 0.024 stabilizing, lobe symmetry exact, L1 flux ~1e-19. Transfer run (perturb=0.05, t=4): drift -0.08%, max|v| 0.059 bounded, lobe1 -1.3e-3 net, bridge +4e-4, L1 flux ~+1e-4 toward star 2 (lobe1->bridge->lobe2 pathway, Coriolis-deflected circulation visible in drift panel). **1v2 ranks bitwise-identical state** (compare_ranks.py cell-wise plotfile diff, worst = 0; printed mass sums differ at 1e-16, reduction order). DEBUG build + 5-step run clean.

**Bugs fixed during bring-up (lessons for future multi-block tests):**
1. Nodal off-by-one in a diagnostics loop: xyv is nodal, mfi.validbox() IS the vertex box; adding +1 read garbage (trivial, caught by bbox report).
2. **FillPhysicalBCs overwrote seam-filled ghosts**: loop order FillGhosts -> FillSeams -> FillPhysicalBCs (StarDisk idiom) is only safe if BC edges and seam edges are disjoint. e1/w2 j-hi are seams (bridge), but the inherited polar j-hi zero-gradient BC overwrote the seam-filled ghosts -> one-sided flux across the seam -> +2.95e-6 mass/step created in the bridge. Fix: per-edge seam flag in BlockGeom, BCs skip seam edges. **Lesson: when adding seams to an existing block layout, audit every BC edge for overlap.**
3. **Zero-gradient outflow runaway in a tenuous hydrostatic atmosphere**: ghost copies of the interior let the discretely-imbalanced boundary shell collapse inward and the BC chases it: +25% mass and Mach>4 tenuous gas by t=4, dt collapsing. Fix: far-field reservoir BC (ghost = initial equilibrium state stored in Uinit) -> exchange driven only by genuine interior relaxation; drift -0.08%/t=4, bounded velocities. **Lesson: zero-gradient outflow is unsafe for steep stratified atmospheres; pin the boundary to equilibrium or use a reservoir.**

### Web server (rho panel viewer, 2026-07-19)

Static site at `~/stardisk-site/` (outside the repo): `index.html` (dark page, one section per test) + `stardisk_rho.png` (StarDisk t=0/t=2pi/drift) + `roche_rho.png` (RocheBinary t=0/t=4 transfer run/drift) + `roche_evolution.mp4` (64-frame video of the perturb=0.05 run to t=4, embedded in the RocheBinary section) + `m2_05.mp4` (q=0.5 pipeline demo). Served by `python3 -m http.server 8000` (nix-shell python via direnv), detached with `setsid nohup`, log `~/stardisk-site/server.log`. Reachable at:

- tailnet: http://100.67.152.108:8000 (machine `blu`)
- LAN: http://178.254.33.110:8000
- local: http://localhost:8000

Restart: `setsid nohup direnv exec /home/cernetic/amrex python3 -m http.server 8000 --directory /home/cernetic/stardisk-site --bind 0.0.0.0 </dev/null >/home/cernetic/stardisk-site/server.log 2>&1 &`. Stop: `pkill -f "http.server 8000"`.

NOTE: `tailscale serve` (proper HTTPS on `blu.<tailnet>.ts.net`) is blocked — serve-config writes need root/operator and sudo is broken in non-interactive shells on this box (`/run/wrappers/bin` has no sudo). Fix once from a real terminal: `sudo tailscale set --operator=$USER`, then `tailscale serve --bg /home/cernetic/stardisk-site`.
