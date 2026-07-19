# PLAN.md — RocheBinary multi-block test: handoff instructions

You are the next agent on the `multiblock-three-block-advection` branch of the
AMReX repo. This document is everything you need to continue: context, rules,
code map, the "do not break" list, verification toolkit, publishing workflow,
and the queued future setups. Read `AGENTS.md` (repo root) first — its
"Current Situation" section is the live scratchpad; keep it updated as you work.

## Mission so far

`Tests/MultiBlock/RocheBinary/` simulates a contact binary (two fixed stars,
gas bridge through L1) on a 9-block curvilinear mesh: 2×4 polar-sector blocks
(annuli around each star) + 1 Coons-patch bridge block through the L1 corridor.
Isothermal Euler in the corotating frame, first-order FV with a well-balanced
Rusanov flux, 20 one-sided NonLocalBC seam fills (no periodicity anywhere).
Verified: exact seam conservation (step-1 drift 0/1e-16), 1v2-rank bitwise
identical, poison test clean, ctest `MultiBlock_RocheBinary_2d` passes, and a
20-orbit q=2 run (226,420 steps, t=580.4) survives with bounded velocities and
sustained L1 mass transfer. Videos/panels: http://100.67.152.108:8000 (site at
`~/stardisk-site/`, see AGENTS.md "Web server" section).

## Hard rules (this box)

- Weak VPS: **max 2 build cores (`make -j2`), max 2 MPI ranks**. No exceptions.
- No system compiler/python with numpy. Everything runs through
  `direnv exec . CMD` (nix env: gcc, openmpi, cmake, plain python3).
- numpy+matplotlib python ONLY via:
  `nix-shell -E 'with import (builtins.fetchTarball { url = "https://github.com/NixOS/nixpkgs/archive/nixos-unstable.tar.gz"; }) {}; mkShell { packages = [ (python3.withPackages (ps: [ ps.numpy ps.matplotlib ])) ]; }' --run "python3 SCRIPT ARGS"`
- ffmpeg/ffprobe are on PATH (system-wide, no wrapper).
- Build: `cd Tests/MultiBlock/RocheBinary && direnv exec . make -j2`
  (`DEBUG=TRUE` for `main2d.gnu.DEBUG.MPI.ex`, ~4.5 min).
- Run: `direnv exec . mpirun -n 2 ./main2d.gnu.MPI.ex roche.key=value ...`
- ctest (configure once): `direnv exec . cmake -S . -B build -DAMReX_SPACEDIM=2
  -DAMReX_ENABLE_TESTS=ON -DAMReX_TEST_TYPE=All` then
  `direnv exec . ctest --test-dir build -R MultiBlock --output-on-failure`.
  Re-run cmake configure after adding any CMakeLists (GLOB_RECURSE at configure time).

## Code map (`Tests/MultiBlock/RocheBinary/`)

- `main.cpp` (single file, ~1150 lines):
  - `RocheParams`: all `roche.*` ParmParse keys + derived (omega, x1, x2, x_l1, phi_l1).
  - `derive()`: omega/COM/L1 (bisection on dPhi/dx) + asserts: r_ann1+r_ann2 < sep,
    L1 inside bridge corridor (x1+r_ann1 < x_L1 < x2-r_ann2).
  - `potential(x,y,p)`: Roche potential (-m1/r1 - m2/r2 - 0.5 omega^2 r^2).
  - `BlockGeom {kind, cx, cy, theta_b, r_ann, jhi_seam}`; `vertex_pos(i,j,g,p)`:
    polar sector mapping (per-block r_ann) + Coons bridge (arcs + straight top/bottom).
  - `rusanov_flux(UL,UR,U0L,U0R,ax,ay,cs,F)`: **well-balanced** — dissipation on
    U-U0, not U. Centered pressure flux unchanged.
  - `RocheBlock` (AmrCore): InitData (vertices, shoelace volumes, hydrostatic IC
    rho=rho_l1*exp(-(Phi-phi_l1)/cs^2), mom=0, Uinit interior copy),
    FillPhysicalBCs (wall at star j-lo; reservoir from Uinit on open edges;
    SKIPS seam edges), Advance (flux + gravity/centrifugal + Coriolis),
    FillGhosts (intra-block FillBoundary), ComputeDt, RhoMinMax, NaNCount,
    TotalMass, TotalVolume, VertexBBox, MaxVel, MidplaneFlux (bridge j-midplane),
    PoisonGhosts/UnfilledGhostCount (poison test), SnapshotInitGhosts.
  - Seam setup in MyMain: 2×8 tangential (offset-only) + bridge jlo<->e1 jhi
    (offset) + bridge jhi<->w2 jhi (dtos sign=(-1,1)). Ghost-slab helpers
    seam_ghost_{lo,hi} (i-edges), seam_ghost_j{lo,hi} (j-edges).
  - Diagnostics per print_int: total mass drift, per-region masses
    (lobe1/bridge/lobe2), L1 flux (MidplaneFlux), max|v|, per-block mass/rho/NaN.
    Final amrex::Abort on any NaN (ctest pass/fail).
- `plot_roche.py`: plotfile reader (`read_block(name,step) -> rho,momx,momy,x,y,t`,
  `BLOCKS` = 9 block names) + 3-panel rho/drift figure.
- `render_frames.py`: argparse single-panel rho frames (star masks from
  --star-x1/x2/--star-a), fixed LogNorm from frame 0.
- `run_roche_video.sh OUTNAME [roche.key=value ...]`: full pipeline (isolated
  `video_runs/OUTNAME/` workdir, frames, ffmpeg h264/yuv420p, publishes to
  `~/stardisk-site/OUTNAME.mp4`). Env knobs FPS, CRF.
- `compare_ranks.py A B step`: cell-wise plotfile diff (1v2-rank bitwise check).
- `make_showcase.py`: bug-showcase figure (data transcribed from run logs).

## Do-not-break list (lessons paid for with bugs)

1. Fill order every step: `FillGhosts()` (intra-block) **then** `FillSeams()`
   then `FillPhysicalBCs()`. The NonLocalBC machinery reads src ghost rings
   (ngrow-grown slivers); they must hold current data.
2. **BCs must skip seam edges.** FillPhysicalBCs runs after the seams; a BC
   edge that overwrites a seam-filled ghost breaks flux cancellation and
   invents mass (bridge gained +3e-6/step). Seam edges are flagged in BlockGeom.
3. **U0's ghosts come only from SnapshotInitGhosts** (copied from U after the
   pipeline ran at t=0). Never fill them analytically: the equilibrium profile
   continues inside the star holes, wall/reservoir ghosts are mirrors/copies of
   the interior — any other U0 ghost makes the well-balanced dissipation a
   boundary mass pump (−82% mass in one step when this was wrong).
4. **cs must be tuned to the potential depth**: rho_surf/rho_l1 =
   exp(DeltaPhi/cs^2), DeltaPhi = Phi_L1 - Phi(surface). Keep the contrast
   ~e^6 (~400). q=1: cs=0.25. q=2: cs=0.38. Wrong cs = instant blow-up.
5. **r_annK must contain each star's lobe** (Eggleton: q=1 → 0.379·sep,
   q=2 → 0.44·sep ≈ 1.76, far side ~1.87) while r_ann1+r_ann2 < sep and L1
   stays in the corridor. derive() asserts both.
6. Corner ghost cells are never filled and never read (4-face first-order
   stencil); the poison test counts only face-adjacent slabs.
7. All mappings are left-handed: shoelace needs abs(), face areas A=(-ey,ex).
   Keep new blocks in the same convention.
8. First-order Rusanov + explicit gravity has truncation-level imbalance;
   expect bounded sloshing (max|v| << cs), small secular drift (~%/orbit).
   The well-balanced flux killed the diffusion runaway — do not reintroduce
   plain-U dissipation.

## Verification toolkit (run after any code change)

```
# 1. Build both: opt + DEBUG=TRUE
# 2. Step-1 invariants (1 and 2 ranks):
direnv exec . mpirun -n 1 ./main2d.gnu.MPI.ex roche.poison_test=1 roche.max_steps=1 roche.print_int=1 roche.plot_int=1000000
#    expect: "Poison test: unfilled face-ghost cells = 0", mass drift ~0 (exactly 0 in equilibrium)
# 3. Static hold (no perturb): roche.stop_time=30 — drift ~1e-3 max, max|v| bounded << cs
# 4. 1v2 ranks bitwise: run both with plot_int, then
#    python3 compare_ranks.py dirA dirB step  -> "worst cell diff = 0"
# 5. ctest --test-dir build -R MultiBlock  (all 4 pass)
```

## Publishing (site at ~/stardisk-site, served on :8000, no restart needed)

1. Panels: nix-shell python `python3 plot_roche.py 0 <step>`; showcase: `make_showcase.py`; cp PNGs to `~/stardisk-site/`.
2. Video: `./run_roche_video.sh OUTNAME roche.key=value ...` (auto-publishes mp4).
3. Edit `~/stardisk-site/index.html`: one <h2>/<p class="sub">/<figure> section
   per artifact (dark style; <img> or <video controls loop muted playsinline>).
4. Verify: `curl -s -o /dev/null -w "%{http_code} %{size_download}\n" http://localhost:8000/FILE` → 200.

## Queued setups (long settling times) — pick in this order

Context: the q=2 20-orbit run settles quasi-steady by t≈50 (~1.7 orbits);
the transient is the interesting part, so stretch it. Each setup is <30 lines
in main.cpp unless noted; run 20 orbits (~12-15 min wall) via run_roche_video.sh.

### A. Semi-detached start (lobe 2 initially empty) — DONE (2026-07-19)
- IC: rho = equilibrium profile for x < x_l1 only, floor beyond (the L1-plane
  discontinuity IS the initial condition — the stream's birth). Gate behind a
  new key, e.g. `roche.ic_mode="semi_detached"` (default "overcontact").
- Physics: stream crosses the corridor, Coriolis-deflects, wraps star 2,
  builds an envelope/disk from nothing. Settling ~10-15 orbits
  (envelope build-up M_env/(rho_l1·cs·A_neck) ~ 7 orbits + circulation).
- Acceptance: lobe2 mass grows monotonically from ~floor to ~1-5% of lobe1;
  L1 flux > 0 sustained and decaying as lobe2 fills; no NaN over 20 orbits;
  video shows deflected stream (not a symmetric bridge).
- RESULT (q=1, cs=0.25, stop_time=70 ≈ 2 orbits, per user "run until settled"): new keys
  `roche.ic_mode` ("overcontact"|"semi_detached") + `roche.vmax` (velocity cap, 0=off).
  IC cut `xc >= x_l1` -> floor. The L1-plane density jump (rho_l1 next to 1e-8 floor) drives
  an expansion-into-vacuum Riemann problem whose floor-cells develop v=mom/rho ~ 5000, collapsing
  dt (~6e-7) and stalling the run. Fix: velocity cap (roche.vmax=1.0, 4*cs) clips the vacuum jets
  after the hydro update while preserving direction; transonic stream untouched. Default vmax=0 so
  overcontact is bitwise unchanged (step-1 drift exactly 0 with vmax=1.0 too — cap is mass-conservative,
  only touches momentum, never triggers in equilibrium). Verified: poison 0 (1+2 ranks), step-1 drift 0,
  1v2 bitwise identical (compare_ranks worst=0), 4 MultiBlock ctests pass, DEBUG 5-step clean, 0 NaN to t=70.
  lobe2 grows monotonically 4.9e-8 -> 7.6e-3 (0.7% of lobe1); L1 flux sustained ~5e-4; drift -2.2% at
  t=70 (gas outflow through accretor's open reservoir edges, ~1%/orbit, do-not-break #8). Vision oracle
  confirms Coriolis-deflected stream, no seam artifacts. Published: ~/stardisk-site/semidetached.mp4 +
  semidetached_rho.png + site section. NOTE: full 1-5% lobe2 fill + L1-flux decay need ~5-10 orbits
  (run was cut to 2 orbits per user); the 20-orbit "settling" target is untested but the scheme is stable.

### B. Differential contact depth (roche.rho_l1_1 != rho_l1_2)
- IC: per-lobe normalization, smoothed over DeltaPhi ~ cs^2 at the neck (avoid
  a hard jump at L1; the bridge is shared). Self-limiting transfer: tau ~
  DeltaM/m_dot ~ 10 orbits of exponentially decaying L1 flux.
- Acceptance: L1 flux(t) decays ~exponentially; lobes equilibrate (depth
  difference halves over the run); plot L1 flux vs t as the key figure.

### C. Counter-rotating start (inertial-rest gas)
- IC: mom = rho * omega * (y, -x) instead of 0 (one line). Mach ~2 at the
  outer edge initially: strong shocks early (dt dips, cost ~2x), circulation
  brakes and re-corotates over ~10 orbits.
- Acceptance: max|v| peaks early then decays; by t ~ 10 orbits the mean flow
  is ~corotating (|v| ~ sloshing level, no systematic counter-rotation).

### D. Driven donor (slow envelope expansion) — needs a source term (~15 lines)
- Add to Advance: relaxation of lobe-1 rho toward (1+eps(t)) * equilibrium,
  eps ramping 0 -> eps_max over tau_drive (params roche.drive_eps,
  roche.drive_tau). Mimics donor expansion; the system never settles — it
  tracks the drive through quasi-steady transfer states.
- Acceptance: L1 flux(t) tracks the ramp (lag < 1 orbit); lobe2 mass gains
  ~the integrated flux; no runaway over 20 orbits.

### Cross-cutting for all four
- Use the q=1 default geometry unless the setup says otherwise; re-verify with
  the toolkit above before the long run; produce the video with
  run_roche_video.sh and a site section; if a run goes wrong, add the failure
  to make_showcase.py + the site showcase section (established practice).
- Commit code + scripts + AGENTS.md scratchpad updates on this branch;
  artifacts (plotfiles, frames, videos) stay untracked.

## Current repo state (2026-07-19)

Branch `multiblock-three-block-advection`, ahead of origin. Last commits:
`1ef2e2efa` (site inventory), `38b822247` (per-star r_ann + well-balanced
Rusanov), `887b507e0` (video pipeline), `49fcdcd10` (RocheBinary test).
All MultiBlock ctests pass. The 20-orbit q=2 run (m1=2 m2=1 cs=0.38
r_ann1=1.85 r_ann2=1.3) is the current production reference.
