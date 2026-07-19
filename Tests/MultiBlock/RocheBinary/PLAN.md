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
- **No 1-vs-2-rank bitwise comparison** (user directive 2026-07-19). Run on 2
  ranks; the dual-run `compare_ranks.py` cell-diff check is dropped. The seam /
  conservation invariants are already rank-independent by construction.
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
# 2. Step-1 invariants (2 ranks):
direnv exec . mpirun -n 2 ./main2d.gnu.MPI.ex roche.poison_test=1 roche.max_steps=1 roche.print_int=1 roche.plot_int=1000000
#    expect: "Poison test: unfilled face-ghost cells = 0", mass drift ~0 (exactly 0 in equilibrium)
# 3. Static hold (no perturb): roche.stop_time=30 — drift ~1e-3 max, max|v| bounded << cs
# 4. ctest --test-dir build -R MultiBlock  (all 4 pass)
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

### B. Differential contact depth (roche.rho_l1_1 != rho_l1_2) &mdash; DONE (2026-07-19)
- IC: per-lobe normalization, smoothed over DeltaPhi ~ cs^2 at the neck (avoid
  a hard jump at L1; the bridge is shared). Self-limiting transfer: tau ~
  DeltaM/m_dot ~ 10 orbits of exponentially decaying L1 flux.
- Acceptance: L1 flux(t) decays ~exponentially; lobes equilibrate (depth
  difference halves over the run); plot L1 flux vs t as the key figure.
- RESULT (q=1, cs=0.25, rho_l1_1=1.5e-2 rho_l1_2=1e-2, contrast 1.5): new keys
  `roche.rho_l1_1` / `roche.rho_l1_2` (sentinel <0 -> rho_l1, so default IC is
  bitwise unchanged), `roche.neck_blend` (multiplier on the blend width), and
  `roche.outer_bc` ("reservoir" default | "closed"). The per-lobe normalization
  blends smoothly across L1 via a tanh in x with half-width
  neck_blend*cs*sqrt(2/|Phi_xx(L1)|) (DeltaPhi(w)=cs^2); equal values take the
  original single-normalization path exactly. Added diagnostics.dat time-series
  (step,t,dt,drift,m_lobe1,m_bridge,m_lobe2,l1_flux,max_v) + plot_l1_flux.py.
- **CRITICAL FINDING (outer_bc):** the inherited far-field reservoir BC pins
  open-edge ghosts to Uinit = the differential profile, which SUSTAINS the
  imbalance. With the reservoir the L1 flux reaches a steady ~0.003 (orbital
  modulation, no decay) and mass leaks -2.3%/orbit through the open edges.
  `roche.outer_bc=closed` switches the polar j-hi (non-seam) + bridge i-edges to
  reflecting walls (new wall_ghost helper); the closed system is exactly
  mass-conservative and the flux decays. Lesson: do-not-break #3 (reservoir for
  the tenuous overcontact atmosphere) does NOT extend to a non-equilibrium IC --
  the reservoir pins the boundary to the initial state and prevents equilibration.
  Use closed walls whenever the IC is not a discrete equilibrium.
- Verified: step-1 drift 0 (1+2 ranks), poison 0 unfilled (1+2 ranks), 1v2
  bitwise identical (closed, compare_ranks worst=0), 4 MultiBlock ctests pass,
  DEBUG 5-step clean (wall_ghost + closed branches under -fcheck=bounds
  -ffpe-trap), mass drift 1.3e-15 over 56,500 steps (7 orbits).
- Physics (contrast 1.5, closed, t=249 ~7 orbits): L1 flux peak 1.9e-3 (t~27)
  -> 3e-5 (t~249), a ~60x decay (self-limiting, the key figure). Lobe mass
  difference drops 0.533 -> overshoot 0.41 (t~47, Coriolis sloshing) -> rebounds
  to a quasi-steady ~0.46. Full convergence to zero does NOT occur: in the
  corotating frame the L1 stream + Coriolis force establish a steady circulation
  (lobe1 -> bridge -> lobe2 -> outer return) that locks in a persistent residual
  -- the physics that sustains real mass-transfer binaries rather than letting
  them instantly equalize. The "depth difference halves" acceptance is only met
  transiently (overshoot to 0.41 = 23% reduction); the honest quasi-steady
  residual is 0.46 (13% net). Contrast 4:1 is too violent (Mach~2, strong
  sloshing); contrast 1.5 is the clean subsonic demonstrator. E-folding of the
  flux envelope ~ a few orbits; the residual circulation is long-lived.
- Published: ~/stardisk-site/differential_{rho,flux,lobes,drift}.png + site section.

### C. Counter-rotating start (inertial-rest gas) &mdash; DONE (2026-07-19)
- IC: mom = rho * omega * (y, -x) instead of 0 (one line in InitData). Mach
  ~3.5 at the outer edge initially: strong shocks early, circulation brakes
  and re-corotates over ~10 orbits.
- Acceptance: max|v| peaks early then decays; by t ~ 10 orbits the mean flow
  is ~corotating (|v| ~ sloshing level, no systematic counter-rotation).
- RESULT (q=1, cs=0.25, stop_time=350 ~ 10 orbits): new key `roche.counter_rotate`
  (int, default 0). When set, InitData gives the gas solid-body counter-rotation
  v = omega*(y,-x) (inertial rest seen in the corotating frame); else branch is
  the old mom=0 path, so default IC is bitwise unchanged. Run with
  `outer_bc=closed` (counter-rotation is not a discrete equilibrium; the
  reservoir would pin the open edges to counter-rotation -- the setup-B lesson)
  and `vmax=2.0` (caps low-density vacuum jets during the violent transient; the
  overcontact min rho ~0.002 drops toward floor as the counter-rotating gas
  dynamically rarefies, and uncapped v=mom/rho runs to 1e9 in <500 steps and
  collapses dt). Added a mass-weighted bulk speed diagnostic
  (BulkSpeedMoments -> diagnostics.dat column `bulk_v`, printed as `bulk|v|`)
  because max|v| is dominated by a few tenuous jet cells and hides the bulk
  braking. New plot script `plot_counter_braking.py` (bulk/max|v| + drift vs
  orbits).
- **CRITICAL (well-balanced reference):** SnapshotInitGhosts copies U (the
  counter-rotating IC, u != 0) into U0 by default. With U0_mom = counter-
  rotating, the Rusanov dissipation acts on (U-U0) and *vanishes* wherever the
  gas counter-rotates -- the strong shear v=omega*r then has no numerical
  dissipation (a central scheme) and goes unstable: v -> 1e9 in <500 steps.
  Fix: when counter_rotate is set, SnapshotInitGhosts zeros U0's momentum
  (`Uinit.setVal(0, UMX, 2, nghost)`), so the reference is the hydrostatic
  atmosphere AT REST. Density dissipation still protects the steep profile
  (acts on rho-rho_eq); momentum dissipation acts on the full field (standard
  Rusanov on momentum), which stabilizes the shear and brakes the counter-
  rotation on the diffusion timescale R^2/(smax*dx) ~ 290 ~ 8 orbits. Lesson:
  **the well-balanced reference is the EQUILIBRIUM atmosphere (u=0); a flowing
  non-equilibrium U0 turns off the dissipation on its own shear and is unstable.**
- Verified: poison 0 unfilled (2 ranks), step-1 mass drift -8.2e-16 (exact,
  closed walls), 4 MultiBlock ctests pass (default IC bitwise unchanged), DEBUG
  5-step clean (-ftrapv/-ffpe-trap, 0 NaN), 0 NaN over 100,500 steps to t=350.
  bulk|v| decays 0.368 -> 0.031 (~12x, to sloshing level ~cs/8) by ~10 orbits;
  max|v| capped at 2.0 during the t<30 transient then releases, settling near
  cs=0.25; mass drift plateaus at 0.25% (closed, mass-conservative cap); 180deg
  point symmetry throughout (e1<->w2, n1<->s2). Vision oracle: clean rho panel,
  no seam artifacts. Published ~/stardisk-site/counter_rot.mp4 (67 frames) +
  counter_rot_rho.png + counter_braking.png + site section. NOTE: a fully clean
  "max|v| -> sloshing" is not reached (worst tenuous cell stays ~cs); the BULK
  flow is the honest corotation metric and it does brake to sloshing.

### D. Driven donor (slow envelope expansion) &mdash; DONE (2026-07-19)
- Add to Advance: relaxation of lobe-1 rho toward (1+eps(t)) * equilibrium,
  eps ramping 0 -> eps_max over tau_drive (params roche.drive_eps,
  roche.drive_tau). Mimics donor expansion; the system never settles — it
  tracks the drive through quasi-steady transfer states.
- Acceptance: L1 flux(t) tracks the ramp (lag < 1 orbit); lobe2 mass gains
  ~the integrated flux; no runaway over 20 orbits.
- RESULT (q=1, cs=0.25, drive_eps=0.05, drive_tau=35, outer_bc=closed, stop_time=710 ~ 20 orbits):
  new keys `roche.drive_eps` (Real, default 0 = off) + `roche.drive_tau` (Real, default 0 -> one
  orbital period). When drive_eps>0, Advance adds a Newtonian-relaxation source on every donor-side
  cell (xc < x_l1): `rho += dt*((1+eps(t))*rho_eq - rho)/tau`, eps(t)=drive_eps*min(1,t/tau) linear
  ramp, rho_eq=rho_l1*exp(-(Phi-Phi_L1)/cs^2). Injected gas at rest in corotating frame (no momentum
  source) -> donor envelope inflates, pressure imbalance drives L1 transfer. drive_eps=0 skips the
  source -> default bitwise unchanged (only changes: `if (p.drive_eps>0)` guard + Advance gaining a
  `Real t` param). New script `plot_driven_donor.py` (10-col diagnostics reader: L1 flux + eps(t)
  overlay, lobe masses, drift).
- **Physics (closed system):** with outer_bc=closed the accretor cannot shed mass, so the system
  fills until L1 pressure-equilibrates. The L1 flux ramps with eps(t) (peak ~1.7e-4 at t~80 ~ 2
  orbits, lag ~ one orbit -> acceptance MET), then DECAYS to a small residual (~5e-6) sustained by
  the ongoing drive — a quasi-steady transfer state, NOT a sustained-flux picture (for that use
  outer_bc=reservoir so the accretor sheds through its open edge and the imbalance persists).
  Donor saturates at +4.6% (= eps minus L1 leak); accretor fills slowly (+1.7% over 20 orbits ->
  acceptance MET: lobe2 gains monotonic, of the order of the integrated flux); total mass drift
  plateaus at +3.6% (drive injects only enough to replace L1 loss once donor saturated); no runaway
  over 20 orbits -> acceptance MET (max|v| plateaus ~0.47 Mach~1.9 in the L1 stream cell; bulk|v|
  ~0.04 subsonic; dt stable; 0 NaN over 174,736 steps).
- Verified: poison 0 unfilled (2 ranks), step-1 mass drift **exactly 0** (default drive_eps=0 ->
  well-balanced dissipation intact), DEBUG 5-step clean (drive-on, -ftrapv/-ffpe-trap, 0 NaN),
  4 MultiBlock ctests pass (default bitwise unchanged).
- Published: `~/stardisk-site/driven_donor.mp4` (875 drift frames, 48.6s) + `driven_donor_rho.png`
  (3-panel t=0/final/drift) + `driven_{flux,lobes,drift}.png` + site section.
- **Lesson:** do-not-break #8 (truncation-level drift ~ %/orbit) is about the UNDRIVEN case; the
  driven setup intentionally injects mass and the drift budget is set by drive_eps * m_donor
  (~+5% here), NOT a bug. All four queued setups (A, B, C, D) are now DONE.

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
