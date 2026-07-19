"""Render one rho panel per RocheBinary plotfile step into frame_%04d.png.

Frames are numbered sequentially (0, 1, 2, ...) in plotfile step order.
The LogNorm range is computed once from step plt0000 and reused for every
frame (fixed scale -> no flicker).

Example:
  python3 render_frames.py --base video_runs/m2_05 --outdir video_runs/m2_05/frames \
      --star-x1 -1.3333 --star-x2 2.6667 --star-a 1.0
"""
import argparse
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, TwoSlopeNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_roche

FIGSIZE = (12, 4.2)
DPI = 100


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base", required=True,
                   help="plotfile tree root (contains e1/, n1/, ..., br/)")
    p.add_argument("--outdir", required=True, help="frame output directory")
    p.add_argument("--star-x1", type=float, default=-2.0,
                   help="x position of star 1 (default -2)")
    p.add_argument("--star-x2", type=float, default=2.0,
                   help="x position of star 2 (default +2)")
    p.add_argument("--star-a", type=float, default=1.0,
                   help="star mask radius (default 1)")
    p.add_argument("--drift", action="store_true",
                   help="render (rho-rho0)/rho0 relative to plt0000 instead of "
                        "rho (for subtle-transfer runs where raw rho looks static)")
    return p.parse_args()


def available_steps(base):
    steps = []
    for d in os.listdir(os.path.join(base, "e1")):
        m = re.fullmatch(r"plt(\d+)", d)
        if m and os.path.isdir(os.path.join(base, "e1", d)):
            steps.append(int(m.group(1)))
    return sorted(steps)


def tree_root(base):
    # The exe writes a RocheBinary/ subdir in its cwd; accept either level.
    for cand in (base, os.path.join(base, "RocheBinary")):
        if os.path.isdir(os.path.join(cand, "e1")):
            return cand
    raise FileNotFoundError(f"no plotfile tree (e1/...) under {base}")


def main():
    args = parse_args()
    plot_roche.BASE = tree_root(os.path.abspath(args.base))
    BLOCKS = plot_roche.BLOCKS
    read_block = plot_roche.read_block

    steps = available_steps(plot_roche.BASE)
    assert steps and steps[0] == 0, f"expected plt0000, got {steps[:3]}"
    os.makedirs(args.outdir, exist_ok=True)

    # Color scale + field setup. Drift mode renders (rho-rho0)/rho0 (relative
    # change from plt0000) -- the informative field for subtle-transfer runs
    # where raw rho barely moves. Vacuum cells (rho0 below the floor) are
    # masked so the diverging scale isn't dominated by division noise.
    rho0 = {b: read_block(b, 0)[0] for b in BLOCKS} if args.drift else None
    if args.drift:
        floor = 1e-4
        last = steps[-1]
        allabs = []
        for b in BLOCKS:
            r0 = rho0[b]
            r = read_block(b, last)[0]
            m = r0 >= floor
            if m.any():
                allabs.append(np.abs((r[m] - r0[m]) / r0[m]))
        vmax = max(float(np.percentile(np.concatenate(allabs), 98)),
                   0.05) if allabs else 0.5
        cmap = plt.get_cmap("RdBu_r")
        cmap.set_bad("#15151f")
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
        field_name, units = r"$(\rho-\rho_0)/\rho_0$", floor
        print(f"drift mode: {len(steps)} steps; scale +/-{vmax:.3g}; "
              f"mask rho0 < {floor:.0e}; stars at x={args.star_x1:.6g}, "
              f"{args.star_x2:.6g}, a={args.star_a:.6g}")
    else:
        vmin, vmax = np.inf, -np.inf
        for b in BLOCKS:
            rho = read_block(b, 0)[0]
            vmin = min(vmin, float(rho.min()))
            vmax = max(vmax, float(rho.max()))
        cmap, norm, field_name = "inferno", LogNorm(vmin=vmin, vmax=vmax), r"$\rho$"
        print(f"{len(steps)} steps; rho scale [{vmin:.6g}, {vmax:.6g}]; "
              f"stars at x={args.star_x1:.6g}, {args.star_x2:.6g}, a={args.star_a:.6g}")

    for iframe, step in enumerate(steps):
        fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
        mesh = None
        time = None
        for b in BLOCKS:
            rho, _mx, _my, x, y, time = read_block(b, step)
            if args.drift:
                r0 = rho0[b]
                field = np.where(r0 >= units, (rho - r0) / r0, np.nan)
            else:
                field = rho
            mesh = ax.pcolormesh(x, y, field, norm=norm, cmap=cmap,
                                 shading="auto")
        for sx in (args.star_x1, args.star_x2):
            ax.add_patch(plt.Circle((sx, 0.0), args.star_a,
                                    color="0.25", zorder=3))
        ax.set_aspect("equal")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(f"{field_name}, t = {time:.6g}")
        fig.colorbar(mesh, ax=ax, fraction=0.025, pad=0.02, aspect=30)
        fig.tight_layout()
        out = os.path.join(args.outdir, f"frame_{iframe:04d}.png")
        fig.savefig(out)
        plt.close(fig)
        if iframe % 10 == 0 or iframe == len(steps) - 1:
            print(f"{out}  (step {step}, t={time:.6g})")
    print(f"wrote {len(steps)} frames to {args.outdir}/")


if __name__ == "__main__":
    main()
