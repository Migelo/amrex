"""Assemble the 9 RocheBinary blocks from AMReX plotfiles and plot rho.

Usage: python3 plot_roche.py [step0 stepf]   (defaults: 0 0)

Plotfile layout used here (HyperCLaw-V1.1):
  Header            text: ncomp, var names, dim, time, ...
  Level_0/Cell_H    text: fab boxes + FabOnDisk file/offset list
  Level_0/Cell_D_*  per fab: ascii "FAB (...) box ncomp\n" then
                    ncomp*ncells little-endian doubles, comp-major, i-fastest
"""
import re
import struct
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm

BASE = "RocheBinary"
BLOCKS = ["e1", "n1", "w1", "s1", "e2", "n2", "w2", "s2", "br"]
STAR_A = 1.0          # star radius (roche.a)
STAR_X = (-2.0, 2.0)  # star positions (roche.sep/2)


def read_header(plt):
    with open(f"{plt}/Header") as f:
        lines = [l.strip() for l in f]
    ncomp = int(lines[1])
    names = lines[2:2 + ncomp]
    dim = int(lines[2 + ncomp])
    time = float(lines[3 + ncomp])
    return names, dim, time


def read_cell_h(plt):
    with open(f"{plt}/Level_0/Cell_H") as f:
        lines = [l.strip() for l in f]
    boxes, i = [], 0
    # box list: a "(<tag>" opener line, box lines "((lo) (hi) (type))", then ")"
    while not lines[i].startswith("(("):
        i += 1
    while lines[i] != ")":
        nums = [int(v) for v in re.findall(r"-?\d+", lines[i])]
        boxes.append((nums[0:2], nums[2:4]))
        i += 1
    i += 1
    nfabs = int(lines[i])
    assert nfabs == len(boxes)
    fabs = []
    for k in range(nfabs):
        parts = lines[i + 1 + k].split()
        assert parts[0] == "FabOnDisk:"
        fabs.append((parts[1], int(parts[2])))
    return boxes, fabs


def read_fab(plt, fname, offset, lo, hi, ncomp):
    nx = hi[0] - lo[0] + 1
    ny = hi[1] - lo[1] + 1
    with open(f"{plt}/Level_0/{fname}", "rb") as f:
        f.seek(offset)
        header = f.readline().decode()
        assert header.startswith("FAB"), header[:40]
        raw = f.read(8 * ncomp * nx * ny)
    data = np.frombuffer(raw, dtype="<f8")
    assert data.size == ncomp * nx * ny
    return data.reshape(ncomp, ny, nx)


def read_block(name, step):
    plt = f"{BASE}/{name}/plt{step:04d}"
    names, dim, time = read_header(plt)
    ncomp = len(names)
    boxes, fabs = read_cell_h(plt)
    nx = max(hi[0] for lo, hi in boxes) + 1
    ny = max(hi[1] for lo, hi in boxes) + 1
    grid = np.full((ncomp, ny, nx), np.nan)
    for (lo, hi), (fname, off) in zip(boxes, fabs):
        fab = read_fab(plt, fname, off, lo, hi, ncomp)
        grid[:, lo[1]:hi[1] + 1, lo[0]:hi[0] + 1] = fab
    assert not np.isnan(grid).any(), f"holes in {plt}"
    idx = {v: k for k, v in enumerate(names)}
    return grid[idx["rho"]], grid[idx["momx"]], grid[idx["momy"]], grid[idx["x"]], grid[idx["y"]], time


def panel(ax, step, step0, kind, **kw):
    for b in BLOCKS:
        rho, mx, my, x, y, t = read_block(b, step)
        if kind == "rho":
            field = rho
        else:  # relative drift from the t=0 state
            rho0 = read_block(b, step0)[0]
            field = (rho - rho0) / rho0
        ax.pcolormesh(x, y, field, **kw)
    return t


def main():
    step0 = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    stepf = int(sys.argv[2]) if len(sys.argv) > 2 else step0
    fig, axes = plt.subplots(1, 3, figsize=(19.5, 4.6), constrained_layout=True)

    vmin, vmax = np.inf, -np.inf
    for b in BLOCKS:
        for s in (step0, stepf):
            r = read_block(b, s)[0]
            vmin, vmax = min(vmin, r.min()), max(vmax, r.max())

    t0 = panel(axes[0], step0, step0, "rho", norm=LogNorm(vmin, vmax), cmap="inferno")
    axes[0].set_title(rf"$\rho$, $t={t0:.3f}$")
    tf = panel(axes[1], stepf, step0, "rho", norm=LogNorm(vmin, vmax), cmap="inferno")
    axes[1].set_title(rf"$\rho$, $t={tf:.3f}$")
    panel(axes[2], stepf, step0, "drift",
          norm=TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=0.5), cmap="RdBu_r")
    axes[2].set_title(rf"$(\rho-\rho_0)/\rho_0$ at $t={tf:.3f}$")

    for ax, label in zip(axes, [r"$\rho$", r"$\rho$", r"$(\rho-\rho_0)/\rho_0$"]):
        ax.set_aspect("equal")
        ax.set_xlabel("x"); ax.set_ylabel("y")
        for sx in STAR_X:
            ax.add_patch(plt.Circle((sx, 0.0), STAR_A, color="0.25", zorder=3))
        fig.colorbar(ax.collections[0], ax=ax, label=label, shrink=0.85)

    out = "roche_rho.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}; t_final={tf:.10f}, rho in [{vmin:.6g}, {vmax:.6g}]")


if __name__ == "__main__":
    main()
