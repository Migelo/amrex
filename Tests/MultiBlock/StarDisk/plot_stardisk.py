"""Assemble the 4 StarDisk blocks from AMReX plotfiles and plot rho.

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

BASE = "StarDisk"
BLOCKS = ["e", "n", "w", "s"]


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
    # box list: "(2 0", box lines..., ")"
    while not lines[i].startswith("(2"):
        i += 1
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
    # comp-major, i fastest -> [comp, j, i]
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
    step0, stepf = 0, 1116
    rho_f, *_ , t_f = None, None, None, None, None, None
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2), constrained_layout=True)

    vmin, vmax = np.inf, -np.inf
    for b in BLOCKS:
        for s in (step0, stepf):
            r = read_block(b, s)[0]
            vmin, vmax = min(vmin, r.min()), max(vmax, r.max())

    t0 = panel(axes[0], step0, step0, "rho", norm=LogNorm(vmin, vmax), cmap="inferno")
    axes[0].set_title(rf"$\rho$, $t={t0:.3f}$")
    tf = panel(axes[1], stepf, step0, "rho", norm=LogNorm(vmin, vmax), cmap="inferno")
    axes[1].set_title(rf"$\rho$, $t={tf:.3f}\;(2\pi)$")
    panel(axes[2], stepf, step0, "drift",
          norm=TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=0.5), cmap="RdBu_r")
    axes[2].set_title(r"$(\rho-\rho_0)/\rho_0$ at $t=2\pi$")

    for ax, label in zip(axes, [r"$\rho$", r"$\rho$", r"$(\rho-\rho_0)/\rho_0$"]):
        ax.set_aspect("equal")
        ax.set_xlabel("x"); ax.set_ylabel("y")
        # mask the star disk r < 1
        ax.add_patch(plt.Circle((0, 0), 1, color="0.25", zorder=3))
        fig.colorbar(ax.collections[0], ax=ax, label=label, shrink=0.85)

    out = "stardisk_rho.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}; t_final={tf:.10f}, rho in [{vmin:.6g}, {vmax:.6g}]")


if __name__ == "__main__":
    main()
