#!/usr/bin/env python3
"""Plot the counter-rotation braking: bulk|v| and max|v| vs time (orbits),
plus the mass drift, from a RocheBinary diagnostics.dat.

Usage: python3 plot_counter_braking.py <diagnostics.dat> <period> <out.png>
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

path, period, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
d = np.loadtxt(path)
step, t = d[:, 0], d[:, 1]
drift, maxv, bulkv = d[:, 3], d[:, 8], d[:, 9]
orbits = t / period

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6), sharex=True,
                               gridspec_kw={"height_ratios": [2, 1]})
ax1.semilogy(orbits, bulkv, "C2", lw=2, label=r"bulk $|v|$ (mass-weighted)")
ax1.semilogy(orbits, maxv, "C3", lw=1.2, label=r"max $|v|$ (worst cell)")
ax1.axhline(0.25, color="k", ls=":", lw=1, alpha=0.6, label=r"$c_s = 0.25$")
ax1.axvline(0, color="k", lw=0.5, alpha=0.3)
ax1.set_ylabel("speed")
ax1.set_title("Counter-rotating start: circulation brakes to corotation")
ax1.legend(loc="upper right")
ax1.grid(True, which="both", alpha=0.3)
ax1.set_ylim(None, 3.0)

ax2.plot(orbits, drift * 100, "C0", lw=1.5)
ax2.set_xlabel("orbits")
ax2.set_ylabel("mass drift (%)")
ax2.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(out, dpi=130)
print("wrote", out)
