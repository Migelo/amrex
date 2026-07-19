"""RocheBinary bring-up bug showcase: mass drift and max|v| vs time for the
failure modes hit during the q=2 (20-orbit) bring-up, vs the fixed scheme.

Data: logged Step lines of the respective runs (transcribed from the run
logs; the plain-Rusanov q=2 run log was overwritten by the fixed run).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "bug_showcase.png"

# (t, mass drift, max|v|) series
runaway = [  # plain Rusanov dissipation, q=2 cs=0.38 (diffusion runaway)
    (5.0, -2.3e-3, 8.7e-2), (10.0, -9.0e-3, 1.67e-1), (20.0, -4.8e-2, 3.54e-1),
    (40.0, -2.83e-1, 6.76e-1), (113.0, -8.08e-1, 8.53e-1),
]
wallbug = [  # well-balanced v1: analytic U0 ghosts inconsistent with BCs
    (3.6, -8.15e-1, 6.6), (3.73, -8.45e-1, 1.55e1), (4.58, -9.82e-1, 2.14e3),
]
fixed = [  # well-balanced + t=0 ghost snapshot, q=2 cs=0.38 r_ann1=1.85 r_ann2=1.3
    (20.6, -8.6e-4, 3.85e-2), (30.8, -7.7e-3, 5.90e-2),
    (41.0, -1.81e-2, 7.13e-2), (51.2, -2.98e-2, 7.68e-2),
]

fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.6), constrained_layout=True)

ax = axes[0]
for data, label, color, marker in [
        (runaway, "plain Rusanov: diffusion runaway", "tab:red", "o"),
        (wallbug, "well-balanced v1: BC-ghost leak", "tab:orange", "s"),
        (fixed, "well-balanced + ghost snapshot (fixed)", "tab:green", "^")]:
    t = [d[0] for d in data]
    d = [d[1] for d in data]
    ax.plot(t, d, marker=marker, color=color, label=label)
ax.set_xlabel("t"); ax.set_ylabel("total mass drift")
ax.set_title("Mass drift vs time")
ax.legend(fontsize=8, loc="lower left")
ax.grid(alpha=0.3)
ax.plot([84], [-0.54], marker="D", color="0.25", zorder=5)
ax.annotate("cs=0.25 mistune (rho contrast ~1e6):\n-54% by t=84, before any fix",
            xy=(84, -0.54), xytext=(40, -0.75),
            arrowprops=dict(arrowstyle="->", color="0.3", lw=0.8), fontsize=8,
            color="0.25")

ax = axes[1]
for data, label, color, marker in [
        (runaway, "plain Rusanov: v doubles every ~10 t", "tab:red", "o"),
        (wallbug, "well-balanced v1: v ~ 2000 (dead)", "tab:orange", "s"),
        (fixed, "fixed: bounded sloshing", "tab:green", "^")]:
    t = [d[0] for d in data]
    v = [d[2] for d in data]
    ax.semilogy(t, v, marker=marker, color=color, label=label)
ax.set_xlabel("t"); ax.set_ylabel("max |v|")
ax.set_title("max|v| vs time (cs = 0.38)")
ax.legend(fontsize=8, loc="upper left")
ax.grid(alpha=0.3, which="both")

fig.suptitle("RocheBinary bring-up bug showcase (q = m1/m2 = 2, 20-orbit target)",
             fontsize=11)
fig.savefig(OUT, dpi=150)
print(f"wrote {OUT}")
