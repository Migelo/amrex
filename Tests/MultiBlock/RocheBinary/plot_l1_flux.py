#!/usr/bin/env python3
"""Plot L1 mass-transfer flux and lobe equilibration from RocheBinary/diagnostics.dat.

The key acceptance figure for the differential-contact-depth (setup B) run is
the L1 mass flux as a function of time: it should start positive (gas flowing
from the denser lobe 1 toward lobe 2) and decay ~exponentially as the lobes
equilibrate, while the lobe masses converge.

Usage:
    python3 plot_l1_flux.py [path/to/diagnostics.dat] [out_prefix]

Reads the 9-column diagnostics file written every print_int by main.cpp:
    # step t dt mass_drift m_lobe1 m_bridge m_lobe2 l1_flux max_v
Writes <out_prefix>_flux.png, <out_prefix>_lobes.png, <out_prefix>_drift.png
(default out_prefix: l1).
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else "RocheBinary/diagnostics.dat"
out = sys.argv[2] if len(sys.argv) > 2 else "l1"

d = np.loadtxt(path, comments="#")
step, t, dt, mdrift, ml1, mbr, ml2, flux, maxv = d.T

# Orbit time for reference (q=1 default: period = 2*pi/omega, omega = sqrt(2)/4).
period = 2.0 * np.pi / (np.sqrt(2.0) / 4.0) if np.any(t > 0) else 1.0

# --- Panel 1: L1 mass flux vs time (the key figure). ---
fig, ax = plt.subplots(figsize=(6.2, 4.2))
ax.plot(t, flux, "b-", lw=1.0)
ax.axhline(0.0, color="k", lw=0.5)
ax.set_xlabel(r"time $t$  (period $\approx %.1f$)" % period)
ax.set_ylabel(r"L1 mass flux $\dot M_{L1}$  ($>0$: lobe 1 $\to$ lobe 2)")
ax.set_title("L1 mass-transfer rate")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(out + "_flux.png", dpi=130)
plt.close(fig)

# --- Panel 2: lobe mass difference (the equilibration signal). The absolute
# lobe masses barely move (transferred mass ~ 3% of a lobe), but their
# difference shows the rapid initial transfer, the Coriolis-driven rebound
# (overshoot below the eventual plateau), and the quasi-steady residual. ---
fig, ax = plt.subplots(figsize=(6.2, 4.2))
diff = ml1 - ml2
ax.plot(t, diff, "k-", lw=1.2, label="lobe1 $-$ lobe2")
ax.axhline(diff[-1], color="r", ls="--", lw=0.8,
           label=r"plateau $\approx %.3f$" % diff[-1])
ax.set_xlabel(r"time $t$")
ax.set_ylabel("lobe mass difference")
ax.set_title("Equilibration (differential depth)")
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(out + "_lobes.png", dpi=130)
plt.close(fig)

# --- Panel 3: total mass drift (conservation sanity). ---
fig, ax = plt.subplots(figsize=(6.2, 4.2))
ax.plot(t, mdrift, "r-", lw=1.0)
ax.axhline(0.0, color="k", lw=0.5)
ax.set_xlabel(r"time $t$")
ax.set_ylabel(r"$(M - M_0) / M_0$")
ax.set_title("Total mass drift")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(out + "_drift.png", dpi=130)
plt.close(fig)

print("wrote %s_flux.png, %s_lobes.png, %s_drift.png" % (out, out, out))
