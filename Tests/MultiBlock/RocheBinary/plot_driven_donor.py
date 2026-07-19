#!/usr/bin/env python3
"""Plot the driven-donor (setup D) diagnostics from RocheBinary/diagnostics.dat.

Setup D relaxes the donor (lobe 1) density toward (1+eps(t))*rho_eq, with
eps(t) a linear ramp 0 -> drive_eps over drive_tau. The key figures are:

  * L1 mass flux vs time, overlaid with the eps(t) ramp on a twin axis -- the
    flux should track the ramp with lag ~drive_tau and then settle to a
    sustained quasi-steady value (unlike setup B, it does NOT decay).
  * Lobe masses vs time -- the donor inflates and holds, the accretor (lobe 2)
    and bridge fill as transferred mass accumulates (closed outer walls).
  * Total mass drift vs time -- grows roughly linearly while the donor ramps,
    then settles to the sustained-drive injection rate.

Usage:
    python3 plot_driven_donor.py [path/to/diagnostics.dat] [out_prefix] \\
        [drive_eps] [drive_tau]

Reads the 10-column diagnostics file written every print_int by main.cpp:
    # step t dt mass_drift m_lobe1 m_bridge m_lobe2 l1_flux max_v bulk_v
Writes <out_prefix>_flux.png, <out_prefix>_lobes.png, <out_prefix>_drift.png
(default out_prefix: driven).
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else "RocheBinary/diagnostics.dat"
out = sys.argv[2] if len(sys.argv) > 2 else "driven"
# Optional drive overlay: eps(t) = drive_eps * min(1, t/drive_tau).
drive_eps = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
# drive_tau <= 0 -> one orbital period (the exe's default).
period = 2.0 * np.pi / (np.sqrt(2.0) / 4.0)  # q=1, sep=4 default
drive_tau = float(sys.argv[4]) if len(sys.argv) > 4 else period

d = np.loadtxt(path, comments="#")
step, t, dt, mdrift, ml1, mbr, ml2, flux, maxv, bulkv = d.T

# --- Panel 1: L1 mass flux vs time + eps(t) ramp (twin axis). ---
fig, ax = plt.subplots(figsize=(6.6, 4.2))
ax.plot(t, flux, "b-", lw=1.0, label=r"$\dot M_{L1}$")
ax.axhline(0.0, color="k", lw=0.5)
ax.set_xlabel(r"time $t$  (period $\approx %.1f$)" % period)
ax.set_ylabel(r"L1 mass flux  ($>0$: donor $\to$ accretor)", color="b")
ax.tick_params(axis="y", labelcolor="b")
ax.grid(True, alpha=0.3)
ax2 = ax.twinx()
eps_of_t = drive_eps * np.minimum(1.0, t / drive_tau)
ax2.plot(t, eps_of_t, "r--", lw=1.2, label=r"$\varepsilon(t)$")
ax2.set_ylabel(r"drive amplitude $\varepsilon(t)$", color="r")
ax2.tick_params(axis="y", labelcolor="r")
ax.set_title(r"Driven donor: L1 flux responds to $\varepsilon(t)$, then equilibrates")
fig.tight_layout()
fig.savefig(out + "_flux.png", dpi=130)
plt.close(fig)

# --- Panel 2: lobe + bridge masses vs time. Donor inflates and holds;
# accretor and bridge accumulate (closed outer walls). ---
fig, ax = plt.subplots(figsize=(6.6, 4.2))
ax.plot(t, ml1, "b-", lw=1.2, label="lobe 1 (donor)")
ax.plot(t, ml2, "g-", lw=1.2, label="lobe 2 (accretor)")
ax.plot(t, mbr, "k-", lw=1.0, label="bridge")
ax.set_xlabel(r"time $t$")
ax.set_ylabel("region mass")
ax.set_title("Driven donor: lobe masses")
ax.legend(loc="best")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(out + "_lobes.png", dpi=130)
plt.close(fig)

# --- Panel 3: total mass drift. Grows while the donor ramps up, then tracks
# the sustained-drive injection rate (closed walls: no mass leaves). ---
fig, ax = plt.subplots(figsize=(6.6, 4.2))
ax.plot(t, mdrift, "r-", lw=1.0)
ax.axhline(0.0, color="k", lw=0.5)
ax.set_xlabel(r"time $t$")
ax.set_ylabel(r"$(M - M_0) / M_0$")
ax.set_title("Driven donor: total mass (drive injects mass)")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(out + "_drift.png", dpi=130)
plt.close(fig)

print("wrote %s_flux.png, %s_lobes.png, %s_drift.png" % (out, out, out))
