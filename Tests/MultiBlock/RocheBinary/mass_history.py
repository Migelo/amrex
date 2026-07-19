"""Parse a RocheBinary stdout log and plot region masses + mass rates vs time.

Usage: python3 mass_history.py <logfile> [out.png]   (defaults: /tmp/semidet_video.log, mass_history.png)

Reads the per-print_int diagnostic lines
    Step #N, t = ..., dt = ..., mass drift = ..., m_lobe1 = ..., m_bridge = ...,
        m_lobe2 = ..., L1 flux = ..., max|v| = ...
and produces a 2-panel figure:
    top:    total mass of star 1 (donor), the bridge, and star 2 (accretor) vs t (log y)
    bottom: dM/dt of each region vs t, with the L1 mass-transfer flux for reference.
"""
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STEP_RE = re.compile(
    r"Step #\d+, t = (?P<t>[-\d.eE]+),.*?"
    r"m_lobe1 = (?P<m1>[-\d.eE]+), m_bridge = (?P<mb>[-\d.eE]+), "
    r"m_lobe2 = (?P<m2>[-\d.eE]+), L1 flux = (?P<f>[-\d.eE]+)")


def parse(path):
    t, m1, mb, m2, f = [], [], [], [], []
    with open(path) as fh:
        for line in fh:
            m = STEP_RE.search(line)
            if m:
                t.append(float(m["t"]))
                m1.append(float(m["m1"]))
                mb.append(float(m["mb"]))
                m2.append(float(m["m2"]))
                f.append(float(m["f"]))
    if not t:
        sys.exit(f"mass_history: no Step lines parsed from {path}")
    return np.array(t), np.array(m1), np.array(mb), np.array(m2), np.array(f)


def smooth(x, w):
    """Moving-average smooth, edge-padded and centered, length-preserving."""
    if w < 2 or len(x) <= w:
        return x
    k = np.ones(w) / w
    xp = np.pad(x, (w // 2, w - 1 - w // 2), mode="edge")
    return np.convolve(xp, k, mode="valid")


def main():
    log = sys.argv[1] if len(sys.argv) > 1 else "/tmp/semidet_video.log"
    out = sys.argv[2] if len(sys.argv) > 2 else "mass_history.png"
    t, m1, mb, m2, f = parse(log)

    dm1 = np.gradient(m1, t)
    dmb = np.gradient(mb, t)
    dm2 = np.gradient(m2, t)
    win = max(1, len(t) // 120)

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(9.5, 7.4), sharex=True,
                                   constrained_layout=True)
    ax0.semilogy(t, m1, label="star 1 (donor)", lw=1.7)
    ax0.semilogy(t, mb, label="bridge", lw=1.7)
    ax0.semilogy(t, m2, label="star 2 (accretor)", lw=1.7)
    ax0.set_ylabel("region mass")
    ax0.legend(loc="lower right", framealpha=0.9)
    ax0.grid(True, which="both", alpha=0.3)
    ax0.set_title("Semi-detached mass transfer (star 1 -> star 2)")

    ax1.plot(t, smooth(dm1, win), lw=1.6, label=r"$dM_1/dt$ (donor)")
    ax1.plot(t, smooth(dm2, win), lw=1.6, label=r"$dM_2/dt$ (accretor)")
    ax1.plot(t, smooth(dmb, win), lw=1.6, label=r"$dM_{\rm bridge}/dt$")
    ax1.plot(t, smooth(f, win), "k--", lw=1.2, alpha=0.7, label="L1 flux")
    ax1.axhline(0.0, color="0.5", lw=0.8)
    ax1.set_ylabel("mass rate  $dM/dt$")
    ax1.set_xlabel(r"$t$   (orbital period $\approx 35.5$)")
    ax1.legend(loc="best", ncol=2, fontsize=8.8, framealpha=0.9)
    ax1.grid(True, alpha=0.3)
    # Annotate why dM2/dt << L1 flux: the accretor retains only a fraction of
    # the mass crossing L1; the rest leaves through its open far-field edges.
    crossed = float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(t)))
    retained = 100.0 * (m2[-1] - m2[0]) / crossed if crossed > 0 else 0.0
    ax1.text(0.50, 0.82,
             f"star 2 retains {retained:.0f}% of the L1 flux\n"
             r"$-$ the rest exits via its open edges",
             transform=ax1.transAxes, fontsize=8.7, color="0.88", ha="center",
             va="center",
             bbox=dict(boxstyle="round,pad=0.35", fc="0.15", ec="0.45", alpha=0.9))

    fig.savefig(out, dpi=150)
    print(f"wrote {out}; {len(t)} samples, t in [{t[0]:.4g}, {t[-1]:.4g}], "
          f"m1[{m1[0]:.5g}->{m1[-1]:.5g}], m2[{m2[0]:.3g}->{m2[-1]:.3g}]")


if __name__ == "__main__":
    main()
