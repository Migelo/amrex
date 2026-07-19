"""Cell-wise compare of two RocheBinary plotfile trees (fab-layout agnostic)."""
import sys
import numpy as np
sys.path.insert(0, ".")
from plot_roche import read_block, BLOCKS
import plot_roche

step = int(sys.argv[3]) if len(sys.argv) > 3 else 10
worst = 0.0
for b in BLOCKS:
    plot_roche.BASE = sys.argv[1]
    r1 = read_block(b, step)
    plot_roche.BASE = sys.argv[2]
    r2 = read_block(b, step)
    for name, a, c in zip(("rho", "momx", "momy", "x", "y"), r1[:5], r2[:5]):
        d = np.abs(a - c).max()
        worst = max(worst, d)
        if d != 0.0:
            print(f"{b}/{name}: max abs diff {d:.3e}")
print(f"step {step}: worst cell diff = {worst:.3e}",
      "BITWISE IDENTICAL" if worst == 0.0 else "DIFFERS")
