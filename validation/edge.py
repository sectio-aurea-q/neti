"""Grenzfälle gegen Arb (gleiche Implementierung wie das Referenzgitter)."""
import sys, subprocess, time
sys.path.insert(0, "/tmp/claude-0/-home-claude/909a837b-f0d9-594e-93b9-6de46e6b8ddf/scratchpad/audit")
import flint
from kappa_sq import kappa_sq as kappa
cases = [
 ("tiny barriers, cos branch", 1e-3, 1e-3, 1.0),
 ("near q=0 from below", None, None, None),
 ("strong asymmetry", 0.05, 200.0, 8.0),
 ("very low T (u=60)", 40.0, 45.0, 60.0),
 ("very high T (u=0.01)", 0.5, 0.7, 0.01),
 ("wide high barrier", 300.0, 300.0, 5.0),
 ("v=1 symmetric, u=3", 1.0, 1.0, 3.0),
]
import math
# q = (2pi/u)^2 v1 v2 - pi^2/4 = 0  <=> v1 v2 = u^2/16 ; choose u=4, v1=v2=1-1e-9 / 1+1e-9
cases[1] = ("q just below 0", 1 - 1e-7, 1 - 1e-7, 4.0)
cases.insert(2, ("q just above 0", 1 + 1e-7, 1 + 1e-7, 4.0))
for name, v1, v2, u in cases:
    t = time.time()
    try:
        k = kappa(flint.arb(v1), flint.arb(v2), flint.arb(u))
        if k is None: raise RuntimeError("no enclosure")
    except Exception as e:
        print(name, "Arb failed", e); continue
    out = subprocess.run(['../target/release/neti', 'reduced', repr(v1), repr(v2), repr(u)], capture_output=True, text=True).stdout
    import json; r = json.loads(out)
    mid = float(k.mid()); rad = float(k.rad())
    print(f"{name:26} ref={k.mid().str(17)} ±{rad/abs(mid):.1e}  neti={r['kappa']:.15e}  rel={abs(r['kappa']/mid-1):.1e}  cos={r['cos_branch']}  ({time.time()-t:.0f}s)")
