"""Second, wider random sweep against the independent 50-digit reference (potential-based)."""
import random, sys, csv, math
from multiprocessing import Pool
sys.path.insert(0, '.')
from independent_mpmath_ref import kappa as ref
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
def lu(a, b): return math.exp(random.uniform(math.log(a), math.log(b)))
pts = []
for i in range(200):
    v1 = lu(1e-4, 3000); v2 = v1 * lu(1, 50) if random.random() < 0.6 else lu(1e-4, 3000)
    u = lu(0.003, 300)
    if min(v1, v2) + 80 > 1400: continue   # kappa beyond float64 range is tested via ln
    pts.append((v1, v2, u))
def job(p):
    try:
        k, e = ref(*p); return (*p, str(k), float(e))
    except Exception as ex:
        return (*p, "ERR " + str(ex)[:50], 0.0)
with Pool(2) as pool:
    out = pool.map(job, pts)
with open('independent_mpmath_sweep2.csv', 'w', newline='') as f:
    w = csv.writer(f)
    for r in out: w.writerow(r)
print(len(out))
