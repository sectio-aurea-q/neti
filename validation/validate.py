"""Vergleicht eckart-ref mit dem zertifizierten Arb-Referenzgitter (Kugelarithmetik)."""
import csv, subprocess, statistics, mpmath
mpmath.mp.dps = 40
A = '/tmp/claude-0/-home-claude/909a837b-f0d9-594e-93b9-6de46e6b8ddf/scratchpad/audit/'
ref = {}
for f in ('eckart_reference.csv', 'eckart_reference_repair.csv'):
    for r in csv.DictReader(open(A + f)):
        if r['kappa_mid'] in ('FAIL', '') or not r['rel_radius']: continue
        if f == 'eckart_reference.csv' and float(r['rel_radius']) > 1e-10: continue
        ref[(r['alpha1'], r['alpha2'], r['ustar'])] = r
rows = []
with open('grid_in.csv', 'w') as g:
    for (a1, a2, u), r in ref.items():
        v1 = mpmath.mpf(a1) * mpmath.mpf(u) / (2 * mpmath.pi); v2 = mpmath.mpf(a2) * mpmath.mpf(u) / (2 * mpmath.pi)
        g.write(f"{mpmath.nstr(v1, 25)},{mpmath.nstr(v2, 25)},{u}\n"); rows.append(r)
out = subprocess.run(['../target/release/eckart-ref', 'grid', 'grid_in.csv'], capture_output=True, text=True).stdout.strip().split('\n')
res = []; errs = 0
for r, line in zip(rows, out):
    p = line.split(',')
    if p[3].startswith('ERR'): errs += 1; print('ERR', r['alpha1'], r['alpha2'], r['ustar'], p[3]); continue
    k = float(p[3]); kr = mpmath.mpf(r['kappa_mid'])
    res.append((float(abs(k / kr - 1)), r['domain'], r['alpha1'], r['alpha2'], r['ustar'], r['cos_branch'], float(kr), k, float(r['rel_radius'])))
res.sort(reverse=True)
print("reference points:", len(ref), "evaluated:", len(res), "errors:", errs)
print("max rel err:", res[0][0], res[0][1:6])
for x in res[:6]: print("  ", x)
print("median rel err:", statistics.median(x[0] for x in res))
print("count > 1e-13:", sum(x[0] > 1e-13 for x in res), " > 1e-12:", sum(x[0] > 1e-12 for x in res))
for b in ('0', '1'):
    s = [x[0] for x in res if x[5] == b]; print("cos_branch", b, len(s), "max", max(s))
with open('validation_result.csv', 'w') as f:
    f.write("rel_err,domain,alpha1,alpha2,ustar,cos_branch,kappa_ref,kappa_eckart_ref,ref_rel_radius\n")
    for x in sorted(res, key=lambda x: (x[1], float(x[2]), float(x[3]), float(x[4]))): f.write(",".join(map(str, x)) + "\n")
