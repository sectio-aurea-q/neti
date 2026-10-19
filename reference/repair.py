"""Recompute every reference point that failed or has relative radius > 1e-10 with
kappa_sq (t = s^2 substitution). Same columns; writes eckart_reference_repair.csv."""
import csv, os, sys, time
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from kappa_sq import kappa_sq
from flint import arb
SRC=os.path.join(HERE,"eckart_reference.csv"); OUT=os.path.join(HERE,"eckart_reference_repair.csv")
done=set()
if os.path.exists(OUT):
    done={(r["alpha1"],r["alpha2"],r["ustar"]) for r in csv.DictReader(open(OUT))}
rows=[r for r in csv.DictReader(open(SRC)) if r["kappa_mid"]=="FAIL" or float(r["rel_radius"])>1e-10]
new=not os.path.exists(OUT); f=open(OUT,"a",newline=""); w=csv.DictWriter(f,list(rows[0].keys()) if rows else [])
if new and rows: w.writeheader()
for r in rows:
    if (r["alpha1"],r["alpha2"],r["ustar"]) in done: continue
    a1,a2,u=float(r["alpha1"]),float(r["alpha2"]),float(r["ustar"])
    v1=arb(a1)*u/(2*arb.pi()); v2=arb(a2)*u/(2*arb.pi())
    t=time.time(); k=kappa_sq(v1,v2,arb(u)); dt=time.time()-t
    r=dict(r); r["seconds"]=round(dt,2)
    if k is None or not k.is_finite():
        r.update(kappa_mid="FAIL",kappa_lo="",kappa_hi="",rel_radius="")
    else:
        r.update(kappa_mid=k.mid().str(25,radius=False),kappa_lo=k.lower().str(25,radius=False),
                 kappa_hi=k.upper().str(25,radius=False),rel_radius=f"{float(k.rad()/abs(k.mid())):.2e}")
    w.writerow(r); f.flush(); print(a1,a2,u,r["rel_radius"] or "FAIL",flush=True)
