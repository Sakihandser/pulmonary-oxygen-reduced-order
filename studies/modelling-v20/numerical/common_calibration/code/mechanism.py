"""Controlled dissolved-oxygen ablation, with one calibration per law and condition."""
from pathlib import Path
import csv,json,importlib.util
import numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parents[1]/'results'
def write(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    spec=importlib.util.spec_from_file_location('common',ROOT/'common_calibration/code/run.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    OUT.mkdir(parents=True,exist_ok=True);rows=[];summ=[];cal=[]
    q0=m.Q0;a=17.2
    for b in [0.,.0031]:
        def c(p):return 20.1*m.S(p)+b*p
        def sat(pv,q,pa=100):
            s=m.S(pa)-(m.S(pa)-m.S(pv))*np.exp(-a/q)
            return c(26.8*(s/(1-s))**(1/2.7))
        target=sat(40,q0)
        for match in ['endpoint','inlet_slope']:
            ac=(-q0*np.log((c(100)-target)/(c(100)-c(40))) if match=='endpoint' else a*(20.1+b/m.SP(40))*(m.S(100)-m.S(40))/(c(100)-c(40)))
            cal.append(dict(dissolved_coefficient=b,match=match,a_saturation=a,a_content=ac,reference_saturation=target,reference_content=c(100)-(c(100)-c(40))*np.exp(-ac/q0)))
            group=[]
            for pv in np.arange(30,51):
                for q in np.linspace(3,5,21):
                    ys=sat(pv,q);yc=c(100)-(c(100)-c(pv))*np.exp(-ac/q)
                    r=dict(dissolved_coefficient=b,match=match,Pv=int(pv),Q=float(q),PA=100,saturation=ys,content=yc,signed_difference=yc-ys,absolute_difference=abs(yc-ys))
                    rows.append(r);group.append(r)
            errors=np.array([r['absolute_difference'] for r in group]);worst=max(group,key=lambda r:r['absolute_difference'])
            summ.append(dict(dissolved_coefficient=b,match=match,states=len(group),maximum=float(errors.max()),median=float(np.median(errors)),p95=float(np.quantile(errors,.95)),Pv_at_max=worst['Pv'],Q_at_max=worst['Q'],a_content=ac))
            if b==0:assert errors.max()<5e-13 and abs(ac-a)<1e-11
    write('dissolved_grid.csv',rows);write('dissolved_summary.csv',summ);write('dissolved_calibration.csv',cal)
    # Signed state maps retain both operators, not just an absolute residual.
    pars=json.loads((ROOT/'common_calibration/results/summary.json').read_text())['parameters'];maps=[];candidates=[]
    for match,p in pars.items():
        group=[]
        for pv in np.arange(30,51):
            for q in np.linspace(3,5,21):
                ys=m.end('saturation',q,pv,p['saturation']);yd=m.end('constant_D',q,pv,p['constant_D']);delta=yd-ys
                row=dict(match=match,Pv=int(pv),Q=float(q),PA=100,saturation=ys,constant_D=yd,signed_difference=delta,absolute_difference=abs(delta),**{f'exceeds_{e:.2f}':int(abs(delta)>2*e+2e-7) for e in [.01,.05,.10]})
                maps.append(row);group.append(row)
        # Candidate ranking is deterministic over the supplied grid, not a clinical optimum.
        for rank,r in enumerate(sorted(group,key=lambda r:r['absolute_difference'],reverse=True)[:5],1):candidates.append(dict(rank=rank,**r))
    write('signed_state_map.csv',maps);write('ranked_grid_states.csv',candidates)
    # Explicit equality at equal exposure with b=0, and exact dissolved-term identity at b>0.
    identity=0.
    for pv,q in [(30,3),(40,q0),(50,5)]:
        rr=np.exp(-a/q);ss=m.S(100)-(m.S(100)-m.S(pv))*rr;ps=26.8*(ss/(1-ss))**(1/2.7)
        cs=20.1*ss+.0031*ps;cc=m.C(100)-(m.C(100)-m.C(pv))*rr
        identity=max(identity,abs((cs-cc)-.0031*(ps-(1-rr)*100-rr*pv)))
    assert identity<1e-12
    (OUT/'checks.json').write_text(json.dumps(dict(zero_dissolved_identity_passed=True,identity_max_residual=identity,grid_rows=len(rows),state_map_rows=len(maps),all_checks_passed=True),indent=2))
    print(json.dumps(summ,indent=2))
if __name__=='__main__':main()
