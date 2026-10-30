"""Deterministic target envelopes and grid checks; no statistical uncertainty model."""
from pathlib import Path
import csv, json, importlib.util, itertools, platform, time
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parents[1]/'results'
SPEC=importlib.util.spec_from_file_location('common_benchmark',ROOT/'common_calibration/code/run.py')
m=importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(m)
LAWS=('saturation','content','constant_D')
EPS=1e-7

def endpoint(law,q,pv,parameter):
    if law!='constant_D':return m.end(law,q,pv,parameter)
    # Gap coordinate avoids unphysical pressure trials at the high-target fits.
    sol=solve_ivp(lambda x,u:[parameter/(10*q*m.CP(100-(100-pv)*np.exp(-u[0])))],
                  (0,1),[0.],method='DOP853',rtol=2e-10,atol=2e-11,max_step=.05,first_step=.001)
    assert sol.success
    return m.C(100-(100-pv)*np.exp(-sol.y[0,-1]))

def write(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def fit(target):
    """Invert each reference endpoint independently of ODE fitting."""
    if not m.C(40)<target<m.C(100):
        raise ValueError('Target must be strictly between inlet and alveolar content')
    p=brentq(lambda p:m.C(p)-target,40,100,xtol=1e-12)
    return {'saturation':-m.Q0*np.log((m.S(100)-m.S(p))/(m.S(100)-m.S(40))),
            'content':-m.Q0*np.log((m.C(100)-target)/(m.C(100)-m.C(40))),
            'constant_D':10*m.Q0*quad(lambda z:m.CP(z)/(100-z),40,p,epsabs=1e-12,epsrel=1e-12)[0]}

def independent_endpoint(pv,q,d):
    # Solve in gap coordinates to keep the near-equilibrium quadrature regular.
    # P = PA - (PA-Pv)*exp(-u), integral CP(P)/(PA-P) dP = integral CP(P(u)) du.
    exposure=d/(10*q)
    def integral(u):
        return quad(lambda v:m.CP(100-(100-pv)*np.exp(-v)),0,u,epsabs=1e-12,epsrel=1e-12)[0]
    hi=1.
    while integral(hi)<exposure:hi*=2
    u=brentq(lambda u:integral(u)-exposure,0,hi,xtol=1e-12)
    return m.C(100-(100-pv)*np.exp(-u))

def main():
    start=time.perf_counter();OUT.mkdir(parents=True,exist_ok=True)
    target=m.end('saturation',m.Q0,40,17.2)
    offsets=[-.05,-.025,-.01,0.,.01,.025,.05]
    fits={o:fit(target+o) for o in offsets}
    calibration=[]
    for offset,pars in fits.items():
        for law,parameter in pars.items():
            achieved=endpoint(law,m.Q0,40,parameter)
            calibration.append(dict(target_offset=offset,target_content=target+offset,law=law,parameter=parameter,
                parameter_units='mL/min/mmHg' if law=='constant_D' else 'L/min',
                reference_output=achieved,reference_residual=achieved-target-offset,
                parameter_ratio_to_nominal=parameter/fits[0.][law]))
    assert max(abs(r['reference_residual']) for r in calibration)<EPS
    grid=list(itertools.product(np.linspace(30,50,21),np.linspace(3,5,21)))
    predictions={(o,law):np.array([endpoint(law,q,pv,pars[law]) for pv,q in grid]) for o,pars in fits.items() for law in LAWS}
    envelope_quadrature_error=max(abs(predictions[o,'constant_D'][i]-independent_endpoint(pv,q,pars['constant_D']))
        for o,pars in fits.items() for i,(pv,q) in enumerate(grid))
    assert envelope_quadrature_error<EPS
    # Endpoint monotonicity follows from increasing fitted coefficients and the scalar comparison principle.
    assert all(np.all(predictions[o1,l] <= predictions[o2,l]+1e-10) for o1,o2 in zip(offsets,offsets[1:]) for l in LAWS)
    envelopes=[];summary=[]
    for c in [0.,.01,.025,.05]:
        for l1,l2 in itertools.combinations(LAWS,2):
            lo1,hi1=predictions[-c,l1],predictions[c,l1]
            lo2,hi2=predictions[-c,l2],predictions[c,l2]
            # Marginal envelopes are deliberately conservative for a shared unknown target.
            gap=np.maximum(lo1-hi2,lo2-hi1)
            for i,(pv,q) in enumerate(grid):
                envelopes.append(dict(calibration_halfwidth=c,Pv=pv,Q=q,PA=100,law1=l1,law2=l2,
                    lower1=lo1[i],upper1=hi1[i],lower2=lo2[i],upper2=hi2[i],signed_envelope_gap=gap[i],
                    **{f'robust_e_{e:.2f}':int(gap[i]>2*e+2*EPS) for e in [.01,.05,.10]}))
            for e in [.01,.05,.10]:
                best=int(np.argmax(gap))
                summary.append(dict(calibration_halfwidth=c,output_error=e,law1=l1,law2=l2,
                    robust_states=int(np.sum(gap>2*e+2*EPS)),states=len(grid),maximum_envelope_gap=gap[best],
                    Pv_at_max=grid[best][0],Q_at_max=grid[best][1]))
    write('calibration_target_scan.csv',calibration)
    write('calibration_envelopes.csv',envelopes)
    write('calibration_envelope_summary.csv',summary)
    original=json.loads((ROOT/'common_calibration/results/summary.json').read_text())['parameters']
    checks=[]
    for match,pars in original.items():
        for pv,q in grid:
            ordinary=m.end('constant_D',q,pv,pars['constant_D'])
            refined=m.end('constant_D',q,pv,pars['constant_D'],rtol=2e-12)
            independent=independent_endpoint(pv,q,pars['constant_D'])
            checks.append(dict(match=match,Pv=pv,Q=q,ordinary=ordinary,refined=refined,quadrature=independent,
                refinement_error=abs(ordinary-refined),quadrature_error=abs(ordinary-independent)))
    assert max(max(r['refinement_error'],r['quadrature_error']) for r in checks)<EPS
    write('full_grid_verification.csv',checks)
    resolution=[]
    for n in [21,41,81]:
        for match,pars in original.items():
            values=[]
            for pv,q in itertools.product(np.linspace(30,50,n),np.linspace(3,5,n)):
                delta=m.end('constant_D',q,pv,pars['constant_D'])-m.end('saturation',q,pv,pars['saturation'])
                values.append((pv,q,delta))
            best=max(values,key=lambda r:abs(r[2]))
            for e in [.01,.05,.1]:
                count=sum(abs(v[2])>2*e+2*EPS for v in values)
                resolution.append(dict(grid_side=n,match=match,error_bound=e,qualifying_states=count,states=n*n,
                    grid_fraction=count/(n*n),maximum_absolute_difference=abs(best[2]),signed_difference_at_max=best[2],
                    Pv_at_max=best[0],Q_at_max=best[1]))
    write('grid_resolution.csv',resolution)
    report=dict(reference_target=target,alveolar_content=m.C(100),upper_headroom=m.C(100)-target,
        calibration_halfwidths=[0,.01,.025,.05],infeasible_symmetric_halfwidth_example=.1,
        calibration_cases=len(calibration),envelope_rows=len(envelopes),verification_states=len(checks),
        maximum_refinement_error=max(r['refinement_error'] for r in checks),
        maximum_quadrature_error=max(r['quadrature_error'] for r in checks),
        envelope_quadrature_states=len(fits)*len(grid),maximum_envelope_quadrature_error=envelope_quadrature_error,
        interpretation='Deterministic bounded target and output scenarios; marginal envelope separation is sufficient and conservative for a shared target. No confidence level or clinical power.',
        all_checks_passed=True,python=platform.python_version(),elapsed_seconds=time.perf_counter()-start)
    (OUT/'summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
    print(json.dumps([r for r in summary if r['law1']=='saturation' and r['law2']=='constant_D' and r['output_error']==.05],indent=2))

if __name__=='__main__':main()
