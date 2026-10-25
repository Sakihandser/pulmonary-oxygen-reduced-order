"""Common calibration and deterministic discrimination; no patient validation."""
from pathlib import Path
import csv, json, itertools, platform
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
import scipy
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parents[1]/'results'
OUT.mkdir(exist_ok=True,parents=True)
Q0=3.9464411557434818
def S(p): return p**2.7/(p**2.7+26.8**2.7)
def SP(p): return 2.7*S(p)*(1-S(p))/p
def C(p): return 20.1*S(p)+.0031*p
def CP(p): return 20.1*SP(p)+.0031
def end(law,q,pv,param,pa=100,rtol=2e-10):
    if law=='saturation':
        s=S(pa)-(S(pa)-S(pv))*np.exp(-param/q)
        return C(26.8*(s/(1-s))**(1/2.7))
    if law=='content': return C(pa)-(C(pa)-C(pv))*np.exp(-param/q)
    sol=solve_ivp(lambda x,y:[param/(10*q)*(pa-y[0])/CP(y[0])],(0,1),[pv],method='DOP853',rtol=rtol,atol=rtol/10)
    assert sol.success
    return C(sol.y[0,-1])
def write(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    a=.06*86/.3
    target=end('saturation',Q0,40,a)
    dc=brentq(lambda d:end('constant_D',Q0,40,d)-target,.1,40,xtol=1e-11)
    ac=-Q0*np.log((C(100)-target)/(C(100)-C(40)))
    slope=a/Q0*CP(40)/SP(40)*(S(100)-S(40))
    pars={'endpoint':{'saturation':a,'content':ac,'constant_D':dc},
          'inlet_slope':{'saturation':a,'content':slope*Q0/(C(100)-C(40)),'constant_D':slope*10*Q0/60}}
    cal=[]
    for match,p in pars.items():
        for law,v in p.items():
            sl=(v/Q0*CP(40)/SP(40)*(S(100)-S(40)) if law=='saturation' else v/Q0*(C(100)-C(40)) if law=='content' else v/(10*Q0)*60)
            cal.append(dict(match=match,law=law,parameter=v,units='L/min' if law!='constant_D' else 'mL/min/mmHg',reference_endpoint=end(law,Q0,40,v),inlet_content_slope=sl))
    assert max(abs(r['reference_endpoint']-target) for r in cal if r['match']=='endpoint')<1e-8
    assert max(abs(r['inlet_content_slope']-slope) for r in cal if r['match']=='inlet_slope')<1e-10
    write('calibration.csv',cal)
    rows=[];maxerr=0
    states=[(40,Q0,100),(30,Q0,100),(50,Q0,100),(40,3,100),(40,5,100),(40,Q0,80),(40,Q0,120)]
    for match,p in pars.items():
        for pv,q,pa in states:
            ref=end('saturation',q,pv,a,pa)
            for law,v in p.items():
                value=end(law,q,pv,v,pa)
                err=abs(value-end(law,q,pv,v,pa,rtol=2e-12));maxerr=max(maxerr,err)
                rows.append(dict(match=match,Pv=pv,Q=q,PA=pa,law=law,content=value,residual=value-ref,relative_percent=100*(value-ref)/target,numerical_refinement=err,role='calibration' if (pv,q,pa)==states[0] else 'withheld_synthetic'))
    write('transfer.csv',rows)
    # Dense common-boundary inlet/flow grid. Error bounds are illustrative, not inferred precision.
    dense=[]
    for match,p in pars.items():
        for pv,q in itertools.product(np.linspace(30,50,21),np.linspace(3,5,21)):
            vals={law:end(law,q,pv,v) for law,v in p.items()}
            for law1,law2 in itertools.combinations(p,2):
                delta=abs(vals[law1]-vals[law2])
                dense.append(dict(match=match,Pv=pv,Q=q,law1=law1,law2=law2,absolute_difference=delta,relative_percent=100*delta/target,exceeds_0_01=int(delta>2*.01+2e-7),exceeds_0_05=int(delta>2*.05+2e-7),exceeds_0_10=int(delta>2*.10+2e-7)))
    write('discrimination_grid.csv',dense)
    summary=[]
    for match in pars:
        for l1,l2 in itertools.combinations(pars[match],2):
            rr=[r for r in dense if (r['match'],r['law1'],r['law2'])==(match,l1,l2)]
            for e,key in [(.01,'exceeds_0_01'),(.05,'exceeds_0_05'),(.1,'exceeds_0_10')]:
                yes=[r for r in rr if r[key]]
                summary.append(dict(match=match,law1=l1,law2=l2,error_bound=e,distinguishable=len(yes),states=len(rr),Pv_min=min((r['Pv'] for r in yes),default=''),Pv_max=max((r['Pv'] for r in yes),default=''),Q_min=min((r['Q'] for r in yes),default=''),Q_max=max((r['Q'] for r in yes),default=''),maximum_difference=max(r['absolute_difference'] for r in rr)))
    write('discrimination_summary.csv',summary)
    hyd=list(csv.DictReader((ROOT/'sensitivity/results/hydraulic_scan.csv').open(encoding='utf-8-sig')))
    # Shares are the same normalized CT inputs used by the hydraulic model.
    import importlib.util
    spec=importlib.util.spec_from_file_location('lung_common',ROOT/'hydraulic_network/code/build_lung_model.py');lung=importlib.util.module_from_spec(spec);spec.loader.exec_module(lung)
    w=np.array(list(lung.LOBE_W.values()));lobes=list(lung.LOBE_W)
    scan=[]
    for d in np.linspace(3,30,55):
        for lam in [0,.25,.5,.75,1]:
            for r in hyd:
                f=np.array([float(r['f_'+l]) for l in lobes]);shares=(1-lam)*w+lam*f
                cm=sum(fi*end('constant_D',Q0*fi,40,d*si) for fi,si in zip(f,shares))
                upper=end('constant_D',Q0,40,d)
                assert cm<=upper+1e-7
                scan.append(dict(k=int(r['k']),LUL_fraction=f[lobes.index('LUL')],D=d,lambda_value=lam,content=cm,allocation_deficit=upper-cm))
    base={(r['D'],r['lambda_value']):r['content'] for r in scan if r['k']==1}
    for r in scan:r['redistribution_loss']=base[r['D'],r['lambda_value']]-r['content']
    write('allocation_scan.csv',scan)
    quaderr=[]
    for pv,q in [(30,Q0),(40,Q0),(40,5),(50,3)]:
        p=brentq(lambda p:quad(lambda z:CP(z)/(100-z),pv,p,epsabs=1e-12)[0]-dc/(10*q),pv,99.999999,xtol=1e-11)
        quaderr.append(abs(C(p)-end('constant_D',q,pv,dc)))
    assert maxerr<1e-7 and max(quaderr)<1e-7
    assert max(abs(r['allocation_deficit']) for r in scan if r['lambda_value']==1)<1e-7
    result=dict(reference_content=target,parameters=pars,refinement_max=maxerr,independent_quadrature_max=max(quaderr),calibration_rows=len(cal),transfer_rows=len(rows),grid_pair_rows=len(dense),allocation_rows=len(scan),python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,checks_passed=True,observation_error='Illustrative bounded error per endpoint; intervals disjoint only when delta > 2*e + 2e-7; no probability or clinical threshold.')
    (OUT/'summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
