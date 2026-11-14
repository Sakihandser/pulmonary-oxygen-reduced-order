"""Offline observed-endpoint validation. All fits exclude test outcomes."""
from pathlib import Path
import csv, json, hashlib, math, platform, importlib.util
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'; OUT.mkdir(parents=True,exist_ok=True)
CONFIG=json.loads((ROOT/'model_config.json').read_text(encoding='utf-8'))
BLOOD=CONFIG['blood']; P50=BLOOD['p50_mmHg']; HILL=BLOOD['hill_exponent']
BINDING=BLOOD['binding_ml_g']; DISSOLVED=BLOOD['dissolved_ml_dl_mmHg']
LAWS=['saturation','content','constant_D','pressure','equilibrium_admixture']

def write(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def read_data(p50=P50,n=HILL,rq=CONFIG['alveolar']['respiratory_quotient'],missing_hb=15.):
    with (ROOT/'data/observations.csv').open(encoding='utf-8') as f:rows=list(csv.DictReader(f))
    for r in rows:
        for k in ['flow','Pv_mmHg','Pa_mmHg','PaCO2_mmHg','Hb_g_dl','FiO2','PB_mmHg','PH2O_mmHg','baseline']:
            r[k]=float(r[k])
        r['p50']=p50;r['hill']=n
        if r['Hb_status'].startswith('assumed'):r['Hb_g_dl']=missing_hb
        r['PA_mmHg']=(float(r['PA_mmHg']) if r['PA_mmHg'] else
            r['FiO2']*(r['PB_mmHg']-r['PH2O_mmHg'])-r['PaCO2_mmHg']*(r['FiO2']+(1-r['FiO2'])/rq))
    return rows

def sat(p,p50=P50,n=HILL):return p**n/(p**n+p50**n)
def invsat(s,p50=P50,n=HILL):return p50*(s/(1-s))**(1/n)
def content(p,hb=15,p50=P50,n=HILL):return BINDING*hb*sat(p,p50,n)+DISSOLVED*p
def cp(p,hb,p50,n):
    if p<=1e-12:return DISSOLVED
    s=sat(p,p50,n);return BINDING*hb*n*s*(1-s)/p+DISSOLVED
def invcontent(c,hb,p50,n):
    if c==0:return 0.
    return brentq(lambda p:content(p,hb,p50,n)-c,1e-10,max(1000,c/DISSOLVED+1),xtol=1e-10)

def conductance_endpoint(q,pv,pa,d,hb,p50,n):
    if d==0:return pv
    if pa-pv<1e-10:return pa
    target=d/(10*q)
    def integral(z):
        return quad(lambda u:cp(pa-(pa-pv)*math.exp(-u),hb,p50,n),0,z,epsabs=1e-10,epsrel=1e-10)[0]
    if target>=integral(32):return pa
    z=brentq(lambda z:integral(z)-target,0,32,xtol=1e-11)
    return pa-(pa-pv)*math.exp(-z)

def predict(law,param,r):
    q,pv,pa,hb,p50,n=[r[k] for k in ['flow','Pv_mmHg','PA_mmHg','Hb_g_dl','p50','hill']]
    if law=='saturation':return invsat(sat(pa,p50,n)-(sat(pa,p50,n)-sat(pv,p50,n))*math.exp(-param/q),p50,n)
    if law=='content':return invcontent(content(pa,hb,p50,n)-(content(pa,hb,p50,n)-content(pv,hb,p50,n))*math.exp(-param/q),hb,p50,n)
    if law=='pressure':return pa-(pa-pv)*math.exp(-param/q)
    if law=='equilibrium_admixture':return invcontent((1-param)*content(pa,hb,p50,n)+param*content(pv,hb,p50,n),hb,p50,n)
    return conductance_endpoint(q,pv,pa,param,hb,p50,n)

def fit_one(law,r,target=None):
    y=r['Pa_mmHg'] if target is None else target
    q,pv,pa,hb,p50,n=[r[k] for k in ['flow','Pv_mmHg','PA_mmHg','Hb_g_dl','p50','hill']]
    if not pv<y<pa:raise ValueError('Calibration target outside inlet/equilibrium bounds')
    if law=='saturation':return -q*math.log((sat(pa,p50,n)-sat(y,p50,n))/(sat(pa,p50,n)-sat(pv,p50,n)))
    if law=='content':return -q*math.log((content(pa,hb,p50,n)-content(y,hb,p50,n))/(content(pa,hb,p50,n)-content(pv,hb,p50,n)))
    if law=='pressure':return -q*math.log((pa-y)/(pa-pv))
    if law=='equilibrium_admixture':return (content(pa,hb,p50,n)-content(y,hb,p50,n))/(content(pa,hb,p50,n)-content(pv,hb,p50,n))
    return 10*q*quad(lambda p:cp(p,hb,p50,n)/(pa-p),pv,y,epsabs=1e-11,epsrel=1e-11)[0]

def fit_many(law,rows):
    if law=='equilibrium_admixture':
        res=minimize_scalar(lambda s:sum((predict(law,s,r)-r['Pa_mmHg'])**2 for r in rows),bounds=(0,.8),method='bounded',options={'xatol':1e-10})
        return float(res.x)
    res=minimize_scalar(lambda x:sum((predict(law,math.exp(x),r)-r['Pa_mmHg'])**2 for r in rows),bounds=(-5,6),method='bounded',options={'xatol':1e-9})
    return math.exp(res.x)

def parameter_unit(law,flow_unit):
    if law=='equilibrium_admixture':return 'dimensionless'
    if law=='constant_D':return 'mL_O2/min/mmHg/m2' if flow_unit.endswith('/m2') else 'mL_O2/min/mmHg'
    return flow_unit

def record(r,law,param,protocol,pred):
    return dict(study=r['study'],condition=r['condition'],population=r['population'],protocol=protocol,law=law,
                parameter=param,parameter_basis=r['flow_unit'],parameter_unit=parameter_unit(law,r['flow_unit']),observed_Pa=r['Pa_mmHg'],predicted_Pa=pred,
                error_mmHg=pred-r['Pa_mmHg'],absolute_error_mmHg=abs(pred-r['Pa_mmHg']),
                PA=r['PA_mmHg'],Pv=r['Pv_mmHg'],flow=r['flow'],Hb=r['Hb_g_dl'])

def within(rows,offset=0):
    results=[];fits=[]
    for study in ['Carta2022','Lichtblau2021','Blanco2010','Blanco2011']:
        rr=[r for r in rows if r['study']==study];base=next(r for r in rr if r['baseline'])
        for law in LAWS:
            par=fit_one(law,base,base['Pa_mmHg']+offset)
            fits.append(dict(study=study,law=law,parameter=par,flow_unit=base['flow_unit'],parameter_unit=parameter_unit(law,base['flow_unit']),baseline_residual=predict(law,par,base)-base['Pa_mmHg']-offset))
            for r in rr:
                if not r['baseline']:results.append(record(r,law,par,'baseline_frozen',predict(law,par,r)))
    return results,fits

def metrics(predictions):
    rows=[]
    for protocol,study in sorted(set((r['protocol'],r['study']) for r in predictions)):
        for law in LAWS:
            rr=[r for r in predictions if (r['protocol'],r['study'],r['law'])==(protocol,study,law)]
            e=np.array([r['error_mmHg'] for r in rr])
            rows.append(dict(protocol=protocol,study=study,law=law,test_rows=len(rr),MAE_mmHg=float(np.mean(abs(e))),RMSE_mmHg=float(np.sqrt(np.mean(e*e))),bias_mmHg=float(e.mean())))
    return rows

def main():
    rows=read_data();preds,fits=within(rows)
    dant=[r for r in rows if r['study']=='Dantzker1979'];base=next(r for r in rows if r['study']=='Carta2022' and r['baseline'])
    for i,test in enumerate(dant):
        for law in LAWS:
            p=fit_many(law,dant[:i]+dant[i+1:]);preds.append(record(test,law,p,'leave_one_patient_out',predict(law,p,test)))
            p=fit_one(law,base);preds.append(record(test,law,p,'Carta_frozen_external',predict(law,p,test)))
    write('baseline_fits.csv',fits);write('predictions.csv',preds);mm=metrics(preds);write('metrics.csv',mm)
    sensitivity=[];wins={}
    import itertools
    for p50,n,rq,hb in itertools.product([24.,26.8,30.],[2.5,2.7,3.],[.75,.8,.9],[12.,15.,18.]):
        pp,_=within(read_data(p50,n,rq,hb))
        for m in metrics(pp):sensitivity.append(dict(P50=p50,hill=n,R=rq,missing_Hb=hb,**m))
    write('assumption_sensitivity.csv',sensitivity)
    # Distinct scenario counts; missing Hb has no effect in the two Zurich studies.
    winrows=[]
    for study in ['Carta2022','Lichtblau2021','Blanco2010','Blanco2011']:
        rr=[r for r in sensitivity if r['study']==study and (study.startswith('Blanco') or r['missing_Hb']==15)]
        keys=sorted(set((r['P50'],r['hill'],r['R'],r['missing_Hb']) for r in rr))
        win={law:0 for law in LAWS}
        for k in keys:
            group=[r for r in rr if (r['P50'],r['hill'],r['R'],r['missing_Hb'])==k]
            win[min(group,key=lambda r:r['RMSE_mmHg'])['law']]+=1
        for law,count in win.items():winrows.append(dict(study=study,law=law,winning_scenarios=count,total_scenarios=len(keys),interpretation='Deterministic assumption counts; not probabilities'))
    write('sensitivity_winners.csv',winrows)
    external_sensitivity=[]
    for p50,n,rq in itertools.product([24.,26.8,30.],[2.5,2.7,3.],[.75,.8,.9]):
        rr=read_data(p50,n,rq);bb=next(r for r in rr if r['study']=='Carta2022' and r['baseline'])
        for law in LAWS:
            par=fit_one(law,bb);ee=[predict(law,par,r)-r['Pa_mmHg'] for r in rr if r['study']=='Dantzker1979']
            external_sensitivity.append(dict(P50=p50,hill=n,R=rq,law=law,RMSE_mmHg=float(np.sqrt(np.mean(np.square(ee))))))
    write('external_sensitivity.csv',external_sensitivity)
    rounding=[]
    for offset in [-.5,0,.5]:
        pp,_=within(rows,offset)
        rounding += [dict(baseline_offset_mmHg=offset,**m) for m in metrics(pp)]
    write('baseline_rounding.csv',rounding)
    # Mismatch of the shared blood curve is not hidden in a fitted exchange parameter.
    blood=[]
    for r in rows:
        if r['Sa_percent']!='':
            blood.append(dict(study=r['study'],condition=r['condition'],reported_Sa=float(r['Sa_percent']),
                curve_Sa=100*sat(r['Pa_mmHg']),reported_Sv=float(r['Sv_percent']),curve_Sv=100*sat(r['Pv_mmHg'])))
    write('blood_curve_audit.csv',blood)
    checks=dict(baseline_fit_max_abs=max(abs(r['baseline_residual']) for r in fits),
                predictions_within_inlet_equilibrium=all(r['Pv']-1e-7<=r['predicted_Pa']<=r['PA']+1e-7 for r in preds))
    # Independent original content-state ODE at representative actual conditions.
    spec=importlib.util.spec_from_file_location('legacy_ox',ROOT.parent/'aggregate_analysis/code/microcirculation.py')
    ox=importlib.util.module_from_spec(spec);import sys;sys.modules['legacy_ox']=ox;spec.loader.exec_module(ox)
    errs=[]
    for r in [rows[0],rows[1],rows[2],rows[7],rows[13]]:
        d=fit_one('constant_D',r);p=predict('constant_D',d,r)
        c=ox.end_capillary_content(r['flow'],d,r['PA_mmHg'],content(r['Pv_mmHg'],r['Hb_g_dl']),hb=r['Hb_g_dl'],rtol=1e-10)
        errs.append(abs(content(p,r['Hb_g_dl'])-c))
    checks['independent_content_ODE_max_ml_dl']=max(errs)
    assert checks['baseline_fit_max_abs']<1e-6 and checks['predictions_within_inlet_equilibrium'] and max(errs)<1e-6
    summary=dict(observation_rows=len(rows),studies=5,supporting_studies=2,supporting_conditions=6,
                 distinct_within_cohort_test_conditions=11,individual_LOPO_test_rows=7,
                 independent_external_Carta_to_Dantzker_rows=7,checks=checks,
                 python=platform.python_version(),numpy=np.__version__,
                 warning='No independent pooled participant N; conditional aggregate predictions, possible Zurich overlap; no globally validated best physiological law.')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2));print(json.dumps(mm,indent=2))

if __name__=='__main__':main()
