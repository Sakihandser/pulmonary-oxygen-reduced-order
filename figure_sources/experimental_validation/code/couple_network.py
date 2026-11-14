"""Hypothetical PAH steady states using the unchanged anatomical flow solver.

The cardiac pressure-flow closure is inherited, not clinically validated.
Effective admixture is not a measured anatomical shunt or regional V/Q map.
"""
from pathlib import Path
import importlib.util,sys,json,itertools
import numpy as np
from scipy.optimize import brentq
from validate import ROOT,OUT,read_data,fit_one,content,sat,invcontent,predict,write,LAWS,CONFIG,P50,HILL
from primary_model import solve as solve_primary

spec=importlib.util.spec_from_file_location('existing_network',ROOT.parent/'hydraulic_network/code/build_lung_model.py')
lung=importlib.util.module_from_spec(spec);spec.loader.exec_module(lung)
case=lung.build_cases();lobes=list(lung.LOBE_W);shares=np.array(list(lung.LOBE_W.values()))
base=next(r for r in read_data() if r['study']==CONFIG['calibration_study'] and r['baseline'])
pars={law:fit_one(law,base) for law in LAWS}

def hydraulic(pvr,k=1):
    def solve(scale):return lung.solve(case,lobe_arterial_multiplier={l:scale if l=='LUL' else scale*k for l in lobes})
    a=brentq(lambda s:solve(s)['PVR_WU']-pvr,.001,1000,xtol=1e-10)
    return solve(a)

def closed(h,law,hb=14,pa=None,vo2=250,admixture=None,capacity_scale=1.):
    pa=base['PA_mmHg'] if pa is None else pa
    q=h['CO_L_min'];f=np.array([h['flow_fraction'][l] for l in lobes]);p50=P50;n=HILL
    if law==CONFIG['primary_law']:
        state=solve_primary(q,hb,pa,vo2,pars[law] if admixture is None else admixture)
        return {k:v for k,v in state.items() if k not in ['Q_min_L_min','equilibrium_content_ml_dl']}
    ce=content(pa,hb,p50,n)
    def lungs(cv):
        pv=invcontent(cv,hb,p50,n)
        rr=dict(base,Hb_g_dl=hb,PA_mmHg=pa,Pv_mmHg=pv)
        cc=[]
        for fi,wi in zip(f,shares):
            p=predict(law,pars[law]*capacity_scale*wi,dict(rr,flow=q*fi))
            cc.append(content(p,hb,p50,n))
        return float(f@cc)
    demand_gap=vo2/(10*q)
    upper_cv=ce-demand_gap
    if upper_cv<=0 or lungs(0)<demand_gap:
        return dict(feasible=False,Pa_mmHg='',Sa_percent='',Ca_ml_dl='',Cv_ml_dl='',Sv_percent='',DO2_ml_min='',mass_balance_error='')
    cv=brentq(lambda cv:lungs(cv)-cv-demand_gap,0,upper_cv,xtol=1e-10)
    ca=lungs(cv);p=invcontent(ca,hb,p50,n);pv=invcontent(cv,hb,p50,n)
    return dict(feasible=True,Pa_mmHg=p,Sa_percent=100*sat(p),Ca_ml_dl=ca,Cv_ml_dl=cv,
                Sv_percent=100*sat(pv),DO2_ml_min=10*q*ca,mass_balance_error=10*q*(ca-cv)-vo2)

def main():
    results=[]
    for pvr,k,hb,vo2,pa,law in itertools.product([1.,3.,5.,8.,12.,16.],[1,25],[10.,14.],[250.,500.],
                                                 [base['PA_mmHg'],150.],LAWS):
        h=hydraulic(pvr,k);o=closed(h,law,hb,pa,vo2)
        results.append(dict(PVR_WU=pvr,contrast=k,law=law,Hb=hb,VO2=vo2,PA=pa,
                            Q=h['CO_L_min'],mPAP=h['mPAP_mmHg'],LUL_flow_fraction=h['flow_fraction']['LUL'],**o))
    write('network_scenarios.csv',results)
    sh=[]
    for pvr,s in itertools.product([1.,3.,5.,8.,12.,16.],[.03,pars['equilibrium_admixture'],.3]):
        h=hydraulic(pvr)
        sh.append(dict(PVR_WU=pvr,effective_admixture=s,Q=h['CO_L_min'],mPAP=h['mPAP_mmHg'],**closed(h,'equilibrium_admixture',admixture=s)))
    write('network_admixture_sensitivity.csv',sh)
    # Real independently measured low-V/Q fractions: rough mechanistic control,
    # treating all blood with V/Q<0.1 as fully unoxygenated is intentionally coarse.
    controls=[]
    for r in read_data():
        if r['study']!='Dantzker1979':continue
        for key in ['measured_shunt_fraction','low_vq_fraction']:
            p=predict('equilibrium_admixture',float(r[key]),r)
            controls.append(dict(condition=r['condition'],input_fraction=key,fraction=float(r[key]),
                                observed_Pa=r['Pa_mmHg'],predicted_Pa=p,error=p-r['Pa_mmHg']))
    write('measured_mixing_controls.csv',controls)
    # Prespecified recruitment sensitivities; b is never fitted to test outputs.
    recruitment=[]
    for study in ['Carta2022','Lichtblau2021','Blanco2010','Blanco2011']:
        rr=[r for r in read_data() if r['study']==study];bb=next(r for r in rr if r['baseline'])
        d0=fit_one('constant_D',bb)
        for b in [0,.5,1.]:
            for r in rr:
                if r['baseline']:continue
                d=d0*(r['flow']/bb['flow'])**b;p=predict('constant_D',d,r)
                recruitment.append(dict(study=study,condition=r['condition'],exponent=b,D_effective=d,
                                        observed_Pa=r['Pa_mmHg'],predicted_Pa=p,error=p-r['Pa_mmHg']))
    write('recruitment_sensitivity.csv',recruitment)
    checks=dict(mass_balance_max=max(abs(r['mass_balance_error']) for r in results if r['feasible']),
                hydraulic_max= max(abs(hydraulic(p)['flow_conservation_error_L_min']) for p in [1,3,5,8,12,16]),
                scenarios=len(results),infeasible=sum(not r['feasible'] for r in results),
                admixture=pars['equilibrium_admixture'],PA_reference=base['PA_mmHg'],
                interpretation='Conditional mechanistic scenarios, not outcome-validated PAH predictions. Regional V/Q and cardiac adaptation are unmeasured.')
    assert checks['mass_balance_max']<1e-6 and checks['hydraulic_max']<1e-7
    (OUT/'network_summary.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
    print(json.dumps(checks,indent=2))
    for r in results:
        if r['law']=='equilibrium_admixture' and r['contrast']==1 and r['Hb']==14 and r['VO2']==250 and r['PA']==base['PA_mmHg']:print(r)

if __name__=='__main__':main()
