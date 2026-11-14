"""Outcome-based compatibility diagnostics; no clinical hypothesis tests."""
from pathlib import Path
import sys,csv,json,math,itertools
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[2]; HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'selection_extension/code'))
import analyze as a
v=a.v

def envelope(s0,r):
    if not 0<=s0<=1 or not math.isfinite(r) or r<=0:raise ValueError('invalid state')
    x=s0**(1/r)
    return min(s0,x),max(s0,x)

def interval_compatible(blo,bhi,tlo,thi,rlo,rhi):
    if not 0<=blo<=bhi<=1 or not 0<=tlo<=thi<=1 or not 0<rlo<=rhi:raise ValueError('invalid interval')
    # Both bounds increase in s0. The entire continuous union over r is an interval.
    low=min(envelope(blo,rlo)[0],envelope(blo,rhi)[0])
    high=max(envelope(bhi,rlo)[1],envelope(bhi,rhi)[1])
    return max(low,tlo)<=min(high,thi)+1e-12

def residual_radius(s0,s1,r):
    f=lambda e:interval_compatible(max(0,s0-e),min(1,s0+e),max(0,s1-e),min(1,s1+e),r,r)
    if f(0):return 0.
    lo=0.;hi=1.
    for _ in range(60):
        mid=(lo+hi)/2
        if f(mid):hi=mid
        else:lo=mid
    return hi

def contents(row,convention='hill',shift=0.,pa_offset=0.):
    p50=row['p50']+shift;hb=row['Hb_g_dl'];n=row['hill']
    ce=v.content(row['PA_mmHg'],hb,p50,n)
    if convention=='reported_saturation':
        if not row.get('Sa_percent') or not row.get('Sv_percent'):raise ValueError('reported saturation missing')
        ca=1.34*hb*float(row['Sa_percent'])/100+.0031*row['Pa_mmHg']
        cv=1.34*hb*float(row['Sv_percent'])/100+.0031*row['Pv_mmHg']
    else:
        ca=v.content(row['Pa_mmHg']+pa_offset,hb,p50,n)
        cv=v.content(row['Pv_mmHg'],hb,p50,n)
    if ce<=cv:raise ValueError('Ce must exceed Cv')
    s=(ce-ca)/(ce-cv)
    if not 0<=s<=1:raise ValueError('observed residual outside physical domain')
    return ce,ca,cv,s

def pa_interval(row,e):
    ce,ca,cv,s=contents(row)
    high_pa=min(row['PA_mmHg'],row['Pa_mmHg']+e)
    low_pa=max(row['Pv_mmHg'],row['Pa_mmHg']-e)
    fun=lambda p:(ce-v.content(p,row['Hb_g_dl'],row['p50'],row['hill']))/(ce-cv)
    return max(0.,fun(high_pa)),min(1.,fun(low_pa))

def pa_radius(b,t,r):
    f=lambda e:interval_compatible(*pa_interval(b,e),*pa_interval(t,e),r,r)
    if f(0):return 0.
    lo=0.;hi=1000.
    for _ in range(60):
        mid=(lo+hi)/2
        if f(mid):hi=mid
        else:lo=mid
    return hi

def pairs(rows,secondary=False):
    for unit,rr in a.unit_groups(rows).items():
        b=next(x for x in rr if x['baseline'])
        if not b['is_primary_PH'] and not (secondary and b['study']=='Lichtblau2025'):continue
        for t in rr:
            if not t['baseline']:yield unit,b,t

def record(unit,b,t,convention='hill',shift=0.):
    ce0,ca0,cv0,s0=contents(b,convention)
    ce1,ca1,cv1,s1=contents(t,convention,shift)
    r=t['flow']/b['flow'];lo,hi=envelope(s0,r)
    method='thermodilution' if t['cohort_group']=='Zurich_early' else ('direct_Fick' if t['cohort_group']=='Zurich2021_2023' else 'not_independently_established')
    ok=lo-1e-12<=s1<=hi+1e-12
    beta=math.log(math.log(s1)/math.log(s0))/(-math.log(r)) if abs(r-1)>1e-12 and 0<s0<1 and 0<s1<1 else math.nan
    return dict(unit=unit,study=t['study'],condition=t['condition'],FiO2=t['FiO2'],primary=int(t['is_primary_PH']),flow_method=method,granularity=t['granularity'],r=r,s0=s0,s1=s1,s_lower=lo,s_upper=hi,compatible=int(ok),residual_halfwidth=residual_radius(s0,s1,r),homogeneous_beta=beta,convention=convention,changed_P50_shift=shift)

def write(name,rows):
    dest=HERE/'results';dest.mkdir(exist_ok=True)
    with (dest/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    rows=a.load();nominal=[];contrast=[];alternative=[]
    for unit,b,t in pairs(rows,True):
        rec=record(unit,b,t);r=rec['r']
        rec['Pa_halfwidth_mmHg']=pa_radius(b,t,r)
        rec['compatible_Qratio_plusminus_10percent']=int(interval_compatible(rec['s0'],rec['s0'],rec['s1'],rec['s1'],.9*r,1.1*r))
        nominal.append(rec)
        ce,_,cv,_=contents(t);s0=rec['s0']
        # The baseline-frozen conditional predictions consume current Cv, not current Ca.
        p0=v.invcontent(ce-(ce-cv)*s0,t['Hb_g_dl'],t['p50'],t['hill'])
        p1=v.invcontent(ce-(ce-cv)*s0**(1/r),t['Hb_g_dl'],t['p50'],t['hill'])
        contrast.append(dict(unit=unit,condition=t['condition'],primary=int(t['is_primary_PH']),FiO2=t['FiO2'],flow_method=rec['flow_method'],r=r,Pa_beta0=p0,Pa_beta1=p1,separation_mmHg=abs(p1-p0),equal_error_budget_mmHg=abs(p1-p0)/2))
        if b.get('Sa_percent') and b.get('Sv_percent') and t.get('Sa_percent') and t.get('Sv_percent'):
            alternative.append(record(unit,b,t,'reported_saturation'))
    write('pair_compatibility.csv',nominal);write('prediction_separation.csv',contrast);write('reported_saturation_audit.csv',alternative)
    sensitivity=[]
    for p50,n,rq,hb in itertools.product([24.,26.8,30.],[2.5,2.7,3.],[.75,.8,.9,1.],[12.,15.,18.]):
        for unit,b,t in pairs(a.load(p50,n,rq,hb)):
            z=record(unit,b,t)
            sensitivity.append(dict(P50=p50,Hill=n,R=rq,missing_Hb=hb,unit=unit,condition=t['condition'],FiO2=t['FiO2'],compatible=z['compatible'],residual_halfwidth=z['residual_halfwidth']))
    write('input_sensitivity.csv',sensitivity)
    shifts=[]
    for shift in [-2.,-1.,0.,1.,2.]:
        for unit,b,t in pairs(rows):shifts.append(record(unit,b,t,shift=shift))
    write('state_curve_sensitivity.csv',shifts)
    summary_rows=[]
    for rec in nominal:
        if not rec['primary']:continue
        subset=[x for x in sensitivity if x['unit']==rec['unit'] and x['condition']==rec['condition']]
        summary_rows.append(dict(unit=rec['unit'],condition=rec['condition'],FiO2=rec['FiO2'],compatible_scenarios=sum(x['compatible'] for x in subset),scenarios=len(subset),min_residual_halfwidth=min(x['residual_halfwidth'] for x in subset),max_residual_halfwidth=max(x['residual_halfwidth'] for x in subset)))
    write('sensitivity_summary.csv',summary_rows)
    # With exact s0, distinguish homogeneous beta 0 from 1 at an arbitrary obtainable flow.
    design=[]
    for s0,r,e in itertools.product([.05,.15,.3],[1.,1.05,1.1,1.25,1.5,2.],[0.,.005,.01,.02]):
        sep=s0**(1/r)-s0
        design.append(dict(s0=s0,r=r,residual_halfwidth=e,residual_separation=sep,disjoint=int(sep>2*e+1e-12),interpretation='homogeneous_candidates_exact_baseline_only'))
    write('measurement_separation.csv',design)
    pri=[x for x in nominal if x['primary'] and x['FiO2']<.3]
    therm=[x for x in contrast if x['primary'] and x['FiO2']<.3 and x['flow_method']=='thermodilution']
    alt=[x for x in alternative if x['primary'] and x['FiO2']<.3]
    info=dict(nonhyperoxic_pairs=len(pri),nominal_incompatible=sum(1-x['compatible'] for x in pri),compatible_pairs=[x['unit']+':'+x['condition'] for x in pri if x['compatible']],thermodilution_separations=[x['separation_mmHg'] for x in therm],reported_saturation_nonhyperoxic_pairs=len(alt),reported_saturation_incompatible=sum(1-x['compatible'] for x in alt),sensitivity_summary=summary_rows,secondary=[x for x in nominal if not x['primary']],interpretation='Exploratory aggregate compatibility audit; not statistical rejection, independent validation or new physiological mechanism.')
    (HERE/'results/summary.json').write_text(json.dumps(info,indent=2,allow_nan=True),encoding='utf-8')
    print(json.dumps({k:v for k,v in info.items() if k not in ['sensitivity_summary','secondary']},indent=2))

if __name__=='__main__':main()
