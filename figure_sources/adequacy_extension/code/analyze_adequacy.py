"""Output-specific deterministic adequacy, not clinical model validation."""
from pathlib import Path
import csv, json, math, sys
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'heterogeneity_extension/code'))
import analyze_heterogeneity as old
HERE=Path(__file__).resolve().parents[1]

def residual_bounds(mean,r,beta_lo=0.,beta_hi=1.,lo=0.,hi=1.):
    if not all(math.isfinite(x) for x in (mean,r,beta_lo,beta_hi,lo,hi)) or not 0<=lo<=mean<=hi<=1 or not 0<mean<1 or r<1 or not 0<=beta_lo<=beta_hi<=1:
        raise ValueError('invalid residual class')
    al=r**(-beta_lo); au=r**(-beta_hi)
    low=mean**al if lo==hi else ((hi-mean)*lo**al+(mean-lo)*hi**al)/(hi-lo)
    return low,mean**au

def content_interval(sl,su,q,ce,m):
    if not all(math.isfinite(x) for x in (sl,su,q,ce,m)) or not 0<=sl<=su<1 or q<=0 or ce<=0 or m<0:
        raise ValueError('invalid output domain')
    cap=1-m/(10*q*ce)
    if sl>cap+1e-13:
        return dict(status='infeasible',low=math.nan,high=math.nan,width=math.nan,midpoint=math.nan,radius=math.nan)
    su=min(su,cap)
    f=lambda s:ce-m/(10*q)*s/(1-s)
    low,high=f(su),f(sl)
    return dict(status='feasible',low=low,high=high,width=high-low,midpoint=(low+high)/2,radius=(high-low)/2)

def certificate(low,high,prediction,tolerance):
    if low>high or tolerance<0:raise ValueError('invalid certificate')
    risk=max(abs(prediction-low),abs(prediction-high))
    return risk,risk<=tolerance+1e-12

def two_state_distributions(s0,s1):
    """Moments E[y]=s1 and E[y^2]=s0, y in [0,1]."""
    if not 0<s0<1 or not s0-1e-12<=s1<=math.sqrt(s0)+1e-12:
        raise ValueError('incompatible residual moments')
    s1=min(max(s1,s0),math.sqrt(s0))
    b=s0/s1
    w=s1*s1/s0
    c=(s1-s0)/(1-s1)
    v=(1-s1)/(1-c)
    return (np.array([1-w,w]),np.array([0.,b])),(np.array([v,1-v]),np.array([c,1.]))

def two_state_bounds(s0,s1,r,beta=1.):
    if r<1 or not 0<beta<=1:raise ValueError('r>=1 and known beta in (0,1] required')
    d0,d1=two_state_distributions(s0,s1)
    p=2*r**(-beta)
    values=[float(w @ y**p) for w,y in (d0,d1)]
    return min(values),max(values)

def noisy_two_state_bounds(s0,s1lo,s1hi,r,beta=1.):
    if r<1 or not 0<beta<=1 or s1lo>s1hi:raise ValueError('invalid noisy domain')
    sl=max(s0,s1lo);su=min(math.sqrt(s0),s1hi)
    if sl>su:raise ValueError('empty noisy moment class')
    a=r**(-beta);p=2*a
    if p>=1:return max(s0,sl**p),min(s0**a,su)
    return max(s0,sl),min(s0**a,su**p)

def general_slope(q,ce,m,s,eta,ce_log_slope,s_log_slope):
    return ce_log_slope+m/(10*q)*((1-eta)*s/(1-s)-s_log_slope/(1-s)**2)

def guaranteed_slope(q,ce,m,s,eta,ce_log_slope,beta_max):
    if not 0<s<1 or not 0<=beta_max<=1:raise ValueError('invalid direction class')
    return general_slope(q,ce,m,s,eta,ce_log_slope,beta_max*(-s*math.log(s)))

def write(name,rows):
    dest=HERE/'results';dest.mkdir(exist_ok=True)
    with (dest/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    rows=[]
    for s in [.05,.10,.20,.40]:
        for beta in [0.,.25,.5,1.]:
            for r in np.linspace(1,2,101):
                sl,su=residual_bounds(s,float(r),beta,beta)
                b=content_interval(sl,su,old.Q0*r,old.CE,old.M)
                rows.append(dict(s0=s,beta=beta,r=r,s_lower=sl,s_upper=su,**b))
    write('adequacy_grid.csv',rows)
    examples=[]
    for label,s,bl,bh,lo,hi,ce,m in [
        ('low_residual',.05,1,1,0,1,old.CE,250),
        ('moderate_residual',.2,1,1,0,1,old.CE,250),
        ('proportional_capacity',.2,0,0,0,1,old.CE,250),
        ('weak_flow_dependence',.2,.25,.25,0,1,old.CE,250),
        ('unknown_beta',.2,0,1,0,1,old.CE,250),
        ('narrow_support',.2,1,1,.15,.25,old.CE,250),
        ('changed_gas_and_demand',.2,1,1,0,1,17.,375.)]:
        sl,su=residual_bounds(s,2,bl,bh,lo,hi)
        b=content_interval(sl,su,10,ce,m)
        examples.append(dict(case=label,s0=s,beta_low=bl,beta_high=bh,z_low=lo,z_high=hi,r=2,Ce=ce,M=m,**b,endpoint_ok_0p2=b['width']<=.2,midpoint_ok_0p2=b['radius']<=.2))
    write('adequacy_examples.csv',examples)
    # Closed-form width roots. Full [0,1] support, known beta=1, r=2.
    roots=[]
    for tol in [.1,.2]:
        for factor,name in [(1,'endpoint'),(2,'minimax')]:
            f=lambda s:content_interval(s,math.sqrt(s),10,old.CE,250)['width']-factor*tol
            sol=brentq(f,1e-9,.5)
            roots.append(dict(tolerance=tol,predictor=name,maximum_s0=sol,r=2,beta=1,Q=10,M=250))
    write('tolerance_thresholds.csv',roots)
    demo=[]
    for name,(w,z) in list(old.models().items())[:2]:
        s0=float(w@z);s1=float(w@np.sqrt(z))
        for r in np.linspace(1,3,81):
            actual=old.state(w,z,float(r),1)['Ca']
            for method in ['baseline_only','two_exact','two_noisy']:
                if method=='baseline_only':sl,su=residual_bounds(s0,float(r),1,1)
                elif method=='two_exact':sl,su=two_state_bounds(s0,s1,float(r))
                else:sl,su=noisy_two_state_bounds(s0,s1-.005,s1+.005,float(r))
                b=content_interval(sl,su,old.Q0*r,old.CE,old.M)
                demo.append(dict(model=name,method=method,s0=s0,s_second=s1,r=r,actual=actual,**b))
    write('second_state_bounds.csv',demo)
    slopes=[]
    for s in [.05,.2,.4]:
        for eta in [0.,.5,1.,1.5]:
            for ge in [-.5,0.,.5]:
                for beta in [0.,.25,.5,1.]:
                    bound=guaranteed_slope(5,old.CE,250,s,eta,ge,beta)
                    slopes.append(dict(s=s,eta=eta,Ce_log_slope=ge,beta_max=beta,lower_dCa_dlogQ=bound,all_models_increase=bound>0))
    write('direction_conditions.csv',slopes)
    # Summarize existing observations; this adds no cohort and no refitting.
    with (ROOT/'review_extension/results/matched_metrics.csv').open(encoding='utf-8-sig') as f:metrics=list(csv.DictReader(f))
    metricrows=[r for r in metrics if r['scope'] in ['air_hypoxia','modern','thermodilution'] and r['model'] in ['content_beta0','content_beta1','saturation_beta0','saturation_beta1']]
    write('empirical_model_choices.csv',metricrows)
    e=next(x for x in examples if x['case']=='moderate_residual')
    inv=lambda c:brentq(lambda p:old.prior.v.content(p,14.)-c,0.,1000.)
    pl,pu=inv(e['low']),inv(e['high'])
    pc=inv(e['midpoint']);pp=(pl+pu)/2
    write('output_scale_example.csv',[dict(Ca_low=e['low'],Ca_high=e['high'],Pa_low=pl,Pa_high=pu,content_midpoint_Pa=pc,pressure_midpoint_Pa=pp,pressure_risk_at_content_midpoint=max(pc-pl,pu-pc),pressure_minimax_risk=(pu-pl)/2)])
    summary=dict(examples=examples,tolerance_thresholds=roots,second_state_at_1p5=[d for d in demo if abs(d['r']-1.5)<1e-12],second_state_at_3=[d for d in demo if d['r']==3],interpretation='Attainable model-class intervals and absolute-error minimax guarantees; no patient validation, no clinical tolerance estimated')
    (HERE/'results/summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
