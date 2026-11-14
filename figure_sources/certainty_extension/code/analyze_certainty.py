from pathlib import Path
import sys,csv,json,math,itertools
import numpy as np
from scipy.special import erfc,erfcinv
from scipy.optimize import brentq
from consistent import bounds,margin,content,classify,minimum_mean,GUARD
HERE=Path(__file__).resolve().parents[1];ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'applicability_extension/code'))
import analyze_applicability as ap

def write(name,rows):
 with (HERE/'results'/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def levy_discrete(a,U,n):
 F=float(erfc(a/(2*math.sqrt(U))))
 p=(np.arange(n)+.5)*F/n
 x=a*a/(4*erfcinv(p)**2)
 return np.r_[x,U],np.r_[np.full(n,F/n),1-F],F

def levy_mean(a,U):
 z=a/(2*math.sqrt(U));F=erfc(z)
 return a*math.sqrt(U)/math.sqrt(math.pi)*math.exp(-z*z)-a*a/2*F+U*(1-F)

def levy_bound(a,U,n,rmax=2.):
 F=float(erfc(a/(2*math.sqrt(U))))
 # Censoring and monotone quantile-midpoint quadrature, uniformly over r<=rmax.
 return (1-F)*math.exp(-U/rmax)+F/(2*n)

def structures():
 w=(.8-.2)/(.8-.001);xA=-np.log([.001,.8]);mu=float(np.dot([w,1-w],xA))
 def fun(x):
  ww=(.2-math.exp(-20))/(math.exp(-x)-math.exp(-20))
  return ww*x+(1-ww)*20-mu
 x=brentq(fun,.01,-math.log(.2));ww=(.2-math.exp(-20))/(math.exp(-x)-math.exp(-20))
 return {'A':(xA,np.array([w,1-w])),'B':(np.array([x,20.]),np.array([ww,1-ww]))},mu

def main():
 (HERE/'results').mkdir(exist_ok=True)
 Ce=ap.v.content(100.,14.);M=250.;Q0=5.;target=content(.15,Q0,Ce,M)
 models,mu=structures();rows=[];wit=[];grid=[];curves=[]
 for label,(x,w) in models.items():
  for r in np.linspace(1,2,101):curves.append(dict(model=label,r=r,s=float(np.exp(-x/r)@w),Ca=content(float(np.exp(-x/r)@w),Q0*r,Ce,M)))
  cases=[('baseline',[1.],.001,None),('baseline_capacity',[1.],.001,(mu-1e-4,mu+1e-4)),
         ('one_extra',[1.,1.25],.001,None),('two_extra',[1.,1.25,1.5],.001,None),
         ('two_extra_noise005',[1.,1.25,1.5],.005,None),('two_extra_noise010',[1.,1.25,1.5],.01,None)]
  for name,rr,e,moment in cases:
   vals=np.exp(-x[None,:]/np.array(rr)[:,None])@w
   sizes=[128,512,2048]+([8192,32768] if label=='A' and name=='two_extra' else [])
   for n in sizes:
    out=bounds(rr,np.maximum(0,vals-e),np.minimum(1,vals+e),2.,n=n,moment=moment)
    if out['status']!='bounded':raise RuntimeError((label,name,n,out))
    decision=classify(out,10.,Ce,M,target)
    grid.append(dict(model=label,information=name,n=n,s_lower=out['outer_lower'],s_upper=out['outer_upper'],inner_lower=out['inner_lower'],inner_upper=out['inner_upper'],gap=out['enclosure_gap'],decision=decision,dual_padding=out['dual_padding']))
   low=content(out['outer_upper'],10.,Ce,M);high=content(out['outer_lower'],10.,Ce,M)
   ulow=content(out['outer_upper'],10.,Ce-.1,275.);uhigh=content(out['outer_lower'],10.,Ce+.1,225.)
   uncertain='all_recover' if ulow>target+GUARD else ('none_recover' if uhigh<target-GUARD else 'not_certified')
   rows.append(dict(model=label,information=name,residual_halfwidth=e,observed_ratios=';'.join(map(str,rr)),cells=n,s_lower=out['outer_lower'],s_upper=out['outer_upper'],Ca_lower=low,Ca_upper=high,target=target,decision=decision,CeM_lower=ulow,CeM_upper=uhigh,CeM_decision=uncertain,enclosure_gap=out['enclosure_gap']))
   for witness in out['witnesses']:
    for xx,ww in zip(witness['x'],witness['weights']):wit.append(dict(model=label,information=name,direction=witness['direction'],exposure=xx,weight=ww,query_s=witness['value'],constraint_violation=witness['maximum_constraint_violation']))
 write('decision_bounds.csv',rows);write('finite_witnesses.csv',wit);write('grid_convergence.csv',grid);write('reference_curves.csv',curves)
 # Exact transformation counterexample and finite-capacity approximations.
 a=-math.log(.2);alias=[];paths=[];aliascurves=[]
 for U,n in itertools.product([8.,16.,32.],[256,1024,4096]):
  x,w,F=levy_discrete(a,U,n);rs=np.linspace(.5,2,301)
  truth=np.exp(-a/np.sqrt(rs));pred=np.exp(-x[None,:]/rs[:,None])@w
  bound=levy_bound(a,U,n);err=float(max(abs(truth-pred)))
  assert err<=bound+GUARD
  # Under current measured Cv and common Ce, conditional content error <=(Ce-Cv)*residual error.
  alias.append(dict(exposure_cap=U,paths=len(x),residual_bound=bound,max_grid_error=err,mean_exposure=float(x@w),analytic_censored_mean=levy_mean(a,U),K_total_Q0_5=5*float(x@w),homogeneous_K0=5*a,content_bound_for_Ce_minus_Cv_5=5*bound))
  if n==4096:
   for xx,ww in zip(x,w):paths.append(dict(exposure_cap=U,exposure=xx,weight=ww))
   for r,ss,pp in zip(rs,truth,pred):aliascurves.append(dict(exposure_cap=U,r=r,fractional_s=ss,fixed_mixture_s=pp,error=pp-ss))
 write('mechanism_alias.csv',alias);write('alias_paths.csv',paths);write('alias_curves.csv',aliascurves)
 capacity=[]
 for e in [.0001,.001,.005,.01]:
  rr=np.array([.5,.75,1.,1.25,1.5,2.]);vv=np.exp(-a/np.sqrt(rr))
  for n in [2048,8192]:
   m=minimum_mean(rr,vv-e,vv+e,n=n)
   capacity.append(dict(residual_halfwidth=e,cells=n,K_lower=5*m['mean_lower'],K_upper=5*m['mean_upper'],gap=5*(m['mean_upper']-m['mean_lower']),dual_padding=m['dual_padding'],conditional_capacity_10_excluded=int(5*m['mean_lower']>10+GUARD)))
 write('capacity_exclusion.csv',capacity)
 # All states for one unit share one model; errors here are hypothetical residual halfwidths.
 empirical=[]
 for unit,rr in ap.a.unit_groups(ap.a.load()).items():
  b=next(z for z in rr if z['baseline'])
  if not b['is_primary_PH']:continue
  rr=[z for z in rr if z['FiO2']<.3]
  if len(rr)<2:continue
  ratios=[z['flow']/b['flow'] for z in rr];values=[ap.contents(z)[3] for z in rr]
  pair=max(ap.residual_radius(ap.contents(b)[3],ap.contents(z)[3],z['flow']/b['flow']) for z in rr)
  for n in [512,2048]:
   m=margin(ratios,values,n=n)
   empirical.append(dict(unit=unit,states=len(rr),n=n,minimum_residual_lower=m['lower'],minimum_residual_upper=m['upper'],gap=m['gap'],v30_pairwise_max_lower=pair))
 write('joint_observation_audit.csv',empirical)
 # Each pair passes the broad beta_i in [0,1] test but no shared curve can pass log convexity.
 rr=[1.,1.25,1/.6];ss=[.2,.27,.3]
 pairpass=all(ap.interval_compatible(ss[i],ss[i],ss[j],ss[j],rr[j]/rr[i],rr[j]/rr[i]) for i in range(3) for j in range(i+1,3))
 assert pairpass and ss[1]>math.sqrt(ss[0]*ss[2])
 counter=margin(rr,ss,n=4096)
 info={'target_content':target,'Ce':Ce,'matched_mean_exposure':mu,'decisions':rows,'alias':alias,'capacity_exclusion':capacity,'pairwise_pass_joint_failure':{'ratios':rr,'residuals':ss,'all_pairs_pass':pairpass,'logconvex_upper':math.sqrt(.2*.3),'analytic_error_lower':(.27**2-.2*.3)/(2*.27+.2+.3),'margin':counter},'interpretation':'Conditional model-class consistency and constructive ambiguity; no new patients or independent clinical validation.'}
 (HERE/'results/summary.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
 print(json.dumps({'decisions':[(r['model'],r['information'],r['decision'],round(r['Ca_lower'],4),round(r['Ca_upper'],4)) for r in rows],'alias':alias[-1],'counterexample':info['pairwise_pass_joint_failure']},indent=2))

if __name__=='__main__':main()
