"""Close the capacity-flow law consistently; no new empirical fitting."""
from pathlib import Path
import sys,math,csv,json
import numpy as np
from scipy.optimize import brentq
R=Path(__file__).resolve().parents[1];O=R/'results';O.mkdir(exist_ok=True)
sys.path.insert(0,str(R.parent/'selection_extension/code'))
import analyze as a
v=a.v

def write(name,rows):
 with (O/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def state(q,m,ce,qref,sstar,beta):
 if not all(math.isfinite(z) for z in [q,m,ce,qref,sstar,beta]) or min(q,ce,qref)<=0 or m<0 or not 0<sstar<1 or not 0<=beta<=1:raise ValueError('outside declared domain')
 x=-math.log(sstar)*(q/qref)**(-beta);s=math.exp(-x);gap=-math.expm1(-x)
 cv=ce-m/(10*q*gap);deficit=m/(10*q*math.expm1(x));ca=ce-deficit
 beta_critical=gap/x
 dca=deficit/q*(1-beta/beta_critical)
 dd=10*ce-(m/q)*beta*x*math.exp(-x)/(gap*gap)
 return dict(Q=q,x=x,s=s,Ca=ca,Cv=cv,DO2=10*q*ca,feasible=cv>=0,beta_critical=beta_critical,dCa_dQ=dca,dDO2_dQ=dd)

def recovery(beta,ce):
 q0=5.;m=250.;ref=state(q0,m,ce,q0,.15,0)
 def f(r,key):return state(q0*r,m,ce,q0,.2,beta)[key]-ref[key]
 def root(key):
  if key=='Ca' and beta==1:return math.inf
  hi=2.
  while f(hi,key)<0 and hi<1e8:hi*=2
  if f(hi,key)<0:raise RuntimeError('unresolved recovery root')
  return brentq(lambda r:f(r,key),1.,hi,xtol=1e-12)
 qa,qd=root('Ca'),root('DO2');eps=(ce-ref['Ca'])/ce
 return dict(beta=beta,Qa_ratio=qa,Qd_ratio=qd,arterial_root_in_display=qa<=2,delivery_root_in_display=qd<=2,s_at_arterial_root=state(qa*q0,m,ce,q0,.2,beta)['s'] if math.isfinite(qa) else math.nan,s_at_delivery_root=state(qd*q0,m,ce,q0,.2,beta)['s'],fixed_s_ratio_identity_residual=(qd-1)-eps*(qa-1) if math.isfinite(qa) else math.nan)

def main():
 ce=v.content(100,14);q0=5.;m=250.;ref=state(q0,m,ce,q0,.15,0);rows=[]
 for beta in [0.,.25,.5,.75,1.]:
  for r in np.linspace(.6,2,141):
   st=state(q0*float(r),m,ce,q0,.2,beta)
   rows.append(dict(beta=beta,Q_ratio=float(r),Ca_ratio=st['Ca']/ref['Ca'],DO2_ratio=st['DO2']/ref['DO2'],**st))
 write('flow_response.csv',rows)
 roots=[recovery(beta,ce) for beta in [0.,.25,.5,.75,1.]];write('recovery_thresholds.csv',roots)
 boundary=[dict(s=float(s),x=-math.log(float(s)),beta_critical=(1-float(s))/-math.log(float(s))) for s in np.linspace(.01,.8,160)]
 write('slope_boundary.csv',boundary)
 # Exact conditional sensitivity to beta at beta=0; measured Pv remains fixed.
 leverage=[]
 for unit,rr in a.unit_groups(a.load()).items():
  b=next(r for r in rr if r['baseline']);s0=v.fit_one('equilibrium_admixture',b);x0=-math.log(s0)
  for r in rr:
   if r['baseline'] or not r['is_primary_PH']:continue
   y=a.beta_pred(0,b,r);p50=r['p50'];n=r['hill'];hb=r['Hb_g_dl']
   sat=y**n/(y**n+p50**n);cp=1.34*hb*n*sat*(1-sat)/y+.0031
   grad=-(v.content(r['PA_mmHg'],hb,p50,n)-v.content(r['Pv_mmHg'],hb,p50,n))*s0*x0*math.log(r['flow']/b['flow'])/cp
   h=1e-5;fd=(a.beta_pred(h,b,r)-a.beta_pred(-h,b,r))/(2*h)
   leverage.append(dict(family=r['cohort_group'],unit=unit,condition=r['condition'],FiO2=r['FiO2'],Q_ratio=r['flow']/b['flow'],dPa_dbeta_at_zero=grad,finite_difference=fd,residual=grad-fd))
 write('conditional_scaling_sensitivity.csv',leverage)
 summary=dict(reference=dict(Q=q0,M=m,Ce=ce,s0=.15,sstar=.2,Ca0=ref['Ca'],DO20=ref['DO2']),beta_critical_s015=(1-.15)/-math.log(.15),beta_critical_s020=(1-.2)/-math.log(.2),recovery=roots,conditional_derivative_max_abs_error=max(abs(z['residual']) for z in leverage),flow_grid_feasible=all(z['feasible'] for z in rows),interpretation='Analytic and numerical consequences of the imposed capacity-flow law, not patient response estimates')
 (O/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=True),encoding='utf-8');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
