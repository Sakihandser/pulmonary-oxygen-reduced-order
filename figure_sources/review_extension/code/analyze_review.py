"""Matched scaling controls; exploratory stress tests, never patient uncertainty."""
from pathlib import Path
import sys,csv,json,itertools,math
import numpy as np
R=Path(__file__).resolve().parents[1];ROOT=R.parent;O=R/'results';O.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'selection_extension/code'))
import analyze as a
import mechanisms as m
v=a.v
LAWS=['saturation','content','pressure','constant_D']
def write(name,rows):
 with (O/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def matched(law,b,r,beta=0,flow_factor=1.):
 p=v.fit_one(law,b)
 x=a.allowed_inputs(r)
 # Algebraic cancellation: beta=0 does not even read changed-state Q.
 x['flow']=b['flow'] if beta==0 else b['flow']*(r['flow']*flow_factor/b['flow'])**beta
 return v.predict(law,p,x)
def predictions(rows,betas=(0,1),flow_factor=1.,shift=0.):
 out=[]
 for unit,rr in a.unit_groups(rows).items():
  b=next(r for r in rr if r['baseline'])
  for r in rr:
   if r['baseline'] or not r['is_primary_PH']:continue
   test=dict(r,p50=r['p50']+shift)
   for law,beta in itertools.product(LAWS,betas):
    y=matched(law,b,test,beta,flow_factor)
    out.append(dict(law=law,beta=beta,**a.prediction_record(r,f'{law}_beta{beta}',y,v.fit_one(law,b))))
 return out
def scopes(pp):
 return {k:z for k,z in dict(air_hypoxia=[p for p in pp if p['FiO2']<.3],modern=[p for p in pp if p['FiO2']<.3 and p['family']!='Melot1983'],thermodilution=[p for p in pp if p['FiO2']<.3 and p['family']=='Zurich_early'],direct_Fick=[p for p in pp if p['family']=='Zurich2021_2023'],historical=[p for p in pp if p['family']=='Melot1983'],pure_oxygen=[p for p in pp if p['FiO2']>.9],inclusive=pp).items() if z}
def scores(pp):
 return [a.metric(scope,rr,model) for scope,rr in scopes(pp).items() for model in sorted(set(p['model'] for p in rr))]
def main():
 rows=a.load();pp=predictions(rows);write('matched_predictions.csv',pp);ss=scores(pp);write('matched_metrics.csv',ss)
 sens=[];winners=[]
 for p50,n,rq,hb in itertools.product([24.,26.8,30.],[2.5,2.7,3.],[.75,.8,.9,1.],[12.,15.,18.]):
  zz=predictions(a.load(p50,n,rq,hb),betas=(0,));metrics=scores(zz)
  sens += [dict(P50=p50,Hill=n,R=rq,missing_Hb=hb,**x) for x in metrics]
  for scope in scopes(zz):
   best=min([x for x in metrics if x['scope']==scope],key=lambda x:x['equal_family_RMSE'])
   winners.append(dict(P50=p50,Hill=n,R=rq,missing_Hb=hb,scope=scope,winner=best['model']))
 write('matched_sensitivity.csv',sens);write('matched_winners.csv',winners)
 deletes=[]
 for family in a.FAMILIES:
  for x in scores([p for p in pp if p['family']!=family]):deletes.append(dict(deleted_family=family,**x))
 write('matched_deletion.csv',deletes)
 stress=[]
 for factor in [.8,1.,1.2]:
  stress += [dict(flow_ratio_factor=factor,**x) for x in scores(predictions(rows,flow_factor=factor))]
 write('flow_stress.csv',stress)
 shift=[]
 for delta in [-2.,-1.,0.,1.,2.]:
  shift += [dict(changed_state_P50_shift=delta,**x) for x in scores(predictions(rows,shift=delta))]
 write('state_curve_stress.csv',shift)
 leverage=[]
 for unit,rr in a.unit_groups(rows).items():
  b=next(r for r in rr if r['baseline'])
  for r in rr:
   if r['baseline'] or not r['is_primary_PH']:continue
   for law in LAWS:
    eps=1e-5
    grad=(matched(law,b,r,1,math.exp(eps))-matched(law,b,r,1,math.exp(-eps)))/(2*eps)
    method={'Zurich_early':'thermodilution','Zurich2021_2023':'direct_Fick','Melot1983':'not_established'}[r['cohort_group']]
    leverage.append(dict(family=r['cohort_group'],unit=unit,condition=r['condition'],FiO2=r['FiO2'],Q_ratio=r['flow']/b['flow'],law=law,dPa_dlogQ=grad,flow_method=method))
 write('flow_leverage.csv',leverage)
 ce=v.content(100,14);bounds=[]
 for radius in [0.,.01,.02]:
  for s0,s1 in itertools.product(sorted(set([.15-radius,.15+radius])),sorted(set([.20-radius,.20+radius]))):
   qa,qd,qf=m.flow_thresholds(5,250,s0,ce,250,s1,ce)
   bounds.append(dict(radius=radius,s0=s0,s1=s1,Qa_ratio=qa/5,Qd_ratio=qd/5,Qf_ratio=qf/5,Pa_PVR_reduction_fixed_driving_pressure=1-5/qa))
 write('boundary_stress.csv',bounds)
 # Portable dimensionless rule at unchanged Ce, allowing both M and s to change.
 gain=[]
 for s0,M0 in itertools.product([.05,.15,.3],[150.,250.,400.]):
  state=m.closed(5,M0,s0,ce)
  epsilon=(ce-state['Ca'])/ce
  for s1,M1 in itertools.product([.1,.2,.35],[150.,250.,400.]):
   qa,qd,qf=m.flow_thresholds(5,M0,s0,ce,M1,s1,ce)
   ratio=m.odds(s1)*M1/(m.odds(s0)*M0)
   gain.append(dict(s0=s0,M0=M0,s1=s1,M1=M1,epsilon=epsilon,burden_ratio=ratio,Qa_ratio=qa/5,Qd_ratio=qd/5,identity_residual=qd/5-1-epsilon*(qa/5-1)))
 write('dimensionless_compensation.csv',gain)
 summary=dict(matched_metrics=[x for x in ss if x['scope'] in ['air_hypoxia','modern','historical']],
  air_hypoxia_wins={law:sum(x['winner']==f'{law}_beta0' and x['scope']=='air_hypoxia' for x in winners) for law in LAWS},
  boundary_stress=bounds,dimensionless_rule_max_residual=max(abs(x['identity_residual']) for x in gain),reference_epsilon=m.odds(.15)*250/(10*5*ce),interpretation='Exploratory comparison; no independent cohort or sampling uncertainty added.')
 (O/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
