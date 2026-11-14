"""Outcome-isolated response prediction and cohort-family model selection."""
from pathlib import Path
import sys,csv,json,itertools,math
import numpy as np
from scipy.optimize import minimize_scalar
R=Path(__file__).resolve().parents[1];O=R/'results';O.mkdir(exist_ok=True)
sys.path.insert(0,str(R.parent/'experimental_validation/code'))
import validate as v
MODELS=v.LAWS+['baseline_Pa','baseline_Aa']
FAMILIES=['Zurich_early','Zurich2021_2023','Melot1983']
KPA=7.500616827
def write(name,rows):
 with (O/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def read_csv(path):
 with path.open(encoding='utf-8') as f:return list(csv.DictReader(f))
def load(p50=26.8,hill=2.7,rq=.8,hb_missing=15.):
 rows=[]
 for raw in read_csv(R.parent/'experimental_validation/data/observations.csv'):
  if raw['study']=='Dantzker1979':continue
  r=dict(raw)
  for k in ['flow','Pv_mmHg','Pa_mmHg','PaCO2_mmHg','Hb_g_dl','FiO2','PB_mmHg','PH2O_mmHg','baseline']:r[k]=float(r[k])
  if raw['Hb_status'].startswith('assumed'):r['Hb_g_dl']=hb_missing
  r.update(calibration_unit=r['study'],cohort_group='Zurich_early' if r['study'] in ['Carta2022','Lichtblau2021'] else r['study'],source_set='previous',Pv_status='reported',granularity='condition_mean')
  rows.append(r)
 for raw in read_csv(R/'data/additional_observations.csv'):
  r=dict(raw);scale=KPA if r['pressure_unit']=='kPa' else 1
  for k in ['flow','Hb_g_dl','FiO2','PB_mmHg','baseline']:r[k]=float(r[k])
  r['Pa_mmHg']=float(r['Pa_native'])*scale;r['PaCO2_mmHg']=float(r['PaCO2_native'])*scale;r['PH2O_mmHg']=47.1
  r['Pv_mmHg']=float(r['Pv_native'])*scale if r['Pv_native'] else v.invsat(float(r['Sv_percent'])/100,p50,hill)
  if r['Hb_status'].startswith('assumed'):r['Hb_g_dl']=hb_missing
  r['source_set']='new';rows.append(r)
 for r in rows:
  r.update(p50=p50,hill=hill)
  r['PA_mmHg']=r['FiO2']*(r['PB_mmHg']-r['PH2O_mmHg'])-r['PaCO2_mmHg']*(r['FiO2']+(1-r['FiO2'])/rq)
  r['is_primary_PH']=r['cohort_group'] in FAMILIES and r['Pv_status']=='reported'
 return rows
def unit_groups(rows):
 return {u:[r for r in rows if r['calibration_unit']==u] for u in sorted(set(r['calibration_unit'] for r in rows))}
def allowed_inputs(r):
 return {k:r[k] for k in ['flow','Pv_mmHg','PA_mmHg','Hb_g_dl','p50','hill']}
def base_fit(law,b):
 if law=='baseline_Pa':return b['Pa_mmHg']
 if law=='baseline_Aa':return b['PA_mmHg']-b['Pa_mmHg']
 return v.fit_one(law,b)
def pred(law,p,r):
 x=allowed_inputs(r)
 if law=='baseline_Pa':return p
 if law=='baseline_Aa':return x['PA_mmHg']-p
 return v.predict(law,p,x)
def prediction_record(r,law,y,par,protocol='baseline_frozen'):
 return dict(study=r['study'],unit=r['calibration_unit'],family=r['cohort_group'],condition=r['condition'],source_set=r['source_set'],primary_PH=int(r['is_primary_PH']),Pv_status=r['Pv_status'],FiO2=r['FiO2'],protocol=protocol,model=law,parameter=par,observed_Pa=r['Pa_mmHg'],predicted_Pa=y,error=y-r['Pa_mmHg'],feasible=int(r['Pv_mmHg']<=y<=r['PA_mmHg']),Q=r['flow'],Pv=r['Pv_mmHg'],PA=r['PA_mmHg'])
def frozen(rows):
 out=[]
 for unit,rr in unit_groups(rows).items():
  b=next(r for r in rr if r['baseline'])
  for law in MODELS:
   p=base_fit(law,b)
   for r in rr:
    if not r['baseline']:out.append(prediction_record(r,law,pred(law,p,r),p))
 return out
def scopes(preds):
 return {
  'primary_PH_all':[r for r in preds if r['primary_PH']],
  'primary_PH_nonhyperoxia':[r for r in preds if r['primary_PH'] and r['FiO2']<.3],
  'new_direct_Pv':[r for r in preds if r['primary_PH'] and r['source_set']=='new'],
  'secondary_exercise':[r for r in preds if r['study']=='Lichtblau2025'],
  'external_COPD_IPF':[r for r in preds if r['family'] in ['Blanco2010','Blanco2011']],
 }
def metric(scope,rr,law):
 x=[r for r in rr if r['model']==law]
 fam={g:[r['error'] for r in x if r['family']==g] for g in sorted(set(r['family'] for r in x))}
 return dict(scope=scope,model=law,test_conditions=len(x),families=len(fam),equal_family_RMSE=math.sqrt(np.mean([np.mean(np.square(e)) for e in fam.values()])),equal_family_MAE=float(np.mean([np.mean(np.abs(e)) for e in fam.values()])),equal_family_bias=float(np.mean([np.mean(e) for e in fam.values()])),pooled_RMSE=math.sqrt(np.mean([r['error']**2 for r in x])),infeasible=sum(1-r['feasible'] for r in x))
def scores(preds):
 out=[]
 for scope,rr in scopes(preds).items():
  for law in sorted(set(r['model'] for r in rr)):
   out.append(metric(scope,rr,law))
 for study in sorted(set(r['study'] for r in preds)):
  rr=[r for r in preds if r['study']==study]
  for law in sorted(set(r['model'] for r in rr)):out.append(metric(study,rr,law))
 return out
def beta_pred(beta,b,r):
 s0=v.fit_one('equilibrium_admixture',b)
 s=math.exp(math.log(s0)*(b['flow']/r['flow'])**beta)
 return pred('equilibrium_admixture',s,r)
def loss(beta,rows):
 errors={}
 for unit,rr in unit_groups(rows).items():
  b=next(r for r in rr if r['baseline'])
  for r in rr:
   if not r['baseline']:errors.setdefault(r['cohort_group'],[]).append((beta_pred(beta,b,r)-r['Pa_mmHg'])**2)
 return float(np.mean([np.mean(x) for x in errors.values()]))
def fit_beta(rows):
 grid=np.linspace(-1,2,61);vals=[loss(b,rows) for b in grid];i=int(np.argmin(vals))
 lo,hi=grid[max(0,i-1)],grid[min(len(grid)-1,i+1)]
 res=minimize_scalar(lambda b:loss(b,rows),bounds=(lo,hi),method='bounded',options={'xatol':1e-9})
 candidates=[(float(res.x),float(res.fun)),(-1.,loss(-1,rows)),(2.,loss(2,rows))]
 return min(candidates,key=lambda z:z[1])
def nested(rows):
 primary=[r for r in rows if r['is_primary_PH']];out=[];folds=[];profile=[]
 for family in FAMILIES:
  train=[r for r in primary if r['cohort_group']!=family];test=[r for r in primary if r['cohort_group']==family]
  beta,mse=fit_beta(train)
  folds.append(dict(held_out_family=family,training_families=';'.join(sorted(set(r['cohort_group'] for r in train))),beta=beta,training_RMSE=math.sqrt(mse)))
  for unit,rr in unit_groups(test).items():
   b=next(r for r in rr if r['baseline'])
   for r in rr:
    if not r['baseline']:out.append(prediction_record(r,'flow_adaptive',beta_pred(beta,b,r),beta,'leave_cohort_family_out_beta'))
 for beta in np.linspace(-1,2,121):
  profile.append(dict(beta=float(beta),all_primary_RMSE=math.sqrt(loss(beta,primary)),**{f:math.sqrt(loss(beta,[r for r in primary if r['cohort_group']==f])) for f in FAMILIES}))
 return out,folds,profile
def main():
 rows=load();pp=frozen(rows);nested_pp,folds,profile=nested(rows)
 write('response_predictions.csv',pp+nested_pp);write('response_metrics.csv',scores(pp+nested_pp));write('flow_exponent_folds.csv',folds);write('flow_exponent_profile.csv',profile)
 inputrows=[{k:r[k] for k in ['study','calibration_unit','cohort_group','condition','baseline','source_set','Pv_status','flow','flow_unit','Hb_g_dl','FiO2','PaCO2_mmHg','PA_mmHg','Pv_mmHg','Pa_mmHg','is_primary_PH']} for r in rows]
 write('analysis_inputs.csv',inputrows)
 # Fixed candidate comparison over declared input scenarios, no outcome-based tuning.
 sens=[]
 for p50,hill,rq,hb in itertools.product([24.,26.8,30.],[2.5,2.7,3.],[.75,.8,.9,1.],[12.,15.,18.]):
  z=frozen(load(p50,hill,rq,hb))
  for scope,rr in scopes(z).items():
   for law in MODELS:
    mm=metric(scope,rr,law);sens.append(dict(P50=p50,hill=hill,R=rq,missing_Hb=hb,**mm))
 write('response_sensitivity.csv',sens)
 win=[]
 for scope in scopes(pp):
  ss=[r for r in sens if r['scope']==scope];keys=sorted(set((r['P50'],r['hill'],r['R'],r['missing_Hb']) for r in ss))
  counts={law:0 for law in MODELS}
  for key in keys:
   group=[r for r in ss if (r['P50'],r['hill'],r['R'],r['missing_Hb'])==key];counts[min(group,key=lambda r:r['equal_family_RMSE'])['model']]+=1
  win.extend(dict(scope=scope,model=law,wins=counts[law],scenarios=len(keys)) for law in MODELS)
 write('response_sensitivity_winners.csv',win)
 summary=dict(additional_records=15,additional_direct_Pv=11,new_unique_individual_cases=2,new_direct_test_conditions=7,secondary_exercise_conditions=3,total_response_conditions=len(pp)//len(MODELS),primary_PH_conditions=sum(r['primary_PH'] for r in pp)//len(MODELS),scenario_count=108,folds=folds,nominal_scores=[r for r in scores(pp+nested_pp) if r['scope'] in ['primary_PH_all','primary_PH_nonhyperoxia','new_direct_Pv','secondary_exercise']])
 (O/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
