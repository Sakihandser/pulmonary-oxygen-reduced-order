"""Forward counterfactuals, exact compensation boundaries and transparent attribution."""
from pathlib import Path
import itertools,math,json
import numpy as np
import analyze as a
v=a.v;O=a.O
def odds(s):
 if not 0<=s<1:raise ValueError('s outside [0,1)')
 return s/(1-s)
def closed(q,m,s,ce):
 if q<=0 or m<0 or ce<=0:raise ValueError('nonphysical input')
 odds(s)
 cv=ce-m/(10*q*(1-s));ca=ce-odds(s)*m/(10*q)
 return dict(Ca=ca,Cv=cv,DO2=10*q*ca,feasible=cv>=0)
def flow_thresholds(q0,m0,s0,ce0,m1,s1,ce1):
 """Exact absolute Q thresholds, including impossible arterial target."""
 ca0=closed(q0,m0,s0,ce0)['Ca'];d0=closed(q0,m0,s0,ce0)['DO2']
 delta=ce1-ca0
 burden=odds(s1)*m1
 qa=(0. if delta>=0 else math.inf) if burden==0 else (burden/(10*delta) if delta>0 else math.inf)
 qd=(d0+odds(s1)*m1)/(10*ce1)
 qf=m1/(10*(1-s1)*ce1)
 return qa,qd,qf
def shapley(b,r,s):
 keys=['PA_mmHg','Pv_mmHg','Hb_g_dl'];contributions={k:0. for k in keys}
 for order in itertools.permutations(keys):
  cur=a.allowed_inputs(b);last=v.predict('equilibrium_admixture',s,cur)
  for k in order:
   cur[k]=r[k];now=v.predict('equilibrium_admixture',s,cur);contributions[k]+=(now-last)/6;last=now
 return contributions
def main():
 rows=a.load();primary=[r for r in rows if r['is_primary_PH'] and r['FiO2']<.3]
 # No-hyperoxia exponent tuning is a declared secondary domain analysis.
 folds=[];predictions=[];profile=[]
 for family in a.FAMILIES:
  train=[r for r in primary if r['cohort_group']!=family];test=[r for r in primary if r['cohort_group']==family]
  beta,mse=a.fit_beta(train);folds.append(dict(held_out=family,beta=beta,training_RMSE=math.sqrt(mse)))
  for unit,rr in a.unit_groups(test).items():
   b=next(r for r in rr if r['baseline'])
   for r in rr:
    if not r['baseline']:predictions.append(a.prediction_record(r,'flow_adaptive_nonhyperoxia',a.beta_pred(beta,b,r),beta,'leave_family_out_nonhyperoxia'))
 for beta in np.linspace(-1,2,121):profile.append(dict(beta=float(beta),RMSE=math.sqrt(a.loss(beta,primary))))
 a.write('nonhyperoxia_exponent_folds.csv',folds);a.write('nonhyperoxia_exponent_predictions.csv',predictions);a.write('nonhyperoxia_exponent_profile.csv',profile)
 nested_score=a.metric('primary_PH_nonhyperoxia',predictions,'flow_adaptive_nonhyperoxia')
 # Delete complete cohort families and individual catheterization units.
 base_predictions=[r for r in a.frozen(rows) if r['primary_PH'] and r['FiO2']<.3]
 deleted=[]
 for key in ['family','unit']:
  for target in sorted(set(r[key] for r in base_predictions)):
   rr=[r for r in base_predictions if r[key]!=target]
   for law in a.MODELS:deleted.append(dict(deletion_level=key,deleted=target,**a.metric('nonhyperoxia',rr,law)))
 a.write('group_deletion.csv',deleted)
 # Endpoint transformation: secondary physical-content loss, same raw predictions.
 pp=a.frozen(rows);by_id={(r['study'],r['calibration_unit'],r['condition']):r for r in rows};metric=[]
 for p in pp:
  r=by_id[p['study'],p['unit'],p['condition']]
  p['content_error']=v.content(p['predicted_Pa'],r['Hb_g_dl'])-v.content(p['observed_Pa'],r['Hb_g_dl']) if p['predicted_Pa']>0 else math.nan
  p['fractional_gap_error']=p['error']/max(1,r['PA_mmHg']-r['Pv_mmHg'])
 for scope,rr in a.scopes(pp).items():
  for law in a.MODELS:
   zz=[x for x in rr if x['model']==law];ff=sorted(set(x['family'] for x in zz))
   for endpoint in ['content_error','fractional_gap_error']:
    metric.append(dict(scope=scope,model=law,endpoint=endpoint,RMSE=math.sqrt(np.mean([np.mean([x[endpoint]**2 for x in zz if x['family']==f]) for f in ff]))))
 a.write('endpoint_scale_scores.csv',metric)
 # Counterfactual attribution uses a frozen baseline s; mismatch is explicit residual.
 decomposition=[];diagnostics=[]
 for unit,rr in a.unit_groups(rows).items():
  b=next(r for r in rr if r['baseline']);s=v.fit_one('equilibrium_admixture',b)
  for r in rr:
   se=v.fit_one('equilibrium_admixture',r)
   diagnostics.append(dict(study=r['study'],unit=unit,condition=r['condition'],baseline=r['baseline'],s_baseline=s,s_state_diagnostic=se,flow_ratio=r['flow']/b['flow'],capacity_ratio_content_equivalent=(r['flow']*-math.log(se))/(b['flow']*-math.log(s)),use='outcome_reconstruction_not_prediction'))
   if r['baseline'] or not r['is_primary_PH'] or r['FiO2']>.3:continue
   sh=shapley(b,r,s);yp=a.pred('equilibrium_admixture',s,r)
   decomposition.append(dict(study=r['study'],unit=unit,condition=r['condition'],observed_delta_Pa=r['Pa_mmHg']-b['Pa_mmHg'],predicted_delta_Pa=yp-b['Pa_mmHg'],alveolar_component=sh['PA_mmHg'],venous_component=sh['Pv_mmHg'],Hb_component=sh['Hb_g_dl'],unexplained_residual=r['Pa_mmHg']-yp,closure_error=sum(sh.values())-(yp-b['Pa_mmHg'])))
 a.write('state_diagnostics.csv',diagnostics);a.write('response_decomposition.csv',decomposition)
 # Explicit phase map: predeclared reference, not pseudo-patients.
 q0,m0,s0,hb,pa=5.,250.,.15,14.,100.
 ce=v.content(pa,hb);ref=closed(q0,m0,s0,ce);maprows=[];boundary=[]
 for si in np.linspace(.05,.35,121):
  qa,qd,qf=flow_thresholds(q0,m0,s0,ce,m0,float(si),ce)
  boundary.append(dict(s_after=float(si),s_before=s0,Q_ratio_Pa_preserved=qa/q0,Q_ratio_DO2_preserved=qd/q0,Q_ratio_feasible=qf/q0))
  for ratio in np.linspace(.6,2.6,121):
   state=closed(q0*ratio,m0,float(si),ce)
   if not state['feasible']:cat='infeasible';saturation=math.nan;p=math.nan
   else:
    p=v.invcontent(state['Ca'],hb,26.8,2.7);saturation=100*v.sat(p)
    cat=('Pa_up' if state['Ca']>=ref['Ca'] else 'Pa_down')+'__'+('DO2_up' if state['DO2']>=ref['DO2'] else 'DO2_down')
   maprows.append(dict(s_after=float(si),flow_ratio=float(ratio),Pa=p,Sa_percent=saturation,DO2=state['DO2'],DO2_ratio=state['DO2']/ref['DO2'],region=cat))
 a.write('compensation_boundaries.csv',boundary);a.write('compensation_map.csv',maprows)
 cases=[]
 for si,ratio,mratio in [(.2,1.2,1),(.2,1.5,1),(.25,1.2,1),(.25,2.,1),(.15,1.2,1.5)]:
  st=closed(q0*ratio,m0*mratio,si,ce);qa,qd,qf=flow_thresholds(q0,m0,s0,ce,m0*mratio,si,ce)
  p0=v.invcontent(ref['Ca'],hb,26.8,2.7);p1=v.invcontent(st['Ca'],hb,26.8,2.7)
  cases.append(dict(s_before=s0,s_after=si,flow_ratio=ratio,demand_ratio=mratio,Pa_before=p0,Pa_after=p1,Sa_before=100*v.sat(p0),Sa_after=100*v.sat(p1),DO2_before=ref['DO2'],DO2_after=st['DO2'],DO2_percent=100*(st['DO2']/ref['DO2']-1),Q_ratio_Pa_preserved=qa/q0,Q_ratio_DO2_preserved=qd/q0))
 a.write('compensation_examples.csv',cases)
 pressure=[]
 for ratio in [.6,.8,1.,1.2]:
  for si in [.15,.20,.25]:
   qa,qd,qf=flow_thresholds(q0,m0,s0,ce,m0,si,ce)
   pressure.append(dict(driving_pressure_ratio=ratio,s_after=si,max_PVR_ratio_for_Pa_preservation=ratio/(qa/q0),max_PVR_ratio_for_DO2_preservation=ratio/(qd/q0)))
 a.write('pressure_resistance_thresholds.csv',pressure)
 # Published hemodynamic changes, never converted into a causal PVR->Q law.
 changes=[]
 for unit,rr in a.unit_groups(rows).items():
  if rr[0]['source_set']!='new':continue
  b=next(r for r in rr if r['baseline'])
  for r in rr:
   if r['baseline']:continue
   pvr=lambda x:float(x['PVR_native'])/(80 if x['study']=='Melot1983' else 1)
   changes.append(dict(study=r['study'],unit=unit,condition=r['condition'],PVR_before=pvr(b),PVR_after=pvr(r),delta_PVR=pvr(r)-pvr(b),Pa_before=b['Pa_mmHg'],Pa_after=r['Pa_mmHg'],delta_Pa=r['Pa_mmHg']-b['Pa_mmHg'],Q_ratio=r['flow']/b['flow'],mPAP_before=b['mPAP_mmHg'],mPAP_after=r['mPAP_mmHg'],measured_DO2_ratio=(float(r['DO2_index_ml_min_m2'])/float(b['DO2_index_ml_min_m2'])) if r['DO2_index_ml_min_m2'] else ''))
 a.write('published_hemodynamic_changes.csv',changes)
 summary=dict(nonhyperoxia_nested_score=nested_score,nonhyperoxia_beta_folds=folds,decomposition_max_abs=max(abs(r['closure_error']) for r in decomposition),compensation_examples=cases,reference=dict(Q=q0,M=m0,s=s0,Hb=hb,PA=pa,Ce=ce),warning='Mechanistic boundaries are model consequences; source Fick Q and outcome-reconstructed s are not independent validation.')
 (O/'mechanism_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
