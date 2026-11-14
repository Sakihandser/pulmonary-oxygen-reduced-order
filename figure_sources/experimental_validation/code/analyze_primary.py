"""Frozen subgroup descriptions and one-factor PH counterfactuals; no new fit."""
import csv,json,math
from validate import CONFIG,OUT,write,LAWS
from couple_network import hydraulic,base,pars
from primary_model import solve,delivery_derivatives


def main():
    ref=CONFIG['reference']; cf=CONFIG['counterfactuals']
    h=hydraulic(ref['PVR_WU']); q0=h['CO_L_min']; s=pars[CONFIG['primary_law']]
    inputs=dict(q=q0,hb=ref['Hb_g_dl'],pa=base['PA_mmHg'],demand=ref['VO2_ml_min'],admixture=s)
    sequence=[]
    for pvr in CONFIG['pvr_grid_WU']:
        state=hydraulic(pvr)
        for closure,q in [('inherited_cardiac_closure',state['CO_L_min']),('flow_clamped',q0)]:
            sequence.append(dict(PVR_WU=pvr,closure=closure,Q_L_min=q,
                                 mPAP_mmHg=8+pvr*q,**solve(**dict(inputs,q=q))))
    write('primary_flow_controls.csv',sequence)
    changes=[('reference',{}),('flow_decreased_25pct',dict(q=q0*cf['flow_fraction'])),
             ('Hb_10',dict(hb=cf['Hb_g_dl'])),('PA_150',dict(pa=cf['PA_mmHg'])),
             ('mixing_0_30',dict(admixture=cf['effective_admixture'])),
             ('demand_500',dict(demand=cf['VO2_ml_min']))]
    baseline=solve(**inputs);counter=[]
    for name,change in changes:
        a=dict(inputs,**change);o=solve(**a)
        assert o['feasible']
        counter.append(dict(scenario=name,Q_L_min=a['q'],Hb_g_dl=a['hb'],PA_mmHg=a['pa'],
                            VO2_ml_min=a['demand'],s_eff=a['admixture'],**o,
                            delta_Sa_pp=o['Sa_percent']-baseline['Sa_percent'],
                            delta_DO2_ml_min=o['DO2_ml_min']-baseline['DO2_ml_min'],
                            DO2_change_percent=100*(o['DO2_ml_min']/baseline['DO2_ml_min']-1)))
    write('primary_counterfactuals.csv',counter)
    derivatives=delivery_derivatives(**inputs)
    steps={'Q':.01,'Hb':.1,'PA':1.,'M':10.,'s_eff':.01}
    units={'Q':'mL_O2/min per L/min','Hb':'mL_O2/min per g/dL','PA':'mL_O2/min per mmHg',
           'M':'mL_O2/min per mL_O2/min','s_eff':'mL_O2/min per fraction'}
    write('primary_local_sensitivity.csv',[dict(parameter=k,derivative=v,unit=units[k],
        illustrative_step=steps[k],linearized_DO2_change_ml_min=v*steps[k]) for k,v in derivatives.items()])
    with (OUT/'predictions.csv').open(encoding='utf-8') as f:pred=list(csv.DictReader(f))
    external=[r for r in pred if r['protocol']=='Carta_frozen_external']
    subsets={'all_pulmonary_vascular':{'D01','D02','D03','D04','D05','D06','D07'},
             'historical_primary_PH':{'D01','D02','D03','D04'},'recurrent_emboli':{'D05','D06','D07'}}
    subgroup=[]
    for subset,ids in subsets.items():
        for law in LAWS:
            rr=[r for r in external if r['condition'] in ids and r['law']==law]
            assert len(rr)==len(ids),(subset,law,[r['condition'] for r in external])
            errors=[float(r['error_mmHg']) for r in rr]
            subgroup.append(dict(subset=subset,n=len(errors),law=law,RMSE_mmHg=math.sqrt(sum(x*x for x in errors)/len(errors)),
                                 MAE_mmHg=sum(abs(x) for x in errors)/len(errors),bias_mmHg=sum(errors)/len(errors)))
    write('primary_subgroup_metrics.csv',subgroup)
    summary=dict(primary_law=CONFIG['primary_law'],effective_admixture=s,reference_inputs=inputs,
                 reference_outputs=baseline,maximum_mass_balance_error=max(abs(r['mass_balance_error']) for r in sequence+counter),
                 exact_partial_DO2=derivatives,observations_added=0,
                 design='Descriptive subgroup audit and specified one-factor scenarios added after initial model selection; no new validation cohort or treatment effect.',
                 fixed_flow_result='At fixed Q, Hb, PA, M and s_eff, pressure and regional redistribution have no direct oxygen effect in this closure.',
                 feasibility='Q_min is a nonnegative-Cv algebraic boundary, not a safe flow, exercise limit or clinical cutoff.')
    (OUT/'primary_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
