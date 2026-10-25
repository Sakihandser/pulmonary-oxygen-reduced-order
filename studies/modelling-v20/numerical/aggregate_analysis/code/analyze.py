"""Offline pilot calibration and falsification checks, with explicit data splits.

No parameters are fitted to oxygen saturation or to external studies. These
analyses compare group central values, not subject-level statistical inference.
"""
from pathlib import Path
import csv
import json
import platform
import hashlib
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from data_inputs import main as prepare, write_csv
from microcirculation import OBSERVATION_MATRIX, HYDRAULIC_ROW, NULL_DIRECTION

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'results'


def load_groups():
    with (ROOT/'data/published_groups.csv').open(encoding='utf-8-sig') as f:
        rows=list(csv.DictReader(f))
    for r in rows:
        for k in list(r):
            if k in ('q','pap','pw','dm','vc','dlco','dlno','pao2','paco2','sao2','pvr','rer','aado2') or k.endswith('_spread'):
                r[k]=float(r[k]) if r[k] else None
    return rows


def fit_model(train, predictor, eta=0.4):
    base,exercise=train
    def coordinate(r):
        return r['q'] if predictor=='flow' else r['pw']+eta*(r['pap']-r['pw'])
    x0=coordinate(base)
    dx=np.log(coordinate(exercise)/x0)
    return dict(predictor=predictor,eta=eta,reference=x0,
                intercepts={v:base[v] for v in ('vc','dm','dlco')},
                exponents={v:float(np.log(exercise[v]/base[v])/dx) for v in ('vc','dm','dlco')})


def predict(model,r,anchor=None):
    if model['predictor']=='pressure' and r.get('pap') is None:
        return None
    def coordinate(row):
        return row['q'] if model['predictor']=='flow' else row['pw']+model['eta']*(row['pap']-row['pw'])
    ref=coordinate(anchor) if anchor else model['reference']
    intercept=anchor if anchor else model['intercepts']
    return {v:intercept[v]*(coordinate(r)/ref)**model['exponents'][v] for v in ('vc','dm','dlco')}


def main():
    prepare();OUTPUT.mkdir(exist_ok=True)
    rows=load_groups()
    train=[r for r in rows if r['study']=='Fermoyle2021' and r['group']=='control' and r['stage'] in ('rest','20W')]
    assert [r['stage'] for r in train]==['rest','20W']
    flow=fit_model(train,'flow');pressure=fit_model(train,'pressure')
    beta_v,beta_d=flow['exponents']['vc'],flow['exponents']['dm']
    common=(beta_v+beta_d)/2
    coupled=dict(flow,exponents=dict(flow['exponents'],vc=common,dm=common))
    # This constrained comparator fits Vc/Dm with equal log residual weights.
    # Its DLCO curve remains separately fitted and is not a mechanistic prediction.
    models={'flow_separate':flow,'pressure_separate_eta0.4':pressure,'shared_volume_exchange':coupled}
    predictions=[]
    for name,model in models.items():
        for r in rows:
            if r['study']=='Fermoyle2021' and r['group']=='control' and r['stage'] in ('rest','20W'):
                role='calibration'
            elif r['study']=='Fermoyle2021' and r['group']=='control':role='within_cohort_stage_holdout'
            elif r['study']=='Olson2016':role='cross_paper_population_transfer_unverified_overlap'
            else:role='disease_domain_stress_not_PAH'
            pred=predict(model,r)
            if pred is None:continue
            for variable in ('vc','dm','dlco'):
                if name=='shared_volume_exchange' and variable=='dlco':continue
                predictions.append(dict(model=name,study=r['study'],group=r['group'],stage=r['stage'],
                    role=role,variable=variable,observed=r[variable],predicted=pred[variable],
                    error=pred[variable]-r[variable],relative_error_percent=100*(pred[variable]/r[variable]-1),
                    published_spread=r[variable+'_spread'],spread_type=r['spread_type']))
    write_csv(OUTPUT/'pilot_predictions.csv',predictions)

    # Anchored ratios separate scale transfer from response transfer. The anchor
    # is an external observation and must never be called an unconditioned prediction.
    ratios=[]
    for group in ('control','HFpEF'):
        cohort=[r for r in rows if r['study']=='Olson2016' and r['group']==group]
        anchor=cohort[0]
        for r in cohort[1:]:
            pred=predict(flow,r,anchor=anchor)
            for variable in ('vc','dm','dlco'):
                ratios.append(dict(study='Olson2016',group=group,stage=r['stage'],variable=variable,
                    observed_ratio=r[variable]/anchor[variable],predicted_ratio=pred[variable]/anchor[variable],
                    conditional_prediction=pred[variable],observed=r[variable],
                    relative_error_percent=100*(pred[variable]/r[variable]-1),
                    role='external_rest_anchored_conditional_response'))
    with (ROOT/'data/pressure_challenge_summaries.csv').open(encoding='utf-8-sig') as f:
        challenge=list(csv.DictReader(f))
    con={r['variable']:float(r['value']) for r in challenge if r['stage']=='CON'}
    cuff={r['variable']:float(r['value']) for r in challenge if r['stage']=='CUFF'}
    for variable in ('vc','dm','dlco'):
        pred=(cuff['q']/con['q'])**flow['exponents'][variable]
        ratios.append(dict(study='DSouza2025',group='healthy',stage='CUFF',variable=variable,
            observed_ratio=cuff[variable]/con[variable],predicted_ratio=pred,
            conditional_prediction=con[variable]*pred,observed=cuff[variable],
            relative_error_percent=100*(con[variable]*pred/cuff[variable]-1),
            role='abstract_summary_challenge_subgroup_mismatch_not_joint_validation'))
    write_csv(OUTPUT/'external_response_checks.csv',ratios)

    profiles=[]
    for eta in np.linspace(0.,1.,101):
        mod=fit_model(train,'pressure',float(eta))
        peak=[r for r in rows if r['study']=='Fermoyle2021' and r['group']=='control' and r['stage']=='peak'][0]
        pred=predict(mod,peak)
        profiles.append(dict(eta=float(eta),training_log_error=0.,beta_v=mod['exponents']['vc'],
            beta_dm=mod['exponents']['dm'],peak_vc=pred['vc'],peak_dm=pred['dm']))
    write_csv(OUTPUT/'pressure_partition_equivalent_fits.csv',profiles)

    family=[]
    for b in np.linspace(0.,beta_v/2,11):
        a=beta_v-2*b;c=beta_d-beta_v+b
        family.append(dict(a_recruitment=float(a),b_radius=float(b),c_contact=float(c),
                           beta_v=float(a+2*b),beta_dm=float(a+b+c),beta_g=float(a+4*b)))
    write_csv(OUTPUT/'mechanism_equivalent_family.csv',family)
    matrix=np.vstack([OBSERVATION_MATRIX,HYDRAULIC_ROW])
    pc_cases=[]
    for label,ra,rv in [('upstream_dominant',4.,.5),('downstream_dominant',.5,4.)]:
        rc=.5;q=5.;pw=8.
        pc_cases.append(dict(label=label,ra=ra,rc=rc,rv=rv,q=q,pw=pw,pvr=ra+rc+rv,
                             pap=pw+q*(ra+rc+rv),pc_midpoint=pw+q*(rv+rc/2),
                             evidence='algebraic_counterexample_not_experimental_prediction'))
    write_csv(OUTPUT/'matched_pvr_pressure_counterexample.csv',pc_cases)

    summary=dict(date='2026-10-13',status='pilot_only_not_validated',calibration=flow,
        pressure_calibration=pressure,volume_exchange_common_exponent=common,
        observation_matrix=OBSERVATION_MATRIX.tolist(),observation_rank=int(np.linalg.matrix_rank(OBSERVATION_MATRIX)),
        augmented_capillary_conductance_rank=int(np.linalg.matrix_rank(matrix)),
        null_direction=NULL_DIRECTION.tolist(),hydraulic_change_along_null=float(HYDRAULIC_ROW@NULL_DIRECTION),
        dynamic_time_constants_identified=False,independent_oxygen_validation=False,
        subject_overlap_between_mayo_studies_verified=False)
    (OUTPUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')

    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.1),layout='constrained')
    control=[r for r in rows if r['study']=='Fermoyle2021' and r['group']=='control']
    grid=np.linspace(6.5,13,150)
    for ax,variable,label in zip(axes,('vc','dm'),('Capillary blood volume (mL)','Membrane conductance (mL/min/mmHg)')):
        ax.errorbar([r['q'] for r in control],[r[variable] for r in control],
                    yerr=[r[variable+'_spread'] for r in control],fmt='o',color='black',capsize=4,label='Published mean +/- SEM')
        for name,model,color,style in [('Separate response',flow,'#006d77','-'),('Shared Vc-Dm response',coupled,'#c44536','--')]:
            ax.plot(grid,flow['intercepts'][variable]*(grid/flow['reference'])**model['exponents'][variable],color=color,ls=style,label=name)
        ax.axvspan(10.8,13,color='#eeeeee',zorder=-1)
        ax.set(xlabel='Cardiac output (L/min)',ylabel=label)
    axes[0].legend(fontsize=8);axes[1].text(.97,.04,'Shaded range: peak stage held out\nTwo earlier stages used for fitting',transform=axes[1].transAxes,ha='right',fontsize=8)
    fig.savefig(OUTPUT/'pilot_calibration.png',dpi=190);fig.savefig(OUTPUT/'pilot_calibration.svg');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(10.5,4.1),layout='constrained')
    external=[r for r in rows if r['study']=='Olson2016' and r['group']=='control']
    for ax,variable,label in zip(axes,('vc','dm'),('Vc / resting Vc','Dm / resting Dm')):
        qrat=[r['q']/external[0]['q'] for r in external]
        ax.plot(qrat,[r[variable]/external[0][variable] for r in external],'ko-',label='Olson 2016 group ratios')
        x=np.linspace(1,3,100)
        ax.plot(x,x**flow['exponents'][variable],color='#006d77',ls='--',label='Frozen 2021 flow response')
        ax.set(xlabel='Q / resting Q',ylabel=label)
        ax.text(.03,.97,'Conditional on external resting values\nNo refit; cohort overlap unverified',transform=ax.transAxes,va='top',fontsize=8)
    axes[1].legend(loc='lower right',fontsize=8)
    fig.savefig(OUTPUT/'cross_paper_response.png',dpi=190);fig.savefig(OUTPUT/'cross_paper_response.svg');plt.close(fig)

    # No sampling confidence bands: paired covariance and individual data are absent.
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.1),layout='constrained')
    axes[0].plot([r['eta'] for r in profiles],[r['beta_dm'] for r in profiles],color='#006d77')
    axes[0].set(xlabel='Assumed pressure partition eta',ylabel='Equally fitting Dm pressure exponent')
    axes[0].text(.05,.9,'All curves fit the same two calibration stages\neta is assumed, not measured',transform=axes[0].transAxes,va='top',fontsize=8)
    axes[1].plot([r['b_radius'] for r in family],[r['beta_g'] for r in family],'o-',color='#c44536')
    axes[1].set(xlabel='Unidentified radius exponent b',ylabel='Capillary conductance exponent a + 4b')
    axes[1].text(.05,.87,'Identical Vc and Dm responses\nDifferent hydraulic responses',transform=axes[1].transAxes,fontsize=8)
    fig.savefig(OUTPUT/'identifiability.png',dpi=190);fig.savefig(OUTPUT/'identifiability.svg');plt.close(fig)
    files=[p for folder in ('code','data') for p in (ROOT/folder).glob('*') if p.is_file()]
    run=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
             matplotlib=matplotlib.__version__,random_sampling=False,
             inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
             fitted_studies=['Fermoyle2021 control rest and 20W'],
             oxygen_endpoints_used_for_fitting=False)
    (OUTPUT/'run_manifest.json').write_text(json.dumps(run,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
