"""Reproducible observation mapping and functional-identifiability analysis.

Prior data and results are read-only. All derived oxygen outputs are conditional
transformations of group means, not mean individual estimates or validation data.
"""
from pathlib import Path
import csv
import json
import hashlib
import platform
import xml.etree.ElementTree as ET
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from observations import (gas_conductance,invert_two_gases,oxygen_estimate,
    log_gas_jacobian,gain_threshold,geometry_observation_rows,GEOMETRY_NULL,
    pressure_driver_jacobian)

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT.parent/'aggregate_analysis'
RESULTS=ROOT/'results'
DATA=ROOT/'data'


def write_csv(path,rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)


def read_groups():
    with (PRIOR/'data/published_groups.csv').open(encoding='utf-8-sig') as f:
        rows=list(csv.DictReader(f))
    for row in rows:
        for k in ('q','dm','vc','dlco','dlno'):
            row[k]=float(row[k]) if row[k] else None
    return rows


def main():
    DATA.mkdir(exist_ok=True);RESULTS.mkdir(exist_ok=True)
    rows=read_groups()
    hb={('Fermoyle2021','control'):13.8,('Fermoyle2021','HFpEF'):11.9,
        ('Olson2016','control'):14.0,('Olson2016','HFpEF'):13.1}
    estimates=[]
    for row in rows:
        h=hb[(row['study'],row['group'])]
        estimate=oxygen_estimate(row['dm'],row['vc'],theta_o2=1.8*h/15)
        estimates.append(dict(study=row['study'],group=row['group'],stage=row['stage'],
            q=row['q'],dm_co=row['dm'],vc=row['vc'],baseline_group_hb=h,
            theta_o2_assumed=1.8*h/15,
            eD_Hb15=oxygen_estimate(row['dm'],row['vc'])['effective'],
            eD_baselineHb_assumption=estimate['effective'],
            eD_over_q=estimate['effective']/row['q'],
            membrane_fraction=estimate['membrane_resistance_fraction'],
            blood_fraction=estimate['blood_resistance_fraction'],
            level='nonlinear_transform_of_group_means_not_individual_mean',
            status='normoxic_operational_estimate_not_independent_measurement',
            Hb_status='baseline_group_mean_held_constant_not_measured_each_stage',
            source_doi=row['source_doi'],method_source_doi='10.1113/EP093578'))
    write_csv(RESULTS/'conditional_oxygen_estimates.csv',estimates)

    # Empirical group slopes from aggregate stages, not physiological elasticity estimates.
    prior=json.loads((PRIOR/'results/summary.json').read_text(encoding='utf-8'))
    beta_v=prior['calibration']['exponents']['vc'];beta_d=prior['calibration']['exponents']['dm']
    elasticity=[]
    for e in estimates:
        if e['study']=='Fermoyle2021' and e['group']=='control':
            beta=e['membrane_fraction']*beta_d+e['blood_fraction']*beta_v
            elasticity.append(dict(stage=e['stage'],beta_effective_local=beta,
                beta_effective_over_flow=beta-1,
                condition='constant_theta_and_local_empirical_power_laws_not_a_new_fit'))
    write_csv(RESULTS/'conditional_response_elasticities.csv',elasticity)

    # Examine the exact same reported DLCO and DLNO under unmeasured theta_CO.
    # This range is a diagnostic parameter profile, NOT a population prior or CI.
    base=next(r for r in rows if r['study']=='Fermoyle2021' and r['group']=='control' and r['stage']=='rest')
    profiles=[]
    for theta_co in np.linspace(.3,1.,141):
        out=invert_two_gases(base['dlco'],base['dlno'],theta_co)
        estimate=oxygen_estimate(out['dm'],out['vc'],theta_o2=1.8*13.8/15)
        jac=log_gas_jacobian(out['dm'],out['vc'],theta_co)
        inv=np.linalg.inv(jac)
        target_gradient=np.array([estimate['membrane_resistance_fraction'],estimate['blood_resistance_fraction']])@inv
        gas_eigen=np.linalg.eigvalsh(inv@inv.T)
        profiles.append(dict(theta_co_assumed=float(theta_co),dm=out['dm'],vc=out['vc'],
            eD=estimate['effective'],forward_dlco=gas_conductance(out['dm'],out['vc'],theta_co),
            forward_dlno=gas_conductance(out['dm'],out['vc'],4.5,1.97),
            forward_max_error=max(abs(gas_conductance(out['dm'],out['vc'],theta_co)-base['dlco']),
                abs(gas_conductance(out['dm'],out['vc'],4.5,1.97)-base['dlno'])),
            gas_log_condition_number=float(np.linalg.cond(jac)),
            log_dm_sensitivity_to_log_dlno=float(inv[0,1]),
            log_vc_sensitivity_to_log_dlno=float(inv[1,1]),
            eD_log_sensitivity_to_log_dlco=float(target_gradient[0]),
            eD_log_sensitivity_to_log_dlno=float(target_gradient[1]),
            variance_dm_equal_unit_independent_gas_errors=float((inv@inv.T)[0,0]),
            variance_vc_equal_unit_independent_gas_errors=float((inv@inv.T)[1,1]),
            covariance_dm_vc_equal_unit_independent_gas_errors=float((inv@inv.T)[0,1]),
            eD_variance_equal_unit_independent_gas_errors=float(target_gradient@target_gradient),
            label='diagnostic_profile_not_empirical_uncertainty'))
    write_csv(RESULTS/'same_gas_observation_profiles.csv',profiles)
    # Same-alpha finite/infinite comparison isolates theta_NO rather than
    # conflating coefficient protocols that also use different alpha.
    theta_center=1/(base['vc']*(1/base['dlco']-1/base['dm']))
    methods=[]
    for label,theta_no in [('finite_NO_4.5',4.5),('infinite_NO_same_alpha',float('inf'))]:
        out=invert_two_gases(base['dlco'],base['dlno'],theta_center,theta_no=theta_no)
        eff=oxygen_estimate(out['dm'],out['vc'],theta_o2=1.8*13.8/15)['effective']
        methods.append(dict(assumption=label,alpha=1.97,theta_co=theta_center,
            theta_co_status='aggregate_closure_identity_not_measured_coefficient',dm=out['dm'],vc=out['vc'],eD=eff,
            role='controlled_algebra_comparison_not_reanalysis_of_individual_data'))
    write_csv(RESULTS/'finite_infinite_same_alpha_comparison.csv',methods)

    # Determine which proposed observations add structural information.
    chosen=estimates[0]
    matrices=geometry_observation_rows(chosen['membrane_fraction'])
    rank_rows=[]
    for label,m in matrices.items():
        if label=='oxygen':continue
        sv=np.linalg.svd(m,compute_uv=False)
        rank_rows.append(dict(observations=label,rows=m.shape[0],rank=int(np.linalg.matrix_rank(m)),
            mechanism_parameters=3,null_dimension=3-int(np.linalg.matrix_rank(m)),
            smallest_singular_value=float(sv[-1]),
            condition='fixed_length_and_known_observation_coefficients_only'))
    write_csv(RESULTS/'additional_observation_rank.csv',rank_rows)

    # Transform 2025 abstract central values ONLY; no paired uncertainty or
    # clinical benefit is inferred from these ambiguous summary statistics.
    with (PRIOR/'data/pressure_challenge_summaries.csv').open(encoding='utf-8-sig') as f:
        challenge=list(csv.DictReader(f))
    con={r['variable']:float(r['value']) for r in challenge if r['stage']=='CON'}
    cuff={r['variable']:float(r['value']) for r in challenge if r['stage']=='CUFF'}
    m=1.23*con['dm'];b=1.8*con['vc'];r_m=cuff['dm']/con['dm'];r_b=cuff['vc']/con['vc']
    z=gain_threshold(r_m,r_b)
    start=oxygen_estimate(con['dm'],con['vc'])['effective']
    end=oxygen_estimate(cuff['dm'],cuff['vc'])['effective']
    gain=dict(dm_ratio=r_m,volume_ratio=r_b,flow_ratio=cuff['q']/con['q'],
        critical_blood_to_membrane_conductance_ratio=z,
        critical_theta_o2_for_fixed_abstract_inputs=z*m/con['vc'],
        assumed_theta_o2=1.8,assumed_hb_reference=15.,
        estimated_CON=start,estimated_CUFF=end,estimated_ratio=end/start,
        estimated_D_over_Q_ratio=(end/start)/(cuff['q']/con['q']),
        status='conditional_abstract_transform_not_independent_oxygen_validation')
    (RESULTS/'pressure_challenge_conditional_gain.json').write_text(json.dumps(gain,indent=2),encoding='utf-8')

    # Record primary-source observation coefficients; no raw pressure value is
    # substituted from arterial PaO2 for the unreported diffusion-test PAO2.
    constants=[dict(name='NO_CO_membrane_ratio',value=1.97,unit='dimensionless',source='10.1183/13993003.00962-2016',status='method_constant_not_fitted'),
        dict(name='theta_NO',value=4.5,unit='mL_NO/mL_blood/min/mmHg',source='10.1183/13993003.00962-2016',status='method_constant_not_fitted'),
        dict(name='O2_CO_membrane_ratio',value=1.23,unit='dimensionless',source='10.1113/EP093578',status='operational_conversion'),
        dict(name='normoxic_theta_O2_reference',value=1.8,unit='mL_O2/mL_blood/min/mmHg',source='10.1113/EP093578',status='conditional_reference_not_hypoxaemic_law'),
        dict(name='multi_FiO2_theta_CO_a',value=.0058,unit='inverse_theta/mmHg',source='10.1183/23120541.00889-2025',status='method_choice'),
        dict(name='multi_FiO2_theta_CO_b',value=.73,unit='inverse_theta',source='10.1183/23120541.00889-2025',status='method_choice'),
        dict(name='dual_gas_theta_CO_a',value=.0062,unit='inverse_theta/mmHg',source='10.1183/23120541.00889-2025',status='method_choice'),
        dict(name='dual_gas_theta_CO_b',value=1.16,unit='inverse_theta',source='10.1183/23120541.00889-2025',status='method_choice')]
    write_csv(DATA/'observation_constants.csv',constants)

    # Preserve supplementary observations as independent *literature benchmarks*,
    # never feed the estimated oxygen values back as measured validation targets.
    # Published table transcriptions are supplied as factual input data.
    # Full publisher XML is unnecessary for offline numerical reproduction.
    with (DATA/'published_operational_oxygen_benchmark.csv').open(encoding='utf-8-sig') as f:
        benchmark=list(csv.DictReader(f))
    assert benchmark and all(r['source_doi']=='10.1113/EP093578' for r in benchmark)

    # Directly transcribed selected primary-source table cells. Do not multiply
    # mean(DLCO/Q) by mean(Q) to manufacture a mean DLCO measurement.
    exercise=[]
    for group,stage,q,qsd,ratio,ratiosd,spo,spod in [
        ('UT','rest',6.20,1.35,4.81,.82,98,1),
        ('UT','90percent_peak_VO2',19.23,6.27,2.65,.79,97,1),
        ('T','rest',6.96,3.85,5.80,2.05,98,1),
        ('T','90percent_peak_VO2',23.83,4.41,2.68,.65,96,1)]:
        exercise.append(dict(study='Ehnes2026',group=group,drug='placebo',stage=stage,
            cohort_n=12,q=q,q_SD=qsd,published_mean_dlco_over_q=ratio,
            dlco_over_q_SD=ratiosd,spo2_percent=spo,spo2_SD=spod,
            source_doi='10.1113/JP290438',location='Table 2',
            q_method_status='not_reextracted_here',
            endpoint='pulse_oximetry_not_arterial_cooximetry',
            role='exercise_reserve_benchmark_not_diffusion_validation'))
    write_csv(DATA/'Ehnes2026_selected_exercise_cells.csv',exercise)
    direct_response=[]
    for group in ('UT','T'):
        rest,peak=[e for e in exercise if e['group']==group]
        direct_response.append(dict(group=group,
            ratio_dlco_over_q=peak['published_mean_dlco_over_q']/rest['published_mean_dlco_over_q'],
            dlco_over_q_change_percent=100*(peak['published_mean_dlco_over_q']/rest['published_mean_dlco_over_q']-1),
            rest_spo2=rest['spo2_percent'],peak_spo2=peak['spo2_percent'],
            q_ratio=peak['q']/rest['q'],
            inference='descriptive_no_paired_statistical_test_and_not_DL_O2'))
    write_csv(RESULTS/'direct_exercise_reserve_benchmark.csv',direct_response)

    # Analytical protocols, not numerical patient scenarios or fitted responses.
    driver_jac=pressure_driver_jacobian(1.,1.,1.,1.)
    protocols=[('upstream_resistance_change_at_fixed_Q_Pw',[0.,0.,1.,0.,0.]),
               ('downstream_pressure_change_at_fixed_Q_R',[1.,0.,0.,0.,0.]),
               ('Ra_to_Rv_shift_at_fixed_Q_Pw_total_R',[0.,0.,-1.,0.,1.])]
    pressure_protocols=[]
    for label,delta in protocols:
        response=driver_jac@delta
        pressure_protocols.append(dict(protocol=label,
            PAP_change_per_unit=float(response[0]),Pc_change_per_unit=float(response[1]),
            Q=1.,normalization='unit_flow_algebra_only_not_patient_prediction',
            status='only_valid_if_other_inputs_are_controlled'))
    write_csv(RESULTS/'pressure_driver_discriminating_protocols.csv',pressure_protocols)

    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(10.7,4.1),layout='constrained')
    x=np.array([p['theta_co_assumed'] for p in profiles])
    for name,color,style in [('dm','#007f86','-'),('vc','#b9473e','--'),('eD','#4c5b70',':')]:
        central=profiles[30][name]
        axes[0].plot(x,[p[name]/central for p in profiles],color=color,ls=style,label=name)
    axes[0].set(xlabel='Assumed blood conductance for CO',ylabel='Estimate relative to theta_CO = 0.45')
    axes[0].legend(fontsize=9)
    axes[0].text(.97,.94,'Same DLCO and DLNO at every point\nProfile is not a confidence interval',transform=axes[0].transAxes,ha='right',va='top',fontsize=8)
    axes[1].plot(x,[p['log_dm_sensitivity_to_log_dlno'] for p in profiles],color='#007f86',label='Dm sensitivity')
    axes[1].plot(x,[p['log_vc_sensitivity_to_log_dlno'] for p in profiles],color='#b9473e',ls='--',label='Vc sensitivity')
    axes[1].axhline(0,color='gray',lw=.8)
    axes[1].set(xlabel='Assumed blood conductance for CO',ylabel='Log sensitivity to DLNO at fixed DLCO')
    axes[1].legend(fontsize=9)
    fig.savefig(RESULTS/'observation_dependence.png',dpi=190);fig.savefig(RESULTS/'observation_dependence.svg');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(10.7,4.1),layout='constrained')
    for study,color,style in [('Fermoyle2021','#007f86','o-'),('Olson2016','#b9473e','s--')]:
        values=[e for e in estimates if e['study']==study and e['group']=='control']
        q=np.array([e['q'] for e in values]);d=np.array([e['eD_baselineHb_assumption'] for e in values])
        axes[0].plot(q/q[0],d/d[0],style,color=color,label=study)
        axes[1].plot(q/q[0],(d/q)/(d[0]/q[0]),style,color=color,label=study)
    axes[0].set(xlabel='Q / resting Q',ylabel='Conditional eD_O2 / resting eD_O2')
    axes[1].set(xlabel='Q / resting Q',ylabel='Conditional (eD_O2 / Q) ratio')
    axes[0].legend(fontsize=8);axes[1].axhline(1,color='gray',lw=.8)
    axes[0].text(.03,.96,'Transformation of group means\nNot an independent oxygen measurement',transform=axes[0].transAxes,va='top',fontsize=8)
    fig.savefig(RESULTS/'conditional_functional_response.png',dpi=190);fig.savefig(RESULTS/'conditional_functional_response.svg');plt.close(fig)

    summary=dict(status='observation_mapping_and_conditional_function_analysis',
        geometry_null_oxygen_dot=float(matrices['oxygen']@GEOMETRY_NULL),
        ranks={r['observations']:r['rank'] for r in rank_rows},
        theta_profile={'min':.3,'max':1.,'empirical_prior':False,'confidence_interval':False},
        independent_oxygen_validation=False,
        published_operational_estimate_is_validation_target=False,
        pressure_challenge=gain)
    (RESULTS/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    files=[p for folder in ('sources','code','data','tests') for p in (ROOT/folder).glob('*') if p.is_file()]
    files.append(PRIOR/'data/published_groups.csv');files.append(PRIOR/'results/summary.json')
    manifest=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        matplotlib=matplotlib.__version__,new_fitting_to_oxygen=False,
        files={str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (RESULTS/'run_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
