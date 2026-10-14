"""Measurement-anchored spatial identifiability witness, not clinical validation.

Reuse v08 archived matched-PVR flows and v09 gas observation/oxygen maps.
All new outputs stay outside the manuscript directories.
"""
import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parents[1]
LOBES = ['RUL', 'RML', 'RLL', 'LUL', 'LLL']


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def read_rows(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def main():
    inputs = {
        'network': ROOT / '03_第二轮_机制扫描/代码与数据/goal9_transit.csv',
        'published': ROOT / '06_第三轮_数据约束微循环/data/published_groups.csv',
        'oxygen_code': ROOT / '06_第三轮_数据约束微循环/code/microcirculation.py',
        'observation_code': ROOT / '07_第四轮_观测映射与功能预测/code/observations.py',
    }
    ox = module('joint_oxygen', inputs['oxygen_code'])
    obs = module('joint_observations', inputs['observation_code'])
    source = next(r for r in read_rows(inputs['published'])
                  if (r['study'], r['group'], r['stage']) == ('Fermoyle2021', 'control', 'rest'))
    network = read_rows(inputs['network'])
    selected = [next(r for r in network if float(r['contrast']) == contrast
                     and r['volume_rule'] == '解剖容积86') for contrast in [1., 25.]]
    # The published gas means and this simulated hydraulic operating point are
    # from different sources. Combining them is a witness, not patient calibration.
    constants = dict(theta_co=.45, theta_no=4.5, alpha=1.97,
                     lambda_o2_co=1.23, theta_o2=1.8 * 13.8 / 15.,
                     hb=13.8, alveolar_po2=100., inlet_po2=40., bypass_flow=0.)
    dlco, dlno = float(source['dlco']), float(source['dlno'])
    inv = obs.invert_two_gases(dlco, dlno, constants['theta_co'])
    dt = obs.oxygen_estimate(inv['dm'], inv['vc'], constants['theta_o2'])['effective']
    cv = ox.oxygen_content(constants['inlet_po2'], hb=constants['hb'])
    ceq = ox.oxygen_content(constants['alveolar_po2'], hb=constants['hb'])
    weights = np.array([1038.02, 486.49, 1204.04, 1278.15, 1053.67])
    weights /= weights.sum()
    results, regional = [], []
    checks = []

    def check(name, condition, evidence):
        checks.append(dict(name=name, passed=bool(condition), evidence=evidence))
        if not condition:
            raise AssertionError(name)

    for hydraulic in selected:
        contrast = float(hydraulic['contrast'])
        q = float(hydraulic['CO_L_min'])
        fractions = np.array([float(hydraulic['f_' + l]) for l in LOBES])
        check(f'flow_fractions_sum_{contrast:g}', abs(fractions.sum()-1) < 1e-12,
              float(fractions.sum()))
        for rule, shares in [('CT_fixed', weights), ('flow_proportional', fractions)]:
            dm, vc = shares * inv['dm'], shares * inv['vc']
            d = np.array([obs.oxygen_estimate(m, v, constants['theta_o2'])['effective']
                          for m, v in zip(dm, vc)])
            gas_co = sum(obs.gas_conductance(m, v, constants['theta_co'])
                         for m, v in zip(dm, vc))
            gas_no = sum(obs.gas_conductance(m, v, constants['theta_no'], constants['alpha'])
                         for m, v in zip(dm, vc))
            flows = q * fractions
            cout = np.array([ox.end_capillary_content(flow, conductance, 100., cv,
                                                      hb=constants['hb'])
                             for flow, conductance in zip(flows, d)])
            mixed = ox.mix_beds(flows, cout, cv)
            one_bed = ox.end_capillary_content(q, dt, 100., cv, hb=constants['hb'])
            check(f'whole_lung_gas_invariant_{contrast:g}_{rule}',
                  max(abs(gas_co-dlco), abs(gas_no-dlno), abs(d.sum()-dt),
                      abs(dm.sum()-inv['dm']), abs(vc.sum()-inv['vc'])) < 1e-10,
                  dict(dlco=gas_co, dlno=gas_no, effective_d=float(d.sum())))
            check(f'oxygen_mass_balance_{contrast:g}_{rule}',
                  abs(mixed['balance_residual']) < 1e-9, mixed['balance_residual'])
            check(f'content_bounds_{contrast:g}_{rule}',
                  np.all(cout >= cv-1e-8) and np.all(cout <= ceq+1e-8),
                  dict(inlet=cv, min=float(cout.min()), max=float(cout.max()), equilibrium=ceq))
            check(f'concavity_upper_bound_{contrast:g}_{rule}',
                  mixed['content'] <= one_bed+1e-7,
                  dict(mixed=mixed['content'], homogeneous=one_bed))
            if rule == 'flow_proportional':
                check(f'uniform_exposure_equality_{contrast:g}',
                      abs(mixed['content']-one_bed) < 1e-7,
                      mixed['content']-one_bed)
            results.append(dict(contrast=contrast, allocation=rule, pvr_wu=float(hydraulic['PVR_WU']),
                                q_l_min=q, dlco=gas_co, dlno=gas_no, dm_co=inv['dm'],
                                vc_ml=inv['vc'], effective_d_o2=dt,
                                mixed_content_ml_dl=mixed['content'], uptake_ml_min=mixed['uptake'],
                                deficit_to_homogeneous_ml_dl=one_bed-mixed['content']))
            for i, lobe in enumerate(LOBES):
                regional.append(dict(contrast=contrast, allocation=rule, lobe=lobe,
                                     flow_share=float(fractions[i]), capacity_share=float(shares[i]),
                                     relative_d_over_q=float(shares[i]/fractions[i]),
                                     dm_co=float(dm[i]), vc_ml=float(vc[i]),
                                     effective_d_o2=float(d[i]), end_content_ml_dl=float(cout[i])))

    # Test the exact spatial null over positive simplex points and coefficient
    # alternatives; this is algebraic robustness, not an uncertainty distribution.
    rng = np.random.default_rng(20261013)
    maximum_error = 0.
    for theta in [.3, .45, .6, 1.]:
        candidate = obs.invert_two_gases(dlco, dlno, theta)
        for shares in rng.dirichlet(np.ones(5), size=20):
            co = sum(obs.gas_conductance(s*candidate['dm'], s*candidate['vc'], theta) for s in shares)
            no = sum(obs.gas_conductance(s*candidate['dm'], s*candidate['vc'], 4.5, 1.97) for s in shares)
            maximum_error = max(maximum_error, abs(co-dlco), abs(no-dlno))
    check('spatial_null_across_coefficients', maximum_error < 1e-10, maximum_error)
    check('matched_hydraulics',
          abs(float(selected[0]['PVR_WU'])-float(selected[1]['PVR_WU'])) < 1e-10 and
          abs(float(selected[0]['CO_L_min'])-float(selected[1]['CO_L_min'])) < 1e-10,
          dict(pvr=[float(r['PVR_WU']) for r in selected], q=[float(r['CO_L_min']) for r in selected]))
    changed = [r for r in results if r['contrast'] == 25.]
    difference = abs(changed[0]['mixed_content_ml_dl']-changed[1]['mixed_content_ml_dl'])
    check('identical_global_observations_distinct_oxygen', difference > 1e-4, difference)
    concentrated = selected[1]
    q = float(concentrated['CO_L_min'])
    f = np.array([float(concentrated['f_' + l]) for l in LOBES])
    refined_contents = [ox.end_capillary_content(q*fi, dt*wi, 100., cv,
                                                hb=constants['hb'], rtol=1e-11)
                        for fi, wi in zip(f, weights)]
    refined_mixed = ox.mix_beds(q*f, refined_contents, cv)['content']
    refined_homogeneous = ox.end_capillary_content(q, dt, 100., cv,
                                                   hb=constants['hb'], rtol=1e-11)
    refined_gap = refined_homogeneous-refined_mixed
    check('oxygen_gap_tolerance_refinement', abs(refined_gap-difference) < 1e-7,
          dict(original=difference, refined=refined_gap, absolute_change=abs(refined_gap-difference)))

    result_dir = OUT / 'results'
    result_dir.mkdir(parents=True, exist_ok=True)
    for filename, rows in [('joint_results.csv', results), ('regional_witness.csv', regional)]:
        with (result_dir / filename).open('w', encoding='utf-8-sig', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    manifest = dict(status='structural_identifiability_witness_only', independent_oxygen_validation=False,
                    dynamically_calibrated_regional_coupling=False,
                    constants=constants, gas_anchor=source,
                    hydraulic_source='v08 archived five-lobe steady network; fixed exchange partition witnesses',
                    assumptions=['parallel additive gas conductances with shared coefficients',
                                 'common regional membrane/volume shares; hydraulic network unchanged',
                                 'conditional CO-to-O2 map and prescribed Hill curve',
                                 'common alveolar/inlet oxygen; no bypass; no venous feedback'],
                    spatial_unidentified_dimensions=4,
                    inputs={k: dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                            for k, p in inputs.items()},
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    checks=checks, results=results)
    (result_dir / 'audit_and_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2),
                                                       encoding='utf-8')
    print(json.dumps(dict(checks_passed=len(checks), conditional_gap_ml_dl=difference,
                          results=results), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
