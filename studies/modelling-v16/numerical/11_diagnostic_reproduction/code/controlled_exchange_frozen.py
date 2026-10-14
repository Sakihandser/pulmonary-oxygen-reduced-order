"""v11 controlled operator comparison. Synthetic diagnostics, no outcome validation."""
from pathlib import Path
import csv, hashlib, importlib.util, itertools, json, platform, sys
import numpy as np
from scipy.optimize import brentq
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parents[1] / 'results'
OUT.mkdir(parents=True, exist_ok=True)
def module(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec); sys.modules[name]=m; spec.loader.exec_module(m)
    return m
inputs={
 'legacy':ROOT/'02_baseline_network/code_and_data/build_lung_model.py',
 'hydraulic':ROOT/'03_mechanism_scans/code_and_data/goal9_transit.csv',
 'oxygen':ROOT/'06_data_constrained_microcirculation/code/microcirculation.py',
 'observation':ROOT/'07_observation_mapping/code/observations.py',
 'joint':ROOT/'08_joint_identifiability/results/joint_results.csv'}
lung=module('controlled_legacy',inputs['legacy'])
ox=module('controlled_oxygen',inputs['oxygen'])
obs=module('controlled_obs',inputs['observation'])
LOBES=['RUL','RML','RLL','LUL','LLL']
weights=np.array([lung.LOBE_W[l] for l in LOBES])
def read(p):
    return list(csv.DictReader(p.open(encoding='utf-8-sig',newline='')))
hyd={k:next(r for r in read(inputs['hydraulic']) if float(r['contrast'])==k and r['volume_rule']=='anatomical_volume_86') for k in [1.,25.]}
Q=float(hyd[1.]['CO_L_min']); VC=86.; TAU=.3; PA=100.
def sat(p): return p**2.7/(p**2.7+26.8**2.7)
def sprime(p): return 2.7*sat(p)*(1-sat(p))/p
def cprime(p,hb): return 1.34*hb*sprime(p)+.0031
def end(mode,flow,share,cv,hb,dtotal,rtol=1e-9):
    pv=ox.pressure_from_content(cv,hb); t=60*VC*share/(1000*flow)
    if mode=='saturation':
        s=sat(PA)-(sat(PA)-sat(pv))*np.exp(-t/TAU)
        p=26.8*(s/(1-s))**(1/2.7)
        return ox.oxygen_content(p,hb)
    if mode=='content':
        ceq=ox.oxygen_content(PA,hb)
        return ceq-(ceq-cv)*np.exp(-t/TAU)
    return ox.end_capillary_content(flow,dtotal*share,PA,cv,hb=hb,rtol=rtol)
cv15=ox.oxygen_content(40.,15.)
target=end('saturation',Q,1.,cv15,15.,0.)
DFIT=brentq(lambda d:end('conductance',Q,1.,cv15,15.,d)-target,.01,50.,xtol=1e-10)
inv=obs.invert_two_gases(16.4,48.,.45)
def dconditional(hb): return obs.oxygen_estimate(inv['dm'],inv['vc'],1.8*hb/15)['effective']
def evaluate(k,allocation,mode,hb,feedback,admixture,rtol=1e-9):
    row=hyd[float(k)]; f=np.array([float(row['f_'+l]) for l in LOBES]); flows=Q*f
    shares=weights if allocation=='fixed' else f
    d=DFIT if mode=='conductance_matched' else dconditional(hb)
    actual='conductance' if mode.startswith('conductance') else mode
    shunt=.03 if admixture else 0.
    def lungs(cv):
        ends=np.array([end(actual,q,s,cv,hb,d,rtol) for q,s in zip(flows,shares)])
        cm=float(f@ends); ca=(1-shunt)*cm+shunt*cv
        return ca,cm,ends
    def systemic(ca):
        used=sum(min(250*demand,emax*Q*flow*10*ca) for flow,demand,emax in lung.DEFAULT_ORGANS.values())
        return ca-used/(10*Q),used
    ceq=ox.oxygen_content(PA,hb)
    if feedback:
        ca=brentq(lambda a:lungs(systemic(a)[0])[0]-a,1e-5,ceq,xtol=1e-10)
        cv,used=systemic(ca); ca,cm,ends=lungs(cv)
        residual=10*Q*(ca-cv)-used
    else:
        cv=ox.oxygen_content(40,hb); ca,cm,ends=lungs(cv); used=None; residual=0.
    uptake=10*Q*(1-shunt)*float(f@(ends-cv))
    conservation=10*Q*(ca-cv)-uptake
    p=ox.pressure_from_content(ca,hb)
    return dict(k=k,allocation=allocation,mode=mode,hb=hb,feedback=int(feedback),admixture=int(admixture),
                d_total=d if actual=='conductance' else '',cv_ml_dl=cv,mixed_end_ml_dl=cm,
                ca_ml_dl=ca,sa_fraction=sat(p),do2_ml_min=10*Q*ca,uptake_ml_min=uptake,
                systemic_used_ml_min=used if used is not None else '',
                balance_error=max(abs(residual),abs(conservation)),
                minimum_end_ml_dl=float(ends.min()),maximum_end_ml_dl=float(ends.max()),equilibrium_ml_dl=ceq)
def write(name,rows):
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as file:
        w=csv.DictWriter(file,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
checks=[]
def check(name,condition,evidence):
    checks.append(dict(name=name,passed=bool(condition),evidence=evidence))
    if not condition: raise AssertionError(name)
case=lung.build_cases()
for k,row in hyd.items():
    def solve_scale(scale):
        return lung.solve(case,lobe_arterial_multiplier={l:scale if l=='LUL' else scale*k for l in LOBES})
    scale=brentq(lambda s:solve_scale(s)['PVR_WU']-4.19,.05,200.,xtol=1e-12)
    rebuilt=solve_scale(scale)
    errors=[abs(rebuilt['CO_L_min']-Q),abs(rebuilt['mPAP_mmHg']-float(row['mPAP_mmHg']))]
    errors.extend(abs(rebuilt['flow_fraction'][l]-float(row['f_'+l])) for l in LOBES)
    check('hydraulic_reconstruction_'+str(int(k)),max(errors)<1e-8,max(errors))
rows=[]
for k,a,m,h,f,s in itertools.product([1,25],['fixed','flow'],['saturation','content','conductance_matched','conductance_conditional'],[15.,13.8],[False,True],[False,True]):
    rows.append(evaluate(k,a,m,h,f,s))
write('ablation_matrix.csv',rows)
check('all_matrix_mass_balance',max(r['balance_error'] for r in rows)<1e-7,max(r['balance_error'] for r in rows))
check('all_matrix_bounds',all(r['cv_ml_dl']-1e-7<=r['minimum_end_ml_dl'] and r['maximum_end_ml_dl']<=r['equilibrium_ml_dl']+1e-7 for r in rows),'inlet <= regional end <= equilibrium')
for k in [1,25]:
    old=next(r for r in rows if (r['k'],r['allocation'],r['mode'],r['hb'],r['feedback'],r['admixture'])==(k,'fixed','saturation',15.,1,1))
    check('legacy_reproduction_'+str(k),abs(old['sa_fraction']-float(hyd[float(k)]['SaO2']))<1e-7,old['sa_fraction'])
    new=next(r for r in rows if (r['k'],r['allocation'],r['mode'],r['hb'],r['feedback'],r['admixture'])==(k,'fixed','conductance_conditional',13.8,0,0))
    archived=next(r for r in read(inputs['joint']) if float(r['contrast'])==k and r['allocation']=='CT_fixed')
    check('v10_reproduction_'+str(k),abs(new['mixed_end_ml_dl']-float(archived['mixed_content_ml_dl']))<1e-7,new['mixed_end_ml_dl'])
check('reference_endpoint_match',abs(end('conductance',Q,1,cv15,15,DFIT)-target)<1e-9,DFIT)
refined=evaluate(25,'fixed','conductance_conditional',13.8,False,False,1e-11)
check('integration_refinement',abs(refined['ca_ml_dl']-new['ca_ml_dl'])<1e-7,abs(refined['ca_ml_dl']-new['ca_ml_dl']))
check('zero_conductance',abs(end('conductance',Q,1,cv15,15,0)-cv15)<1e-12,'zero transfer returns inlet')
check('equilibrium_inlet',abs(end('conductance',Q,1,ox.oxygen_content(PA,15),15,DFIT)-ox.oxygen_content(PA,15))<1e-12,'equilibrium fixed point')
for mode in ['saturation','content','conductance_matched','conductance_conditional']:
    a=evaluate(1,'flow',mode,15,False,False);b=evaluate(25,'flow',mode,15,False,False)
    check('proportional_allocation_invariance_'+mode,abs(a['ca_ml_dl']-b['ca_ml_dl'])<1e-7,abs(a['ca_ml_dl']-b['ca_ml_dl']))
path=[('Legacy','saturation',15,True,True),('Content relaxation','content',15,True,True),
      ('Endpoint matched conductance','conductance_matched',15,True,True),
      ('Gas mapped conductance','conductance_conditional',15,True,True),
      ('Hb and blood coefficient','conductance_conditional',13.8,True,True),
      ('Fixed inlet','conductance_conditional',13.8,False,True),
      ('No admixture','conductance_conditional',13.8,False,False)]
pathrows=[]
for label,m,h,f,s in path:
    r=evaluate(25,'fixed',m,h,f,s); base=evaluate(1,'fixed',m,h,f,s)
    r.update(step=label,redistribution_content_loss=base['ca_ml_dl']-r['ca_ml_dl'],
             redistribution_sa_loss_pp=100*(base['sa_fraction']-r['sa_fraction']))
    pathrows.append(r)
for i,r in enumerate(pathrows):
    r['step_change_content_loss']=0 if i==0 else r['redistribution_content_loss']-pathrows[i-1]['redistribution_content_loss']
write('sequential_decomposition.csv',pathrows)
# Exact conductance represented by each prescribed relaxation trajectory.
profiles=[]
for p in np.linspace(40,99.9,180):
    ceq=ox.oxygen_content(PA,15); c=ox.oxygen_content(p,15)
    dc=.6*VC/TAU*(ceq-c)/(PA-p)
    ds=.6*VC/TAU*cprime(p,15)/sprime(p)*(sat(PA)-sat(p))/(PA-p)
    profiles.append(dict(po2=p,saturation_equivalent_d=ds,content_equivalent_d=dc,endpoint_matched_d=DFIT,gas_mapped_d=dconditional(15)))
write('equivalent_conductance.csv',profiles)
near=.6*VC/TAU*cprime(PA,15)
check('relaxation_equilibrium_limit',abs(profiles[-1]['content_equivalent_d']-near)/near<.003,near)
fig,axs=plt.subplots(1,2,figsize=(10,3.6),layout='constrained')
for key,label in [('saturation_equivalent_d','Saturation relaxation'),('content_equivalent_d','Content relaxation')]:
    axs[0].plot([r['po2'] for r in profiles],[r[key] for r in profiles],label=label)
axs[0].axhline(DFIT,color='gray',ls='--',label='Endpoint matched')
axs[0].axhline(dconditional(15),color='#a64b38',ls=':',label='Gas mapped')
axs[0].set(xlabel='Capillary PO2 (mmHg)',ylabel='Equivalent D (mL/min/mmHg)',title='A  Relaxation implies state dependent capacity')
axs[0].legend(fontsize=7,frameon=False)
axs[1].plot(range(7),[r['redistribution_content_loss'] for r in pathrows],'o-',color='#305a7d')
axs[1].set(xticks=range(7),xticklabels=['Legacy','Content','Matched D','Gas D','Hb','Inlet','Mix'],ylabel='k=1 minus k=25 content (mL/dL)',title='B  Sequential changes in redistribution loss')
axs[1].tick_params(axis='x',labelsize=8)
for ax in axs: ax.spines[['top','right']].set_visible(False)
fig.savefig(OUT/'operator_decomposition.png',dpi=200);fig.savefig(OUT/'operator_decomposition.svg');plt.close(fig)
# Quantify boundary interaction rather than assign a unique causal percentage.
interaction=[]
for hb,dmode in itertools.product([15.,13.8],['conductance_matched','conductance_conditional']):
    losses={}
    for feedback,admixture in itertools.product([0,1],[0,1]):
        pair=[next(r for r in rows if (r['k'],r['allocation'],r['mode'],r['hb'],r['feedback'],r['admixture'])==(k,'fixed',dmode,hb,feedback,admixture)) for k in [1,25]]
        losses[feedback,admixture]=pair[0]['ca_ml_dl']-pair[1]['ca_ml_dl']
    interaction.append(dict(hb=hb,mode=dmode,feedback_effect_no_admixture=losses[1,0]-losses[0,0],
      admixture_effect_no_feedback=losses[0,1]-losses[0,0],
      feedback_admixture_interaction=losses[1,1]-losses[1,0]-losses[0,1]+losses[0,0]))
write('boundary_interactions.csv',interaction)
manifest=dict(date='2026-10-14',status='controlled_synthetic_comparison',independent_validation=False,
  python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,
  q_l_min=Q,vc_relaxation_ml=VC,tau_s=TAU,d_endpoint_matched=DFIT,d_conditional_hb15=dconditional(15),
  d_conditional_hb13_8=dconditional(13.8),near_equilibrium_equivalent_d=near,
  diagnostic_fit='one homogeneous bed at Q=archived Q, PA=100, Pv=40, Hb=15; matches saturation endpoint only',
  capacity_comparison='same regional shares; relaxation Vc=86, conditional inversion Vc differs; capacity construction is an explicit factor',
  feedback='six legacy compartments with extraction ceilings and demand 250; admixture is postexchange content mixing',
  checks=checks,path=pathrows,interactions=interaction,
  inputs={key:dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for key,p in inputs.items()},
  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),matrix_rows=len(rows))
(OUT/'audit_and_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({key:manifest[key] for key in ['matrix_rows','d_endpoint_matched','d_conditional_hb15','near_equilibrium_equivalent_d','path','interactions','checks']},ensure_ascii=False,indent=2))
