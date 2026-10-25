"""Offline sensitivity experiments and numerical verification."""
from pathlib import Path
import csv, hashlib, importlib.util, itertools, json, platform, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from io import BytesIO
HERE=Path(__file__).resolve()
PHASE=HERE.parents[1]; ROOT=PHASE.parent
OUT=PHASE/'results'; OUT.mkdir(exist_ok=True)
FIG=PHASE/'figures'; FIG.mkdir(exist_ok=True)
PROTOCOL=json.loads((PHASE/'protocol.json').read_text())
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m; spec.loader.exec_module(m); return m
lung=load('sensitivity_hydraulic',ROOT/'hydraulic_network/code/build_lung_model.py')
obs=load('sensitivity_observation',ROOT/'observation_analysis/code/observations.py')
LOBES=list(lung.LOBE_W); W=np.array(list(lung.LOBE_W.values()))
Q0=PROTOCOL['reference']['Q']; PA=100.; checks=[]
def C(p,h=15):
    p=np.maximum(p,0); return 1.34*h*p**2.7/(p**2.7+26.8**2.7)+.0031*p
def S(p): return p**2.7/(p**2.7+26.8**2.7)
def CP(p,h=15):
    return 1.34*h*2.7*p**1.7*26.8**2.7/(p**2.7+26.8**2.7)**2+.0031
def inv(c,h=15): return brentq(lambda p:C(p,h)-c,0,200,xtol=1e-11)
def end(mode,q,share,pv=40,h=15,D=0,ratio=1,vc=86,tau=.3):
    if mode in ['saturation','content']:
        expo=.06*vc*share/(q*tau)*ratio
        if mode=='content': return C(PA,h)-(C(PA,h)-C(pv,h))*np.exp(-expo)
        s=S(PA)-(S(PA)-S(pv))*np.exp(-expo)
        return C(26.8*(s/(1-s))**(1/2.7),h)
    if D==0: return float(C(pv,h))
    sol=solve_ivp(lambda x,y:[D*share/(10*q)*(PA-y[0])/CP(max(y[0],0),h)],(0,1),[pv],method='DOP853',rtol=2e-10,atol=2e-11)
    if not sol.success: raise RuntimeError(sol.message)
    return float(C(min(PA,sol.y[0,-1]),h))
DM_RECOMPUTED=brentq(lambda d:end('D',Q0,1,D=d)-end('saturation',Q0,1),.1,40,xtol=1e-11)
# Use the full-precision reference calibration from the specified matrix.
DM=float(next(r['d_total'] for r in csv.DictReader((ROOT/'reference_matrix/results/ablation_matrix.csv').open(encoding='utf-8-sig')) if r['mode']=='conductance_matched'))
def gas(h=15,theta=.45):
    r=obs.invert_two_gases(16.4,48,theta)
    return r,obs.oxygen_estimate(r['dm'],r['vc'],1.8*h/15)['effective']
DG=gas()[1]
CASE=lung.build_cases(); HYD={}
def hydraulic(k):
    def get(a): return lung.solve(CASE,lobe_arterial_multiplier={l:a if l=='LUL' else a*k for l in LOBES})
    a=brentq(lambda a:get(a)['PVR_WU']-4.19,.001,400,xtol=1e-11); r=get(a)
    HYD[k]=(a,r,np.array([r['flow_fraction'][l] for l in LOBES])); return HYD[k]
for k in PROTOCOL['k']:hydraulic(k)
def run(k=25,allocation='fixed',mode='matched',h=15,pv=40,Q=Q0,ratio=1,D=None,M=None,scale=1,s=0,closure='capped',vc=86,tau=.3):
    f=HYD[k][2]; shares=W if allocation=='fixed' else f
    if D is None:D=DM if mode=='matched' else gas(h)[1]
    actual=mode if mode in ['saturation','content'] else 'D'
    def lungs(cv):
        pin=inv(cv,h)
        ends=np.array([end(actual,Q*x,y,pin,h,D,ratio,vc,tau) for x,y in zip(f,shares)])
        cm=float(f@ends); return (1-s)*cm+s*cv,cm,ends
    def use(ca):
        return sum(min(M*d,10*Q*fraction*min(1,e*scale)*ca) for fraction,d,e in lung.DEFAULT_ORGANS.values())
    feasible=True
    if M is None:
        cv=float(C(pv,h)); ca,cm,ends=lungs(cv); used=10*Q*(ca-cv); unmet=''; residual=0
    elif closure=='total':
        max_u=10*Q*lungs(0)[0]
        feasible=M<=max_u+1e-8
        if not feasible:
            return dict(k=k,allocation=allocation,mode=mode,hb=h,Q=Q,pv='',ratio=ratio,D=D,M=M,scale=scale,s=s,closure=closure,feasible=False,ca='',cv='',cm='',U='',unmet='',balance='',maximum_u=max_u,limited='')
        cv=brentq(lambda v:10*Q*(lungs(v)[0]-v)-M,0,C(PA,h),xtol=1e-10)
        ca,cm,ends=lungs(cv); used=M;unmet=0;residual=abs(10*Q*(ca-cv)-M)
    else:
        # Solve in venous content; the compartment sum always respects nonnegative Cv.
        cv=brentq(lambda v:lungs(v)[0]-use(lungs(v)[0])/(10*Q)-v,0,C(PA,h),xtol=1e-10)
        ca,cm,ends=lungs(cv); used=use(ca);unmet=M-used;residual=abs(10*Q*(ca-cv)-used)
    limited='' if M is None or closure=='total' else ';'.join(name for name,(fraction,d,e) in lung.DEFAULT_ORGANS.items() if 10*Q*fraction*min(1,e*scale)*ca<M*d-1e-8)
    if not (cv-1e-7<=ends.min()<=ends.max()<=C(PA,h)+1e-7):raise AssertionError('bounds')
    return dict(k=k,allocation=allocation,mode=mode,hb=h,Q=Q,pv=inv(cv,h),ratio=ratio,D=D if actual=='D' else '',M='' if M is None else M,scale=scale,s=s,closure='fixed_inlet' if M is None else closure,feasible=feasible,ca=ca,cv=cv,cm=cm,U=used,unmet=unmet,balance=residual,maximum_u='',limited=limited)
def write(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader();w.writerows(rows)
def check(name,ok,value):
    checks.append(dict(name=name,passed=bool(ok),value=value)); assert ok,name
def scan(name,configs):
    rr=[run(**x) for x in configs];write(name,rr);print(name,len(rr),flush=True);return rr
hydrows=[dict(k=k,a=a,PVR=r['PVR_WU'],Q=r['CO_L_min'],**{'f_'+l:f[i] for i,l in enumerate(LOBES)}) for k,(a,r,f) in HYD.items()]
write('hydraulic_scan.csv',hydrows)
check('matched_PVR_and_flow',max(max(abs(r['PVR']-4.19),abs(r['Q']-Q0)) for r in hydrows)<1e-8,hydrows)
ks=scan('k_scan.csv',[dict(k=k,mode=m,allocation=a) for k,m,a in itertools.product(PROTOCOL['k'],['saturation','content','matched','gas'],['fixed','flow'])])
ds=scan('D_scan.csv',[dict(k=k,D=d,mode='constant_D',allocation=a) for k,d,a in itertools.product([1,25],PROTOCOL['D'],['fixed','flow'])])
rat=scan('ratio_scan.csv',[dict(k=k,mode=m,ratio=r,allocation=a) for k,m,r,a in itertools.product([1,25],['saturation','content'],PROTOCOL['ratio_multipliers'],['fixed','flow'])])
ratio_check=[]
for m in ['saturation','content']:
    x=run(mode=m);y=run(mode=m,vc=172,tau=.6);ratio_check.append(abs(x['ca']-y['ca']))
check('fixed_Vc_tau_ratio',max(ratio_check)<1e-12,max(ratio_check))
fb=scan('feedback_grid.csv',[dict(k=k,mode=m,M=M,scale=e,s=.03) for k,m,M,e in itertools.product([1,25],['matched','gas'],PROTOCOL['demands'],PROTOCOL['ceiling_scales'])])
mix=scan('mixing_scan.csv',[dict(k=k,mode=m,M=250,s=s) for k,m,s in itertools.product([1,25],['saturation','content','matched','gas'],PROTOCOL['mixing'])])
qscan=scan('fixed_Q_scan.csv',[dict(k=k,mode=m,Q=q) for k,m,q in itertools.product([1,25],['saturation','content','matched','gas'],[3,Q0,5])])
alt=scan('total_demand_closure.csv',[dict(k=k,mode=m,M=M,s=.03,closure='total') for k,m,M in itertools.product([1,25],['matched','gas'],PROTOCOL['demands'])])
# Breakpoint scan and exact compartment activation roots, with the same flow and capacity.
dense=scan('demand_path.csv',[dict(mode=m,M=float(M),s=.03) for m,M in itertools.product(['matched','gas'],np.arange(50,401,5))])
critical=[]
for mode in ['matched','gas']:
    for name,(f,d,e) in lung.DEFAULT_ORGANS.items():
        def gap(M):return 10*Q0*f*e*run(mode=mode,M=M,s=.03)['ca']-M*d
        root=brentq(gap,1,1000,xtol=1e-7)
        v=run(mode=mode,M=root,s=.03)
        critical.append(dict(mode=mode,compartment=name,critical_demand=root,ca=v['ca'],cv=v['cv'],critical_ca_per_M=d/(10*Q0*f*e)))
write('extraction_breakpoints.csv',critical)
# Multiple reference matches; each fitted D is frozen before the other inlet/Q states.
transfer=[]
for pref in PROTOCOL['reference_fit_Pv']:
    dfit=brentq(lambda d:end('D',Q0,1,pref,D=d)-end('saturation',Q0,1,pref),.01,50)
    for ptest,qtest in [(30,Q0),(40,Q0),(50,Q0),(40,3),(40,5)]:
        for mode in ['saturation','content','constant_D']:
            value=end(mode if mode!='constant_D' else 'D',qtest,1,ptest,D=dfit)
            ref=end('saturation',qtest,1,ptest)
            transfer.append(dict(fit_Pv=pref,test_Pv=ptest,test_Q=qtest,mode=mode,fitted_D=dfit,content=value,reference_content=ref,residual=value-ref,is_fit_state=ptest==pref and qtest==Q0))
write('cross_state_transfer.csv',transfer)
# Additional fair match: common local equivalent capacity at P=40, not endpoint.
deq=.6*86/.3*CP(40)/ (2.7*S(40)*(1-S(40))/40)*(S(PA)-S(40))/(PA-40)
matched_local=scan('local_capacity_match.csv',[dict(k=k,mode=m,D=deq) for k,m in itertools.product([1,25],['saturation','content','constant_D'])])
gasrows=[]
for theta in [.30,.45,1.00]:
    r,d=gas(13.8,theta)
    for k in [1,25]:
        for feedback in [False,True]:
            value=run(k=k,mode='gas',h=13.8,D=d,M=250 if feedback else None,s=.03 if feedback else 0)
            gasrows.append(dict(thetaCO=theta,reported_DLCO=16.4,reported_DLNO=48,alpha=1.97,thetaNO=4.5,inferred_DM=r['dm'],inferred_Vc=r['vc'],thetaO2=1.8*13.8/15,assumed_D=d,**value))
write('gas_convention_outputs.csv',gasrows)
# External kernel: source lines in Alvin PartialPressure.CalculatePressures, N=50.
# Its pressure recurrence is independent of our content-ODE implementation.
# beta=Patm*K*A/(V*delta); fit beta once at reference transit, then freeze.
N=50
target_c=end('saturation',Q0,1);target_p=inv(target_c)
def alvin(pv,Q,beta,n=N):
    dt=.06*86/(Q*n); p=pv;arr=[p]
    for i in range(1,n):
        gain=max(0,min((PA-p),dt*beta*(PA-p)))
        p+=gain;arr.append(p)
    return np.array(arr)
beta=brentq(lambda b:alvin(40,Q0,b)[-1]-target_p,.001,30)
bench=[];traj=[]
for pv,Q in [(40,Q0),(30,Q0),(50,Q0),(40,3),(40,5)]:
    curve=alvin(pv,Q,beta);truth=float(C(curve[-1]))
    for mode,d in [('saturation',DM),('content',DM),('matched',DM),('gas',DG)]:
        value=end(mode if mode in ['saturation','content'] else 'D',Q,1,pv,D=d)
        bench.append(dict(Pv=pv,Q=Q,mode=mode,reference_fit=(pv==40 and Q==Q0),alvin_content=truth,model_content=value,residual=value-truth,beta_s_inverse=beta,sections=N))
    for i,p in enumerate(curve):traj.append(dict(Pv=pv,Q=Q,node=i,x=i/N,alvin_PO2=p,alvin_content=C(p)))
write('external_alvin_benchmark.csv',bench);write('external_alvin_trajectories.csv',traj)
native=[]
for n in [50,100,500,1000]:
    # Same frozen continuous coefficient, source's last stored node at (n-1)/n.
    native.append(dict(sections=n,endpoint_PO2=alvin(40,Q0,beta,n)[-1],x_end=(n-1)/n))
write('external_grid_refinement.csv',native)
check('external_reference_match',max(abs(r['residual']) for r in bench if r['reference_fit'] and r['mode'] in ['saturation','matched'])<1e-7,beta)
# Structural log sensitivities, with explicit missing direction.
J=np.array([[1,2,0],[1,1,1]],float);J3=np.vstack([J,[1,4,0]])
write('information_singular_values.csv',[dict(mapping='v_m',rank=int(np.linalg.matrix_rank(J)),singular_values=';'.join(map(str,np.linalg.svd(J,compute_uv=False))),null_direction='-2;1;1'),dict(mapping='v_m_g',rank=int(np.linalg.matrix_rank(J3)),singular_values=';'.join(map(str,np.linalg.svd(J3,compute_uv=False))),null_direction='none')])
check('information_rank',np.linalg.matrix_rank(J)==2 and abs(np.linalg.det(J3))>0,float(np.linalg.det(J3)))
check('matrix_oxygen_balance',max(r['balance'] for r in fb+mix+dense)<1e-7,max(r['balance'] for r in fb+mix+dense))
check('flow_allocation_invariance',max(abs(next(r['ca'] for r in ks if r['k']==1 and r['mode']==m and r['allocation']=='flow')-r['ca']) for m in ['saturation','content','matched','gas'] for r in ks if r['mode']==m and r['allocation']=='flow')<1e-7,'all seven flow vectors')
check('constant_D_allocation_bound',all(next(x['ca'] for x in ds if x['k']==r['k'] and x['D']==r['D'] and x['allocation']=='flow')>=r['ca']-1e-7 for r in ds),'Jensen upper bound')
check('capped_reference_use',abs(next(r['U'] for r in fb if r['k']==25 and r['mode']=='matched' and r['M']==250 and r['scale']==1)-184.44990)<2e-5,'reference systemic use reproduced')
check('nonnegative_unmet',min(r['unmet'] for r in fb+mix+dense)>-1e-7,min(r['unmet'] for r in fb+mix+dense))
# Cross-check pressure-state and content-state endpoints for all 128 conditions.
arch=list(csv.DictReader((ROOT/'reference_matrix/results/ablation_matrix.csv').open(encoding='utf-8-sig')))
errors=[]
for a in arch:
    mode={'conductance_matched':'matched','conductance_conditional':'gas'}.get(a['mode'],a['mode'])
    r=run(k=int(float(a['k'])),allocation=a['allocation'],mode=mode,h=float(a['hb']),M=250 if int(a['feedback']) else None,s=.03*int(a['admixture']))
    errors.append(abs(r['ca']-float(a['ca_ml_dl'])))
check('reference_128_pressure_ODE_crosscheck',max(errors)<2e-7,max(errors))
native_path=OUT/'external_alvin_native.csv'
if native_path.exists():
    native_rows=list(csv.DictReader(native_path.open(encoding='utf-8-sig')))
    native_error=max(abs(float(r['native_PO2'])-t['alvin_PO2']) for r,t in zip(native_rows,traj))
    check('compiled_upstream_Csharp_port_agreement',len(native_rows)==len(traj)==250 and native_error<5e-5,native_error)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,3,figsize=(12,3.7),layout='constrained')
colors=['#255b83','#d58e32','#a54545','#327e65']
for m,col,label in zip(['saturation','content','matched','gas'],colors,['Saturation','Content','Matched D','Gas D']):
    r=[x for x in ks if x['mode']==m and x['allocation']=='fixed'];baseline=r[0]['ca']
    ax[0].plot([HYD[x['k']][2][3] for x in r],[baseline-x['ca'] for x in r],'o-',c=col,label=label)
    rr=[x for x in ks if x['mode']==m and x['allocation']=='flow']
    ax[0].plot([HYD[x['k']][2][3] for x in rr],[rr[0]['ca']-x['ca'] for x in rr],'--',c=col,alpha=.6)
ax[0].set(xlabel='LUL flow fraction',ylabel='Content loss from k=1 (mL/dL)',title='A  Flow concentration');ax[0].legend(fontsize=7)
for alloc,ls in [('fixed','-'),('flow','--')]:
    dvals=PROTOCOL['D'];loss=[next(r['ca'] for r in ds if r['k']==1 and r['D']==d and r['allocation']==alloc)-next(r['ca'] for r in ds if r['k']==25 and r['D']==d and r['allocation']==alloc) for d in dvals]
    ax[1].plot(dvals,loss,'o'+ls,label=alloc)
ax[1].set(xlabel='Total D (mL/min/mmHg)',ylabel='k=1 minus k=25 (mL/dL)',title='B  Same law different capacity');ax[1].legend(fontsize=7)
for m,col in [('matched',colors[2]),('gas',colors[3])]:
    rr=[x for x in dense if x['mode']==m]
    ax[2].plot([x['M'] for x in rr],[x['U'] for x in rr],c=col,label=m+' U')
ax[2].plot([50,400],[50,400],':',c='gray',label='specified demand')
ax[2].set(xlabel='Demand (mL/min)',ylabel='Actual use (mL/min)',title='C  Capped feedback at k=25');ax[2].legend(fontsize=7)
def save_figure(fig,stem):
    fig.savefig(FIG/f'{stem}.png',dpi=180)
    fig.savefig(FIG/f'{stem}.svg')
    buffer=BytesIO();fig.savefig(buffer,format='png',dpi=600);buffer.seek(0)
    with Image.open(buffer) as source:rgb=source.convert('RGB')
    with (FIG/f'{stem}.tif').open('wb') as target:rgb.save(target,format='TIFF',compression='tiff_lzw',dpi=(600,600))
save_figure(fig,'Fig5')
plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(9,3.5),layout='constrained')
for m,col in zip(['saturation','content','matched','gas'],colors):
    rr=[r for r in bench if r['mode']==m and not r['reference_fit']]
    ax[0].plot(range(4),[r['residual'] for r in rr],'o-',label=m,c=col)
ax[0].set(xticks=range(4),xticklabels=['Pv30','Pv50','Q3','Q5'],ylabel='Content residual (mL/dL)',title='A  External pressure-kernel qualification');ax[0].legend(fontsize=7)
for m,col in [('matched',colors[2]),('gas',colors[3])]:
    rr=[r for r in rat if r['mode']=='saturation' and r['k']==25 and r['allocation']=='fixed'] if m=='matched' else [r for r in rat if r['mode']=='content' and r['k']==25 and r['allocation']=='fixed']
    ax[1].plot([r['ratio'] for r in rr],[r['ca'] for r in rr],'o-',c=col,label='saturation' if m=='matched' else 'content')
ax[1].set(xlabel='Vc/tau relative to reference',ylabel='Content at k=25 (mL/dL)',title='B  Relaxation capacity sensitivity');ax[1].legend(fontsize=7)
save_figure(fig,'FigS1')
plt.close(fig)
summary=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,DM=DM,independently_refitted_DM=DM_RECOMPUTED,DG=DG,ratio=DG/DM,checks=checks,benchmark_beta=beta,patient_validation=False,external_benchmark_scope='ported Alvin pressure submodule only; harmonized blood curve; no Unity application run',result_rows={p.name:sum(1 for _ in p.open())-1 for p in OUT.glob('*.csv')})
(OUT/'run_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
