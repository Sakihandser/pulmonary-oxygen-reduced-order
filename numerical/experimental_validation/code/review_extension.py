"""Post-selection robustness and observation-design diagnostics, no new observations.

Exact empirical bootstrap holds all fitted parameters and measured inputs fixed.
It is not a post-selection confidence procedure or population superiority test.
"""
import csv, itertools, json, math
import numpy as np
from validate import ROOT, OUT, LAWS, read_data, fit_one, predict, content, write
from primary_model import solve

GROUPS={'all_pulmonary_vascular':list(range(7)), 'historical_primary_PH':list(range(4)), 'recurrent_emboli':list(range(4,7))}

def compositions(total, length):
    if length==1:
        yield (total,); return
    for x in range(total+1):
        for tail in compositions(total-x,length-1): yield (x,)+tail

def empirical_weights(n):
    counts=np.array(list(compositions(n,n)),dtype=int)
    weights=np.array([math.factorial(n)/math.prod(math.factorial(int(k)) for k in c)/n**n for c in counts])
    assert abs(weights.sum()-1)<1e-12
    return counts,weights

def weighted_quantile(x,w,p):
    order=np.argsort(x,kind='stable')
    return float(x[order[min(np.searchsorted(np.cumsum(w[order]),p),len(x)-1)]])

def bootstrap(errors):
    n=len(errors);counts,w=empirical_weights(n)
    rmse=np.sqrt(counts@np.square(errors)/n)
    wins=np.isclose(rmse,rmse.min(axis=1)[:,None],atol=1e-10,rtol=0)
    share=(w[:,None]*wins/wins.sum(axis=1)[:,None]).sum(axis=0)
    ref=LAWS.index('equilibrium_admixture')
    rows=[]
    for j,law in enumerate(LAWS):
        delta=rmse[:,ref]-rmse[:,j]
        rows.append(dict(law=law,n=n,RMSE_mmHg=float(np.sqrt(np.mean(errors[:,j]**2))),
                         admixture_minus_law_RMSE_mmHg=float(np.sqrt(np.mean(errors[:,ref]**2))-np.sqrt(np.mean(errors[:,j]**2))),
                         empirical_delta_q025=weighted_quantile(delta,w,.025),empirical_delta_q975=weighted_quantile(delta,w,.975),
                         empirical_best_share=float(share[j]),unique_count_vectors=len(w),ordered_resamples=n**n))
    return rows

def algebraic_jacobian(q,m,s,ce):
    """Rows Ca/Ce,Cv/Ce,log(Q); columns log Q, log M, logit(s)."""
    b=m/(10*q*(1-s));a=s*b
    return np.array([[a,-a,-a],[b,-b,-s*b],[ce,0,0]])/ce

def recover(ca,cv,ce,q):
    if not 0<=cv<=ca<=ce or ce<=cv or q<=0: raise ValueError('Outside identifiable physical domain')
    return (ce-ca)/(ce-cv),10*q*(ca-cv)

def gas_boundary(r,rq):
    # Deliberately never reads arterial oxygen or A-a gradient.
    return r['FiO2']*(r['PB_mmHg']-r['PH2O_mmHg'])-r['PaCO2_mmHg']*(r['FiO2']+(1-r['FiO2'])/rq)

def main():
    with (OUT/'predictions.csv').open(encoding='utf-8') as f:rr=list(csv.DictReader(f))
    ext=[r for r in rr if r['protocol']=='Carta_frozen_external']
    errors=np.array([[float(next(r['error_mmHg'] for r in ext if r['condition']==f'D{i:02d}' and r['law']==law)) for law in LAWS] for i in range(1,8)])
    boot=[];deletions=[]
    for name,idx in GROUPS.items():
        boot.extend(dict(subset=name,**r) for r in bootstrap(errors[idx]))
        for deleted in idx:
            remain=[i for i in idx if i!=deleted];rmse=np.sqrt(np.mean(errors[remain]**2,axis=0))
            for j,law in enumerate(LAWS):deletions.append(dict(subset=name,deleted=f'D{deleted+1:02d}',law=law,n_remaining=len(remain),RMSE_mmHg=float(rmse[j]),best=bool(abs(rmse[j]-min(rmse))<1e-10)))
    write('review_bootstrap.csv',boot);write('review_case_deletion.csv',deletions)
    gas=[];input_audit=[]
    for p50,hill,rq in itertools.product([24.,26.8,30.],[2.5,2.7,3.],[.75,.8,.9]):
        obs=read_data(p50,hill,rq);base=next(r for r in obs if r['study']=='Carta2022' and r['baseline'])
        external=[r for r in obs if r['study']=='Dantzker1979']
        for mode in ['source_reconstruction','gas_equation_only']:
            test=[dict(r,PA_mmHg=gas_boundary(r,rq)) if mode=='gas_equation_only' else r for r in external]
            for law in LAWS:
                par=fit_one(law,base)
                for group,idx in GROUPS.items():
                    e=np.array([predict(law,par,test[i])-test[i]['Pa_mmHg'] for i in idx])
                    gas.append(dict(boundary=mode,P50=p50,hill=hill,R=rq,subset=group,law=law,n=len(idx),RMSE_mmHg=float(np.sqrt(np.mean(e*e))),MAE_mmHg=float(np.mean(abs(e)))))
            if (p50,hill,rq)==(26.8,2.7,.8):
                input_audit.extend(dict(condition=r['condition'],boundary=mode,PA_mmHg=r['PA_mmHg']) for r in test)
    write('review_gas_sensitivity.csv',gas);write('review_gas_inputs.csv',input_audit)
    ref=json.loads((OUT/'primary_summary.json').read_text())['reference_inputs']
    state=solve(**ref);ce=state['equilibrium_content_ml_dl'];ca=state['Ca_ml_dl'];cv=state['Cv_ml_dl'];s=ref['admixture'];q=ref['q'];m=ref['demand']
    # Same arterial endpoint, with different s,M/Q; and same arterial/venous pair under common Q,M scaling.
    witnesses=[]
    for snew in [.08,s,.25]:
        mnew=10*q*(ce-ca)*(1-snew)/snew
        o=solve(**dict(ref,admixture=snew,demand=mnew))
        witnesses.append(dict(family='same_arterial',Q=q,M=mnew,s_eff=snew,**o))
    for factor in [.6,1.,1.4]:
        o=solve(**dict(ref,q=q*factor,demand=m*factor))
        witnesses.append(dict(family='same_arterial_and_venous',Q=q*factor,M=m*factor,s_eff=s,**o))
    assert all(r['feasible'] for r in witnesses)
    write('review_nonidentifiability_witness.csv',witnesses)
    j=algebraic_jacobian(q,m,s,ce);ranks=[]
    for design,matrix in [('arterial_only',j[:1]),('arterial_and_independent_Q',j[[0,2]]),('arterial_and_venous',j[:2]),('arterial_venous_independent_Q',j)]:
        singular=np.linalg.svd(matrix,compute_uv=False)
        ranks.append(dict(design=design,unknowns=3,rank=int(np.linalg.matrix_rank(matrix,tol=1e-10)),singular_values=';'.join(format(v,'.16g') for v in singular),rank_tolerance=1e-10))
    write('review_observation_rank.csv',ranks)
    # Independent bounded content errors; no inferred probability distribution.
    bounds=[]
    for width in [.05,.1,.2]:
        for ea,ev,ee,eq in itertools.product([-width,width],[-width,width],[-width,width],[-.1,.1]):
            aa,vv,ec,qq=ca+ea,cv+ev,ce+ee,q*(1+eq)
            ss,mm=recover(aa,vv,ec,qq)
            bounds.append(dict(content_halfwidth_ml_dl=width,flow_relative_error=eq,error_Ca=ea,error_Cv=ev,error_Ce=ee,s_eff=ss,M_ml_min=mm,DO2_ml_min=10*qq*aa))
    write('review_measurement_bounds.csv',bounds)
    summary=dict(added_observations=0,exact_resampling=True,calibration_refitted_in_bootstrap=False,
                 design='Retrospective post-selection diagnostics; empirical resampling quantiles have no asserted population coverage.',
                 rank_assumptions='Known Hb, blood curve and PA; interior s and positive M,Q. Independent Q adds rank; Fick-derived Q from the same oxygen observations does not.',
                 reference=dict(q=q,m=m,s=s,ce=ce,ca=ca,cv=cv),bootstrap=boot,ranks=ranks,
                 maximum_same_arterial_Pa_difference=max(abs(r['Pa_mmHg']-state['Pa_mmHg']) for r in witnesses),
                 maximum_same_pair_Pv_content_difference=max(abs(r['Cv_ml_dl']-cv) for r in witnesses if r['family']=='same_arterial_and_venous'))
    (OUT/'review_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
