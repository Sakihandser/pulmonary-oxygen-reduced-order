"""Sharp baseline-constrained bounds and matched-capacity counterexamples.

All simulations are prescribed model scenarios, not patient observations.
"""
from pathlib import Path
import csv, json, math, sys
import numpy as np
from scipy.optimize import brentq

R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R.parent / 'coupling_extension/code'))
import analyze_coupling as prior
Q0, M, S0, SSTAR = 5., 250., .15, .20
CE = prior.v.content(100., 14.)

def g(z):
    z = np.asarray(z, dtype=float)
    out = np.zeros_like(z)
    mask = z > 0
    out[mask] = -z[mask]*np.log(z[mask])
    return out

def validate(weights, z):
    w, z = np.asarray(weights, float), np.asarray(z, float)
    if w.ndim != 1 or z.shape != w.shape or not len(w):
        raise ValueError('matching nonempty one-dimensional vectors required')
    if not np.all(np.isfinite(w)) or not np.all(np.isfinite(z)) or np.any(w < 0) or abs(w.sum()-1) > 1e-12 or np.any(z < 0) or np.any(z > 1):
        raise ValueError('weights must sum to one; residuals must lie in [0,1]')
    return w, z

def content_at(s, r, ce=CE, m=M, q0=Q0):
    if not all(math.isfinite(v) for v in [s,r,ce,m,q0]) or not 0 <= s < 1 or min(r,ce,q0) <= 0 or m < 0:
        raise ValueError('invalid oxygen balance domain')
    q = q0*r
    cv = ce-m/(10*q*(1-s))
    ca = ce-m*s/(10*q*(1-s))
    return ca, cv, 10*q*ca

def state(weights, z, r, beta, ce=CE, m=M, q0=Q0):
    w, z = validate(weights, z)
    if not math.isfinite(beta) or not 0 <= beta <= 1 or not math.isfinite(r) or r <= 0:
        raise ValueError('r>0 and beta in [0,1] required')
    zr = z**(r**(-beta))
    s = float(w @ zr)
    ca, cv, delivery = content_at(s,r,ce,m,q0)
    h = float(w @ g(zr))
    q = q0*r
    dca = m/(10*q*q)*(s/(1-s)-beta*h/(1-s)**2)
    dd = 10*ce-(m/q)*beta*h/(1-s)**2
    return dict(s=s,h=h,Ca=ca,Cv=cv,DO2=delivery,dCa_dQ=dca,dDO2_dQ=dd,feasible=cv>=0)

def bounds(mean, r, beta, lower=0., upper=1., ce=CE, m=M, q0=Q0, target=None):
    if not all(math.isfinite(v) for v in [mean,r,beta,lower,upper,ce,m,q0]) or not 0 <= lower <= mean <= upper <= 1 or not 0 < mean < 1 or r < 1 or not 0 <= beta <= 1 or min(ce,q0)<=0 or m<0:
        raise ValueError('invalid bound domain')
    a = r**(-beta)
    if lower == upper:
        low = high = mean**a
        hlow = hhigh = float(g(mean))
    else:
        w = (mean-lower)/(upper-lower)
        low = (1-w)*lower**a+w*upper**a
        high = mean**a
        hlow = float((1-w)*g(lower)+w*g(upper))
        hhigh = float(g(mean))
    sf = 1-m/(10*q0*r*ce)
    if low > sf+1e-13:
        return dict(s_lower=low,s_upper=high,feasible_s_upper=sf,Ca_lower=math.nan,Ca_upper=math.nan,status='no_feasible_model',h_lower=hlow,h_upper=hhigh)
    hi = min(high,sf)
    ca_low = content_at(hi,r,ce,m,q0)[0]
    ca_high = content_at(low,r,ce,m,q0)[0]
    if target is None:
        status='not_classified'
    elif ca_low >= target-1e-12:
        status='guaranteed'
    elif ca_high < target-1e-12:
        status='impossible'
    else:
        status='unresolved'
    return dict(s_lower=low,s_upper=high,feasible_s_upper=hi,Ca_lower=ca_low,Ca_upper=ca_high,status=status,h_lower=hlow,h_upper=hhigh)

def models():
    # A: finite unequal paths with the prescribed baseline residual.
    za=np.array([.001,.8]); wh=(SSTAR-za[0])/(za[1]-za[0]);wa=np.array([1-wh,wh])
    mean_x=float(wa @ (-np.log(za)))
    # B: match both E[z] and E[x], retaining finite positive capacities.
    xhi=20.
    def wb(xlo):return (xhi-mean_x)/(xhi-xlo)
    def residual(xlo):return wb(xlo)*math.exp(-xlo)+(1-wb(xlo))*math.exp(-xhi)-SSTAR
    xlo=brentq(residual,1e-9,mean_x,xtol=1e-14)
    b_w=np.array([wb(xlo),1-wb(xlo)]); b_z=np.exp(-np.array([xlo,xhi]))
    return {'finite_A':(wa,za),'matched_B':(b_w,b_z),'homogeneous':(np.array([1.]),np.array([SSTAR])),'ideal_limit':(np.array([1-SSTAR,SSTAR]),np.array([0.,1.]))}

def first_recovery(w,z,beta,target,key):
    grid=np.linspace(1,2,1001)
    vals=[state(w,z,float(r),beta)[key]-target for r in grid]
    for i in range(len(grid)-1):
        if vals[i] <= 0 <= vals[i+1]:
            root=brentq(lambda r:state(w,z,r,beta)[key]-target,float(grid[i]),float(grid[i+1]),xtol=1e-13)
            assert state(w,z,root,beta)['feasible']
            return root
    return math.nan

def write(name,rows):
    out=R/'results';out.mkdir(exist_ok=True)
    with (out/name).open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def main():
    ref=content_at(S0,1)
    rows=[];metadata=[];roots=[]
    for name,(w,z) in models().items():
        for i,(wi,zi) in enumerate(zip(w,z)):
            xi=-math.log(zi) if zi>0 else math.inf
            metadata.append(dict(model=name,path=i+1,weight=wi,z=zi,x=xi,K=Q0*wi*xi))
        for r in np.linspace(1,2,101):
            st=state(w,z,float(r),1)
            rows.append(dict(model=name,Q_ratio=r,Ca_ratio=st['Ca']/ref[0],DO2_ratio=st['DO2']/ref[2],**st))
        st=state(w,z,1,1)
        roots.append(dict(model=name,baseline_s=st['s'],total_K=sum(t['K'] for t in metadata if t['model']==name),baseline_h=st['h'],baseline_dCa_dQ=st['dCa_dQ'],Qa_ratio_in_1_to_2=first_recovery(w,z,1,ref[0],'Ca'),Qd_ratio_in_1_to_2=first_recovery(w,z,1,ref[2],'DO2'),Ca_at_twice_flow=state(w,z,2,1)['Ca']))
    write('matched_models.csv',metadata);write('heterogeneous_response.csv',rows);write('bounded_recovery.csv',roots)
    rowsb=[]
    for lo,hi in [(0.,1.),(.001,.8),(.05,.5)]:
        for beta in [0.,.25,.5,.75,1.]:
            for r in [1.,1.25,1.5,2.]:
                rowsb.append(dict(z_lower=lo,z_upper=hi,beta=beta,Q_ratio=r,**bounds(SSTAR,r,beta,lo,hi,target=ref[0])))
    write('sharp_bounds.csv',rowsb)
    local=[]
    for s in [.05,.15,.2,.4,.6,.8]:
        for beta in np.linspace(0,1,21):
            b=bounds(s,1,float(beta));scale=M/(10*Q0**2)
            local.append(dict(s=s,beta=beta,slope_lower=scale*(s/(1-s)-beta*b['h_upper']/(1-s)**2),slope_upper=scale*s/(1-s)))
    write('local_slope_bounds.csv',local)
    summary=dict(reference=dict(Q0=Q0,M=M,Ce=CE,s0=S0,sstar=SSTAR,Ca0=ref[0],DO20=ref[2]),models=roots,matched_total_K_residual=roots[0]['total_K']-roots[1]['total_K'],matched_s_residual=roots[0]['baseline_s']-roots[1]['baseline_s'],bounds_at_twice_flow=[r for r in rowsb if r['Q_ratio']==2],response_states=len(rows),all_response_states_feasible=all(r['feasible'] for r in rows),interpretation='Prescribed structural counterexample and sharp model-class bounds; no physiological distribution estimated')
    (R/'results/summary.json').write_text(json.dumps(summary,indent=2,allow_nan=True),encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
