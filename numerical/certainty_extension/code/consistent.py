"""Continuous cell outer bounds and realizable finite-mixture witnesses.

The proof is mathematical; floating evaluations use a declared numerical guard.
No solver status alone is called a physiological certificate.
"""
import numpy as np
from scipy.optimize import linprog
GUARD=1e-8
OPTIONS={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9}

def validate(r,lo,hi):
 r=np.asarray(r,float);lo=np.asarray(lo,float);hi=np.asarray(hi,float)
 if r.ndim!=1 or r.size==0 or not (r.shape==lo.shape==hi.shape):raise ValueError('shape')
 if not np.all(np.isfinite(r)) or np.any(r<=0) or np.any(lo>hi) or np.any(lo<0) or np.any(hi>1):raise ValueError('domain')
 return r,lo,hi

def matrices(r,lo,hi,n=2048,xmax=32.,moment=None):
 r,lo,hi=validate(r,lo,hi)
 if n<2 or xmax<=1e-6:raise ValueError('support')
 x=np.linspace(1e-6,xmax,n+1)
 f=np.exp(-x[None,:]/r[:,None])
 # On each cell f(right)<=f(x)<=f(left).
 outer=np.vstack([f[:,1:],-f[:,:-1]])
 inner=np.vstack([f,-f]);b=np.r_[hi,-lo]
 if moment is not None:
  ml,mh=moment
  if not 0<=ml<=mh<=xmax:raise ValueError('moment')
  outer=np.vstack([outer,x[:-1],-x[1:]])
  inner=np.vstack([inner,x,-x]);b=np.r_[b,mh,-ml]
 return x,outer,inner,b

def solve(c,A,b):
 return linprog(c,A_ub=A,b_ub=b,A_eq=np.ones((1,len(c))),b_eq=[1.],bounds=(0,None),method='highs',options=OPTIONS)

def dual_bound(c,A,b,result):
 if not result.success:raise ValueError('no dual')
 lam=np.maximum(0.,-result.ineqlin.marginals)
 nu=float(np.min(c+A.T@lam))
 val=nu-float(lam@b)
 # Conservative declared floating guard scales with multiplier amplification.
 padding=GUARD*(1.+float(np.sum(lam)))
 return val-padding,padding

def relaxation(A,b):
 n=A.shape[1]
 AA=np.c_[A,-np.ones(len(b))]
 c=np.r_[np.zeros(n),1.]
 res=linprog(c,A_ub=AA,b_ub=b,A_eq=np.r_[np.ones(n),0.][None,:],b_eq=[1.],bounds=(0,None),method='highs',options=OPTIONS)
 if not res.success:raise RuntimeError(res.message)
 lam=np.maximum(0.,-res.ineqlin.marginals)
 lam/=max(1.,float(lam.sum()))
 val=float(np.min(A.T@lam)-lam@b)-GUARD*(1+lam.sum())
 return max(0.,val),float(res.fun),res

def margin(r,values,n=2048,xmax=32.):
 x,A,I,b=matrices(r,values,values,n,xmax)
 lower,_,_=relaxation(A,b);_,upper,res=relaxation(I,b)
 w=np.maximum(0,res.x[:-1]);w/=w.sum()
 actual=np.exp(-x[None,:]/np.array(r)[:,None])@w
 upper=max(float(np.max(abs(actual-np.array(values)))),upper)+GUARD
 return dict(lower=lower,upper=upper,gap=upper-lower)

def bounds(r,lo,hi,query,n=2048,xmax=32.,moment=None):
 if query<=0 or not np.isfinite(query):raise ValueError('query')
 x,A,I,b=matrices(r,lo,hi,n,xmax,moment)
 necessary,slack,_=relaxation(A,b)
 out={'n':n,'query':query,'infeasibility_lower':necessary,'outer_slack':slack}
 if necessary>GUARD:
  out['status']='incompatible';return out
 q=np.exp(-x/query)
 low=solve(q[1:],A,b);up=solve(-q[:-1],A,b)
 if not low.success or not up.success:
  out['status']='unresolved';return out
 ol,pl=dual_bound(q[1:],A,b,low);nu,pu=dual_bound(-q[:-1],A,b,up)
 out.update(outer_lower=max(0.,ol),outer_upper=min(1.,-nu),dual_padding=max(pl,pu))
 witnesses=[]
 for direction,c in [('lower',q),('upper',-q)]:
  res=solve(c,I,b)
  if not res.success:continue
  w=np.maximum(0,res.x);w/=w.sum()
  violation=float(np.max(I@w-b))
  if violation>GUARD:continue
  take=w>1e-12
  witnesses.append({'direction':direction,'value':float(q@w),'x':x[take].tolist(),'weights':w[take].tolist(),'maximum_constraint_violation':violation})
 out['witnesses']=witnesses
 out['status']='bounded' if len(witnesses)==2 else 'unresolved'
 if len(witnesses)==2:
  out['inner_lower']=witnesses[0]['value'];out['inner_upper']=witnesses[1]['value']
  out['enclosure_gap']=max(out['inner_lower']-out['outer_lower'],out['outer_upper']-out['inner_upper'])
 return out

def content(s,Q,Ce,M):
 if not (0<=s<1 and Q>0 and Ce>0 and M>=0):raise ValueError('state')
 if M/(10*Q)>Ce*(1-s)+GUARD:raise ValueError('negative venous content')
 return Ce-M*s/(10*Q*(1-s))

def minimum_mean(r,lo,hi,n=8192,xmax=32.):
 r,lo,hi=validate(r,lo,hi)
 x=np.linspace(0.,xmax,n+1);f=np.exp(-x[None,:]/r[:,None]);b=np.r_[hi,-lo]
 # The extra cell [xmax,infinity) prevents an unmeasured support cap from
 # creating a false capacity exclusion. Its objective lower bound is xmax.
 A=np.vstack([np.c_[f[:,1:],np.zeros(len(r))],-np.c_[f[:,:-1],f[:,-1]]])
 I=np.vstack([f,-f]);c=np.r_[x[:-1],xmax]
 out=solve(c,A,b);inn=solve(x,I,b)
 if not out.success or not inn.success:raise RuntimeError('minimum capacity unresolved')
 lower,pad=dual_bound(c,A,b,out)
 w=np.maximum(0,inn.x);w/=w.sum()
 assert np.max(I@w-b)<=GUARD
 return {'mean_lower':lower,'mean_upper':float(x@w),'dual_padding':pad}

def classify(out,Q,Ce,M,target):
 if out['status']!='bounded':return out['status']
 low=content(out['outer_upper'],Q,Ce,M);high=content(out['outer_lower'],Q,Ce,M)
 if low>target+GUARD:return 'all_recover'
 if high<target-GUARD:return 'none_recover'
 a=content(out['inner_upper'],Q,Ce,M);b=content(out['inner_lower'],Q,Ce,M)
 if a<target-GUARD and b>target+GUARD:return 'opposite_witnesses'
 return 'unresolved'
