"""Exact content-error propagation and conditional decision margins. No patient errors estimated."""
from pathlib import Path
import sys,csv,json,itertools
import numpy as np
H=Path(__file__).resolve().parents[1];R=H.parent
sys.path.insert(0,str(R/'certainty_extension/code'))
from consistent import bounds,content,classify
from analyze_certainty import structures

def residual(ce,ca,cv):
 if not ce>ca>cv:raise ValueError('strictly ordered contents required')
 return (ce-ca)/(ce-cv)

def box_interval(ce,ca,cv,ee,ea,ev):
 if min(ee,ea,ev)<0 or not ce-ee>ca+ea or not ca-ea>cv+ev:
  raise ValueError('box must preserve strict content ordering')
 return ((ce-ee-ca-ea)/(ce-ee-cv+ev),(ce+ee-ca+ea)/(ce+ee-cv-ev))

def write(name,rows):
 with (H/'results'/name).open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
 (H/'results').mkdir(exist_ok=True)
 old=json.loads((R/'certainty_extension/results/summary.json').read_text())
 ce=old['Ce'];target=old['target_content'];models,_=structures()
 precision=[];decisions=[];witnesses=[]
 for name,(x,w) in models.items():
  ratios=np.array([1.,1.25,1.5]);ss=np.exp(-x[None,:]/ratios[:,None])@w
  states=[]
  for r,s in zip(ratios,ss):
   ca=content(s,5*r,ce,250);cv=ca-250/(50*r);gap=ce-cv
   states.append((ce,ca,cv))
   precision.append(dict(model=name,r=r,s=s,Ce=ce,Ca=ca,Cv=cv,endpoint_gap=gap,equal_content_error_for_001=.001*gap/2,arterial_only_error_for_001=.001*gap))
  for e in [.001,.003,.005,.01,.02]:
   intervals=[box_interval(*st,e,e,e) for st in states]
   out=bounds(ratios,[v[0] for v in intervals],[v[1] for v in intervals],2.,n=32768)
   assert out['status']=='bounded'
   decisions.append(dict(model=name,content_halfwidth=e,cells=32768,s_lower=out['outer_lower'],s_upper=out['outer_upper'],Ca_lower=content(out['outer_upper'],10,ce,250),Ca_upper=content(out['outer_lower'],10,ce,250),target=target,decision=classify(out,10,ce,250,target),enclosure_gap=out['enclosure_gap']))
   for v in out['witnesses']:
    for xx,ww in zip(v['x'],v['weights']):witnesses.append(dict(model=name,content_halfwidth=e,direction=v['direction'],exposure=xx,weight=ww,constraint_violation=v['maximum_constraint_violation']))
 a=next(v for v in old['decisions'] if v['model']=='A' and v['information']=='two_extra')
 margin=a['Ca_lower']-target;deficit=ce-a['Ca_lower']
 frontier=[dict(Ce_down=b,M_relative_up_max=(margin-b)/deficit,margin=margin,deficit=deficit) for b in [0,.005,.01,margin]]
 write('measurement_precision.csv',precision);write('propagated_decisions.csv',decisions);write('measurement_witnesses.csv',witnesses);write('query_margin.csv',frontier)
 summary=dict(content_margin=margin,relative_consumption_margin=margin/deficit,precision=precision,decisions=decisions,assumptions='Content boxes are stipulated, not measured. Exact flow, fixed target, fixed structure class. Query Ce and M uncertainty is separate from observation error propagation. No probability coverage is assigned.')
 (H/'results/summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
