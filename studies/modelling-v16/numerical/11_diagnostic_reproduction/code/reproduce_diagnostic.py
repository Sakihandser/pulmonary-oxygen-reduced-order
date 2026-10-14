"""Rerun 128 original conditions and independent endpoints without article text."""
from pathlib import Path
import csv,hashlib,importlib.util,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
PHASE=Path(__file__).resolve().parents[1]
OUT=PHASE/'results';OUT.mkdir(exist_ok=True)
def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
def read(p):
 with p.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def run():
 code=PHASE/'code/controlled_exchange_frozen.py'
 model=load('diagnostic_production',code)
 audit=load('independent_audit_functions',PHASE/'code/independent_endpoint_functions.py')
 new=read(OUT/'ablation_matrix.csv');old=read(ROOT/'09_controlled_exchange/results/ablation_matrix.csv')
 assert len(new)==len(old)==128
 maxdiff=0.;errors=[]
 for a,b in zip(old,new):
  for k in a:
   try:maxdiff=max(maxdiff,abs(float(a[k])-float(b[k])))
   except ValueError:assert a[k]==b[k]
  f=np.array([float(model.hyd[float(b['k'])]['f_'+l]) for l in model.LOBES]);shares=model.weights if b['allocation']=='fixed' else f
  cv,h=float(b['cv_ml_dl']),float(b['hb']);d=float(b['d_total']) if b['d_total'] else 0
  ends=np.array([audit.independent_end(b['mode'],model.Q*x,s,cv,h,d) for x,s in zip(f,shares)])
  cm=float(f@ends);ca=(1-.03*int(b['admixture']))*cm+.03*int(b['admixture'])*cv
  errors.append(ca-float(b['ca_ml_dl']))
 assert maxdiff<1e-7 and max(abs(e) for e in errors)<2e-7
 out=dict(matrix_rows=128,production_checks=len(model.checks),all_checks_pass=all(c['passed'] for c in model.checks),
          archived_max_difference=maxdiff,independent_max_error=max(abs(e) for e in errors),source_cell_audit='102 cells separately audited in preserved v12 archive',physiological_validation=False)
 (OUT/'diagnostic_rerun_audit.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
if __name__=='__main__':run()
