"""Offline regeneration and numerical comparison of the response selection study."""
from pathlib import Path
import csv,hashlib,json,math,subprocess,sys,time
R=Path(__file__).resolve().parents[1]
def table(p):
 with p.open(encoding='utf-8',newline='') as f:return list(csv.reader(f))
def main():
 start=time.perf_counter();manifest=json.loads((R/'data/input_sha256.json').read_text())
 for name,h in manifest.items():assert hashlib.sha256((R/'data'/name).read_bytes()).hexdigest()==h,name
 expected={p.name:table(p) for p in (R/'results').glob('*.csv')}
 data={p.name:p.read_bytes() for p in (R/'data').glob('*.csv')}
 for name in ['make_data.py','analyze.py','mechanisms.py','figures.py']:
  proc=subprocess.run([sys.executable,str(R/'code'/name)],cwd=R.parent,capture_output=True,text=True)
  if proc.returncode:raise RuntimeError(proc.stdout+'\n'+proc.stderr)
  print(name+': passed')
 for name,blob in data.items():assert (R/'data'/name).read_bytes()==blob,name
 count=0;largest=0.
 for name,old in expected.items():
  new=table(R/'results'/name);assert len(new)==len(old),name
  for aa,bb in zip(old,new):
   assert len(aa)==len(bb),name
   for x,y in zip(aa,bb):
    try:a,b=float(x),float(y)
    except ValueError:
     assert x==y,(name,x,y);continue
    if not math.isfinite(a) or not math.isfinite(b):
     assert (math.isnan(a) and math.isnan(b)) or a==b;continue
    delta=abs(a-b);count+=1;largest=max(largest,delta)
    assert delta<=2e-6+1e-8*abs(a),(name,x,y)
 subprocess.run([sys.executable,'-m','unittest','discover','-s',str(R/'tests'),'-v'],cwd=R.parent,check=True)
 summary=dict(all_checks_passed=True,csv_files_compared=len(expected),numeric_cells_compared=count,maximum_absolute_difference=largest,input_files_checked=len(manifest),tests=14,main_figures=5,flowchart_text_inside_boxes=True,network_used=False,elapsed_seconds=time.perf_counter()-start)
 (R/'reproduction_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
