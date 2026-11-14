from pathlib import Path
import sys,subprocess,csv,json,math,time
H=Path(__file__).resolve().parents[1]
def read(p):
 with p.open(encoding='utf-8') as f:return list(csv.DictReader(f))
old={p.name:read(p) for p in (H/'results').glob('*.csv')};start=time.perf_counter()
for f in ['analyze_certainty.py','figures_certainty.py']:subprocess.run([sys.executable,str(H/'code'/f)],check=True)
count=0;maximum=0.
for name,rows in old.items():
 new=read(H/'results'/name);assert len(rows)==len(new),name
 for a,b in zip(rows,new):
  assert a.keys()==b.keys()
  for key in a:
   try:x=float(a[key]);y=float(b[key])
   except ValueError:assert a[key]==b[key];continue
   if not math.isfinite(x):assert a[key]==b[key];continue
   count+=1;maximum=max(maximum,abs(x-y));assert abs(x-y)<=1e-7+1e-8*abs(x),(name,key,x,y)
subprocess.run([sys.executable,'-m','unittest','discover','-s',str(H/'tests'),'-v'],check=True)
s=dict(csv_files=len(old),numeric_cells=count,maximum_absolute_difference=maximum,tests=13,all_checks_passed=True,elapsed_seconds=time.perf_counter()-start,interpretation='Numerical verification of specified classes, not clinical validation')
(H/'reproduction_summary.json').write_text(json.dumps(s,indent=2),encoding='utf-8');print(json.dumps(s,indent=2))
