from pathlib import Path
import subprocess,sys,csv,json,math
H=Path(__file__).resolve().parents[1]
def read(p):
 with p.open(encoding='utf8') as f:return list(csv.DictReader(f))
old={p.name:read(p) for p in (H/'results').glob('*.csv')}
subprocess.run([sys.executable,str(H/'code/analyze_measurement.py')],check=True)
count=0;maximum=0
for name,rows in old.items():
 new=read(H/'results'/name);assert len(rows)==len(new)
 for a,b in zip(rows,new):
  assert a.keys()==b.keys()
  for k in a:
   try:x,y=float(a[k]),float(b[k])
   except ValueError:assert a[k]==b[k];continue
   count+=1;maximum=max(maximum,abs(x-y));assert math.isclose(x,y,rel_tol=1e-8,abs_tol=1e-7)
subprocess.run([sys.executable,'-m','unittest','discover','-s',str(H/'tests'),'-v'],check=True)
v=dict(csv_files=len(old),numeric_cells=count,max_absolute_difference=maximum,tests=6,all_checks_passed=True)
(H/'reproduction_summary.json').write_text(json.dumps(v,indent=2),encoding='utf8')
print(json.dumps(v))
