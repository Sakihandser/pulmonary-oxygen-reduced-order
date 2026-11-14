from pathlib import Path
import sys,subprocess,csv,json,math
R=Path(__file__).resolve().parents[1]
def read(p):
 with p.open(encoding='utf-8',newline='') as f:return list(csv.reader(f))
expected={p.name:read(p) for p in (R/'results').glob('*.csv')}
for script in ['analyze_coupling.py','figures_coupling.py']:subprocess.run([sys.executable,str(R/'code'/script)],check=True,stdout=subprocess.DEVNULL)
n=0;maximum=0.
for name,old in expected.items():
 new=read(R/'results'/name);assert len(old)==len(new),name
 for aa,bb in zip(old,new):
  assert len(aa)==len(bb),name
  for x,y in zip(aa,bb):
   try:a,b=float(x),float(y)
   except ValueError:assert x==y,(name,x,y);continue
   if math.isnan(a) or math.isnan(b):assert math.isnan(a) and math.isnan(b);continue
   if math.isinf(a) or math.isinf(b):assert a==b;continue
   err=abs(a-b);assert err<=2e-7+1e-9*abs(a),(name,x,y);n+=1;maximum=max(maximum,err)
subprocess.run([sys.executable,'-m','unittest','discover','-s',str(R/'tests'),'-v'],check=True)
summary=dict(result_csv_files=len(expected),numeric_cells_compared=n,maximum_absolute_difference=maximum,tests=8,network_used=False)
(R/'reproduction_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
