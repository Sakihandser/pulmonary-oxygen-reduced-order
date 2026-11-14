from pathlib import Path
import sys,subprocess,csv,json,math,hashlib
R=Path(__file__).resolve().parents[1]
def read(p):
 with p.open(encoding='utf-8',newline='') as f:return list(csv.reader(f))
expected={p.name:read(p) for p in (R/'results').glob('*.csv')}
for script in ['analyze_review.py','figures_review.py']:
 subprocess.run([sys.executable,str(R/'code'/script)],cwd=R.parent,check=True,stdout=subprocess.DEVNULL)
count=0;maxerr=0.
for name,old in expected.items():
 new=read(R/'results'/name);assert len(old)==len(new),name
 for a,b in zip(old,new):
  assert len(a)==len(b),name
  for x,y in zip(a,b):
   try:xx,yy=float(x),float(y)
   except ValueError:assert x==y,(name,x,y);continue
   if math.isnan(xx) or math.isnan(yy):assert math.isnan(xx) and math.isnan(yy);continue
   if math.isinf(xx) or math.isinf(yy):assert xx==yy;continue
   err=abs(xx-yy);maxerr=max(maxerr,err);count+=1
   assert err<=2e-7+1e-9*abs(xx),(name,x,y)
subprocess.run([sys.executable,'-m','unittest','discover','-s',str(R/'tests'),'-v'],cwd=R.parent,check=True)
summary=dict(result_csv_files=len(expected),numeric_cells_compared=count,maximum_absolute_difference=maxerr,tests=9,network_used=False,interpretation='Deterministic computational verification, not prospective physiological validation')
(R/'reproduction_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
