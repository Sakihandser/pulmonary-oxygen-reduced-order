from pathlib import Path
import csv,json,math,subprocess,sys,time
R=Path(__file__).resolve().parents[1]
def rows(p):
    with p.open(encoding='utf-8') as f:return list(csv.DictReader(f))
def main():
    start=time.perf_counter();expected={p.name:rows(p) for p in (R/'results').glob('*.csv')}
    subprocess.run([sys.executable,str(R/'code/analyze_applicability.py')],check=True,stdout=subprocess.DEVNULL)
    count=0;maximum=0.
    for n,old in expected.items():
        new=rows(R/'results'/n);assert len(old)==len(new),n
        for x,y in zip(old,new):
            assert x.keys()==y.keys(),n
            for key in x:
                try:v=float(x[key]);w=float(y[key])
                except ValueError:assert x[key]==y[key];continue
                if math.isnan(v):assert math.isnan(w);continue
                assert math.isfinite(v) and math.isfinite(w)
                count+=1;err=abs(v-w);maximum=max(maximum,err)
                assert err<=2e-9+1e-10*abs(v),(n,key,v,w)
    subprocess.run([sys.executable,'-m','unittest','discover','-s',str(R/'tests'),'-v'],check=True)
    subprocess.run([sys.executable,str(R/'code/figures_applicability.py')],check=True)
    report=dict(csv_files=len(expected),numeric_cells=count,maximum_absolute_difference=maximum,tests=10,all_checks_passed=True,elapsed_seconds=time.perf_counter()-start,interpretation='Numerical verification only')
    (R/'reproduction_summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
