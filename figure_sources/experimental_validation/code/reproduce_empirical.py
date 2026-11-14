"""Reproduce the exploratory empirical extension offline in a disposable copy."""
from pathlib import Path
import csv,hashlib,json,math,platform,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]

def table(path):
    with path.open(encoding='utf-8',newline='') as f:return list(csv.reader(f))

def main():
    start=time.perf_counter()
    manifest=json.loads((ROOT/'input_sha256.json').read_text())
    for name,digest in manifest.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,('Input hash',name)
    expected={p.name:table(p) for p in (ROOT/'results').glob('*.csv')}
    data={p.name:p.read_bytes() for p in (ROOT/'data').glob('*observations.csv')}
    for name in ['make_data.py','validate.py','couple_network.py','analyze_primary.py','review_extension.py','figures.py','review_figures.py']:
        result=subprocess.run([sys.executable,str(ROOT/'code'/name)],cwd=ROOT.parent,text=True,capture_output=True)
        if result.returncode:raise RuntimeError(result.stdout+'\n'+result.stderr)
        print(name+': passed')
    for name,original in data.items():assert (ROOT/'data'/name).read_bytes()==original,('Transcription changed',name)
    cells=0;largest=0.
    for name,old in expected.items():
        new=table(ROOT/'results'/name)
        assert len(old)==len(new),name
        for a,b in zip(old,new):
            assert len(a)==len(b),name
            for x,y in zip(a,b):
                try:xx,yy=float(x),float(y)
                except ValueError:
                    assert x==y,(name,x,y);continue
                if not math.isfinite(xx) or not math.isfinite(yy):
                    assert (math.isnan(xx) and math.isnan(yy)) or xx==yy;continue
                delta=abs(xx-yy);largest=max(largest,delta);cells+=1
                assert delta<=2e-6+1e-8*abs(xx),(name,xx,yy)
    subprocess.run([sys.executable,'-m','unittest','discover','-s',str(ROOT/'tests'),'-v'],cwd=ROOT.parent,check=True)
    summary=dict(all_checks_passed=True,input_files_checked=len(manifest),csv_files_compared=len(expected),
                 numeric_cells_compared=cells,maximum_absolute_difference=largest,python=platform.python_version(),
                 figures=['Fig8','Fig9','Fig10','Fig11','Fig12'],elapsed_seconds=time.perf_counter()-start,network_used=False,
                 interpretation='Exploratory conditional endpoint comparison; not clinical validation.')
    (ROOT/'reproduction_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
