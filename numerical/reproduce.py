"""Reproduce the article's deterministic numerical calculations offline."""
from pathlib import Path
import csv,hashlib,json,math,platform,shutil,subprocess,sys,time
ROOT=Path(__file__).resolve().parent
PHASE=ROOT/'sensitivity'
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def table(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.reader(f))
def execute(p):subprocess.run([sys.executable,str(ROOT/p)],check=True,cwd=ROOT)
def comparable_cells(left,right):
    """Compare scalar/vector numeric cells; preserve exact nonnumeric text."""
    try:
        x=[float(v) for v in left.split(';')]
        y=[float(v) for v in right.split(';')]
    except ValueError:
        assert left==right,('text mismatch',left,right)
        return []
    assert len(x)==len(y),('vector length',left,right)
    return list(zip(x,y))
def reproduce():
    start=time.perf_counter()
    manifest=json.loads((ROOT/'release_sha256.json').read_text(encoding='utf-8'))
    for name,digest in manifest.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,('File integrity',name)
    expected={p.relative_to(ROOT).as_posix():read(p) for p in ROOT.rglob('*.csv') if p.parent.name=='results'}
    for code in ['aggregate_analysis/code/analyze.py','observation_analysis/code/analyze_observations.py',
                 'identifiability/code/joint_identifiability.py','controlled_exchange/code/reproduce_diagnostic.py',
                 'sensitivity/code/run_experiments.py','common_calibration/code/run.py','common_calibration/code/mechanism.py','robustness/code/run.py']:execute(code)
    count=0;maximum=0.
    for name,old in expected.items():
        new=read(ROOT/name);assert len(old)==len(new),name
        for a,b in zip(old,new):
            assert a.keys()==b.keys(),name
            for key in a:
                for x,y in comparable_cells(a[key],b[key]):
                    if math.isnan(x) or math.isnan(y):assert math.isnan(x) and math.isnan(y);continue
                    if math.isinf(x) or math.isinf(y):assert x==y;continue
                    err=abs(x-y);maximum=max(maximum,err);count+=1
                    assert err<=2e-6+1e-8*abs(x),(name,key,x,y)
    source=read(ROOT/'aggregate_analysis/data/source_cells.csv')
    assert len(source)==102 and source==read(ROOT/'comparison_data/results/source_cells_reextracted.csv')
    report=PHASE/'reporting';result=PHASE/'results'
    krows=read(result/'k_scan.csv');fb=read(result/'feedback_grid.csv');total=read(result/'total_demand_closure.csv')
    round4=lambda x:f'{0.0 if round(float(x),4)==0 else float(x):.4f}'
    def ca(k,mode,alloc):return float(next(r['ca'] for r in krows if r['k']==str(k) and r['mode']==mode and r['allocation']==alloc))
    bodies={4:[],6:[],7:[],8:[]}
    for label,mode in zip(['Saturation relaxation','Content relaxation','Endpoint matched D','Gas mapped D'],['saturation','content','matched','gas']):
        v=lambda k,a:ca(k,mode,a)
        bodies[4].append([label,round4(v(1,'fixed')),round4(v(25,'fixed')),round4(v(1,'fixed')-v(25,'fixed')),round4(v(25,'flow')-v(25,'fixed'))])
    path=read(ROOT/'controlled_exchange/results/sequential_decomposition.csv');printed=table(report/'main_table7.csv')
    for i,r in enumerate(path):
        bodies[6].append([printed[i+1][0],round4(r['ca_ml_dl']),round4(r['redistribution_content_loss']),round4(r['redistribution_sa_loss_pp']),printed[i+1][-1]])
    for mode in ['matched','gas']:
        for closure in ['capped','total']:
            r=next(r for r in (fb if closure=='capped' else total) if r['mode']==mode and r['k']=='25' and r['M']=='250' and r['scale']=='1')
            bodies[7].append([mode,closure,*[f'{float(r[key]):.3f}' for key in ['ca','cv','pv']],*[f'{float(r[key]):.2f}' for key in ['U','unmet']]])
    for mode,label in [('matched','Endpoint matched'),('gas','Gas mapped')]:
        for k in [1,25]:
            bodies[8].append([label,str(k),round4(ca(k,mode,'fixed')),round4(ca(k,mode,'flow')),round4(ca(k,mode,'flow')-ca(k,mode,'fixed'))])
    for n,body in bodies.items():assert body==table(report/f'main_table{ {4:5,6:7,7:8,8:6}[n]}.csv')[1:],('Main table',n)
    display=read(report/'main_table4.csv')
    transfers=read(ROOT/'common_calibration/results/transfer.csv')
    for row in display:
        for law,col in [('saturation','Saturation'),('content','Content'),('constant_D','Constant D')]:
            raw=next(x for x in transfers if x['match']=='endpoint' and x['law']==law and float(x['Pv'])==float(row['Pv']) and abs(float(x['Q'])-float(row['Q']))<5.1e-6 and float(x['PA'])==float(row['PA']))
            assert abs(float(raw['content'])-float(row[col]))<=5.1e-6
    # Refresh figure sources with the results calculated in this invocation.
    for name in ['k_scan.csv','D_scan.csv','demand_path.csv','hydraulic_scan.csv','external_alvin_benchmark.csv','ratio_scan.csv']:
        shutil.copy2(result/name,report/name)
    for name in ['ablation_matrix.csv','equivalent_conductance.csv','sequential_decomposition.csv']:
        shutil.copy2(ROOT/'controlled_exchange/results'/name,report/name)
    for code in ['plot_main_figures.py','plot_sensitivity_figures.py']:execute(Path('sensitivity/reporting/tools')/code)
    shutil.copy2(ROOT/'common_calibration/results/transfer.csv',report/'common_transfer.csv')
    shutil.copy2(ROOT/'common_calibration/results/allocation_scan.csv',report/'allocation_scan.csv')
    for name in ['dissolved_grid.csv','dissolved_summary.csv','dissolved_calibration.csv','signed_state_map.csv','ranked_grid_states.csv']:
        shutil.copy2(ROOT/'common_calibration/results'/name,report/name)
    execute('sensitivity/reporting/tools/plot_calibration_maps.py')
    assert len(list((report/'figures').glob('*.tif')))==9
    import numpy,scipy,matplotlib
    summary=dict(numerical_csv_files_compared=len(expected),numeric_cells_compared=count,
                 maximum_absolute_csv_difference=maximum,transcribed_source_cells_checked=len(source),
                 numerical_main_tables_verified=[4,5,6,7,8],frozen_transfer_rows_checked=len(display),
                 main_tables_supplied=10,figures_exported=9,all_checks_passed=True,
                 mechanism_checks=json.loads((ROOT/'common_calibration/results/checks.json').read_text()),
                 robustness_checks=json.loads((ROOT/'robustness/results/summary.json').read_text()),
                 native_kernel_check='comparison with supplied compiled C# nodes; native recompilation is a separate optional command',
                 python=platform.python_version(),numpy=numpy.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,
                 elapsed_seconds=time.perf_counter()-start,network_used=False,physiological_validation=False)
    (ROOT/'reproduction_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':reproduce()
