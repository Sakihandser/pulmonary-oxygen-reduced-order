"""Regenerate Fig5 and FigS1 from saved full-precision CSV files."""
from pathlib import Path
import csv
from io import BytesIO
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
PHASE=Path(__file__).resolve().parents[1]
RESULT=PHASE;FIG=PHASE/'figures';FIG.mkdir(exist_ok=True)
def read(name):return list(csv.DictReader((RESULT/name).open(encoding='utf-8')))
def save(fig,stem):
    fig.savefig(FIG/(stem+'.png'),dpi=220);fig.savefig(FIG/(stem+'.svg'))
    buf=BytesIO();fig.savefig(buf,format='png',dpi=600);buf.seek(0)
    with Image.open(buf) as source:rgb=source.convert('RGB')
    with (FIG/(stem+'.tif')).open('wb') as f:rgb.save(f,format='TIFF',compression='tiff_lzw',dpi=(600,600))
    plt.close(fig)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
ks=read('k_scan.csv');ds=read('D_scan.csv');dense=read('demand_path.csv')
hyd={r['k']:float(r['f_LUL']) for r in read('hydraulic_scan.csv')}
fig=plt.figure(figsize=(7.4,6.1),layout='constrained');fig.suptitle('Synthetic stress tests: declared parameter slices',fontsize=11);grid=fig.add_gridspec(2,2)
ax=[fig.add_subplot(grid[0,0]),fig.add_subplot(grid[0,1]),fig.add_subplot(grid[1,:])]
colors=['#255b83','#d58e32','#a54545','#327e65']
for mode,col,label in zip(['saturation','content','matched','gas'],colors,['Saturation','Content','Matched D','Gas D']):
    for alloc,ls in [('fixed','-'),('flow','--')]:
        rr=[r for r in ks if r['mode']==mode and r['allocation']==alloc]
        baseline=float(rr[0]['ca'])
        ax[0].plot([hyd[r['k']] for r in rr],[baseline-float(r['ca']) for r in rr],ls,marker='o' if alloc=='fixed' else None,c=col,label=label if alloc=='fixed' else None)
ax[0].set(xlabel='LUL flow fraction',ylabel='Loss from k=1 (mL/dL)',title='A  Flow concentration');ax[0].legend(fontsize=8)
for alloc,ls in [('fixed','-'),('flow','--')]:
    dvals=list(dict.fromkeys(r['D'] for r in ds));loss=[]
    for val in dvals:
        pair=[next(float(r['ca']) for r in ds if r['k']==k and r['D']==val and r['allocation']==alloc) for k in ['1','25']];loss.append(pair[0]-pair[1])
    ax[1].plot(list(map(float,dvals)),loss,ls,marker='o',label=alloc)
ax[1].set(xlabel='D (mL/min/mmHg)',ylabel='k=1 minus k=25 (mL/dL)',title='B  Same law different capacity');ax[1].legend(fontsize=8)
for mode,col in [('matched',colors[2]),('gas',colors[3])]:
    rr=[r for r in dense if r['mode']==mode]
    ax[2].plot([float(r['M']) for r in rr],[float(r['U']) for r in rr],c=col,label=mode+' actual use')
ax[2].plot([50,400],[50,400],':',c='gray',label='specified demand')
ax[2].set(xlabel='Demand (mL/min)',ylabel='Actual use (mL/min)',title='C  Capped feedback at k=25');ax[2].legend(fontsize=8)
save(fig,'Fig6')
bench=read('external_alvin_benchmark.csv');rat=read('ratio_scan.csv')
fig,ax=plt.subplots(1,2,figsize=(7.4,3.4),layout='constrained')
for mode,col in zip(['saturation','content','matched','gas'],colors):
    rr=[r for r in bench if r['mode']==mode and r['reference_fit']=='False']
    ax[0].plot(range(4),[float(r['residual']) for r in rr],'o-',label=mode,c=col)
ax[0].set(xticks=range(4),xticklabels=['Pv30','Pv50','Q3','Q5'],ylabel='Residual (mL/dL)',title='A  External pressure kernel');ax[0].legend(fontsize=8)
for mode,col in [('saturation',colors[0]),('content',colors[1])]:
    rr=[r for r in rat if r['mode']==mode and r['k']=='25' and r['allocation']=='fixed']
    ax[1].plot([float(r['ratio']) for r in rr],[float(r['ca']) for r in rr],'o-',c=col,label=mode)
ax[1].set(xlabel='Vc/tau relative to reference',ylabel='Content (mL/dL)',title='B  Relaxation capacity');ax[1].legend(fontsize=8)
save(fig,'FigS1')
