"""Generate state-selection maps, allocation surface, and dissolved-oxygen control."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from PIL import Image
ROOT=Path(__file__).resolve().parents[1];FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
def read(n):return list(csv.DictReader((ROOT/n).open(encoding='utf-8-sig')))
def save(fig,n):
    fig.savefig(FIG/f'Fig{n}.png',dpi=300);fig.savefig(FIG/f'Fig{n}.svg')
    import io
    buf=io.BytesIO();fig.savefig(buf,format='png',dpi=600);buf.seek(0)
    with Image.open(buf) as im:im.convert('RGB').save(FIG/f'Fig{n}.tif',compression='tiff_lzw',dpi=(600,600))
    plt.close(fig)
plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'path'})
rows=read('common_transfer.csv');maps=read('signed_state_map.csv')
fig,axs=plt.subplots(2,2,figsize=(7.4,6.2),layout='constrained')
for ax,match,title in zip(axs[0],['endpoint','inlet_slope'],['A  Common endpoint match','B  Common inlet slope match']):
    for law,label,color in [('content','Content relaxation','#315b78'),('constant_D','Constant D','#a9633c')]:
        rr=[r for r in rows if r['match']==match and r['law']==law]
        ax.plot(range(7),[float(x['residual']) for x in rr],'o-',label=label,color=color,ms=4)
    ax.axhline(0,color='gray',lw=.7);ax.set(xticks=range(7),xticklabels=['Ref','Pv30','Pv50','Q3','Q5','PA80','PA120'],title=title,ylabel='Difference from saturation (mL/dL)');ax.tick_params(axis='x',labelrotation=45,labelsize=8)
axs[0,0].legend(fontsize=8,frameon=False)
styles=[':','--','-'];levels=[.0200002,.1000002,.2000002]
for ax,match,title in zip(axs[1],['endpoint','inlet_slope'],['C  Endpoint match state map','D  Inlet slope match state map']):
    rr=[r for r in maps if r['match']==match];pv=sorted({float(r['Pv']) for r in rr});qs=sorted({float(r['Q']) for r in rr})
    z=np.array([[float(next(r['signed_difference'] for r in rr if float(r['Pv'])==p and float(r['Q'])==q)) for p in pv] for q in qs])
    im=ax.pcolormesh(pv,qs,z,shading='nearest',cmap='RdBu_r',norm=TwoSlopeNorm(0,vmin=-1.5,vmax=1.5))
    ax.contour(pv,qs,abs(z),levels=levels,colors='black',linestyles=styles,linewidths=.9)
    ax.plot(40,3.9464411557434818,'*',mfc='#f4cf37',mec='black',ms=9)
    ax.set(xlabel='Inlet oxygen pressure (mmHg)',ylabel='Flow (L/min)',title=title,xlim=(30,50),ylim=(3,5))
fig.colorbar(im,ax=list(axs[1]),label='Constant D minus saturation (mL/dL)',shrink=.8)
fig.legend(handles=[Line2D([],[],c='black',ls=s,label=f'e = {e:.2f} mL/dL') for s,e in zip(styles,[.01,.05,.10])],loc='outside lower center',ncol=3,frameon=False,fontsize=8)
save(fig,2)
scan=read('allocation_scan.csv');fig,axs=plt.subplots(1,2,figsize=(7.4,3.2),layout='constrained')
for ax,lam,title in zip(axs,[0,.5],['A  Fixed CT shares','B  Halfway allocation interpolation']):
    rr=[r for r in scan if float(r['lambda_value'])==lam];xs=sorted({float(r['LUL_fraction']) for r in rr});ds=sorted({float(r['D']) for r in rr})
    z=np.array([[float(next(r['redistribution_loss'] for r in rr if float(r['D'])==d and float(r['LUL_fraction'])==x)) for x in xs] for d in ds])
    im=ax.pcolormesh(xs,ds,z,shading='nearest',cmap='viridis',vmin=0,vmax=2.1)
    ax.axvline(xs[1],color='white',lw=.7,ls=':');ax.axvline(xs[2],color='white',lw=.7,ls=':')
    ax.set(xlabel='LUL flow fraction',ylabel='Total D (mL/min/mmHg)',title=title)
fig.colorbar(im,ax=axs,label='k=1 minus k content (mL/dL)',shrink=.9);save(fig,5)
grid=read('dissolved_grid.csv');fig,axs=plt.subplots(1,2,figsize=(7.4,3.0),layout='constrained')
for ax,match,title in zip(axs,['endpoint','inlet_slope'],['A  Endpoint match','B  Inlet slope match']):
    for b,ls,col in [(0.,'--','#444444'),(.0031,'-','#315b78')]:
        rr=[r for r in grid if r['match']==match and float(r['dissolved_coefficient'])==b and float(r['Q'])==5]
        ax.plot([float(r['Pv']) for r in rr],[float(r['signed_difference']) for r in rr],ls,c=col,label=f'b = {b:g}')
    ax.set(xlabel='Inlet oxygen pressure (mmHg)',ylabel='Content minus saturation (mL/dL)',title=title);ax.axhline(0,lw=.5,c='gray');ax.legend(frameon=False,fontsize=8)
save(fig,'S2')
