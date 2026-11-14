from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1]
with (R/'results/heterogeneous_response.csv').open(encoding='utf-8') as f:rows=list(csv.DictReader(f))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'baseline_heterogeneity'})
fig,axs=plt.subplots(1,2,figsize=(9.2,3.6))
styles={'finite_A':('#0072B2','-','Finite model A'),'matched_B':('#D55E00','-','Matched model B'),'homogeneous':('#555555','--','Homogeneous bound'),'ideal_limit':('#999999',':','Ideal upper bound')}
for name,(color,ls,label) in styles.items():
    rr=[r for r in rows if r['model']==name];q=[float(r['Q_ratio']) for r in rr]
    for ax,key in zip(axs,['Ca_ratio','DO2_ratio']):ax.plot(q,[float(r[key]) for r in rr],color=color,ls=ls,label=label,lw=1.8)
for ax,title,ylabel in zip(axs,['A  Arterial content','B  Oxygen delivery'],['Content / original target','Delivery / original target']):
    ax.axhline(1,color='black',lw=.8,alpha=.65);ax.set(xlabel='Flow / reference flow',ylabel=ylabel,title=title,xlim=(1,2));ax.grid(alpha=.16)
axs[0].legend(frameon=False,fontsize=8,loc='lower left')
axs[0].text(.02,.97,'A and B: same baseline\nand total fixed capacity',transform=axs[0].transAxes,ha='left',va='top',fontsize=8)
fig.tight_layout(pad=1.1,w_pad=1.7)
out=R/'figures';out.mkdir(exist_ok=True)
for ext in ['png','tif','svg']:
    kwargs={'metadata':{'Date':None}} if ext=='svg' else {}
    fig.savefig(out/f'Fig7.{ext}',dpi=300,**kwargs)
plt.close(fig)
