from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];F=R/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','savefig.facecolor':'white'})
with (R/'results/flow_response.csv').open() as f:rows=list(csv.DictReader(f))
fig,axs=plt.subplots(1,3,figsize=(10,3.2),layout='constrained')
colors=['#007F7A','#407AA1','#7669AD','#DA8A2E','#C95454']
for beta,color in zip([0,.25,.5,.75,1],colors):
 z=[r for r in rows if float(r['beta'])==beta];x=[float(r['Q_ratio']) for r in z]
 for ax,key in zip(axs[:2],['Ca_ratio','DO2_ratio']):ax.plot(x,[float(r[key]) for r in z],color=color,label=f'{beta:g}',lw=1.8)
for ax in axs[:2]:ax.axhline(1,color='#777777',ls=':',lw=1);ax.axvline(1,color='#bbbbbb',ls=':',lw=.7);ax.set(xlabel='Flow / original flow',xlim=(.6,2))
axs[0].set(title='A  Arterial content',ylabel='Content / original content');axs[1].set(title='B  Oxygen delivery',ylabel='Delivery / original delivery')
axs[1].legend(title='Scaling exponent β',loc='upper left',fontsize=7,title_fontsize=8,frameon=False)
s=np.linspace(.01,.8,160);crit=(1-s)/(-np.log(s));ax=axs[2];ax.fill_between(s,0,crit,color='#DDF0EB');ax.fill_between(s,crit,1,color='#F2E3D5');ax.plot(s,crit,color='#233E55')
ax.scatter([.2],[(1-.2)/-np.log(.2)],s=18,color='black');ax.annotate('s = 0.20\nβcrit = 0.497',xy=(.2,(1-.2)/-np.log(.2)),xytext=(.32,.36),fontsize=8,arrowprops={'arrowstyle':'-','lw':.7})
ax.text(.16,.12,'Content rises',fontsize=8);ax.text(.12,.91,'Content falls',fontsize=8)
ax.set(title='C  Local response to flow',xlabel='Current effective fraction s',ylabel='Scaling exponent β',xlim=(.01,.8),ylim=(0,1))
fig.savefig(F/'Fig6.png',dpi=180,bbox_inches='tight');fig.savefig(F/'Fig6.tif',dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'});fig.savefig(F/'Fig6.svg',bbox_inches='tight');plt.close(fig)
print('Exported final Figure 6')
