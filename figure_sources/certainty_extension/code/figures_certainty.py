from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parents[1]
def read(n):
 with (HERE/'results'/n).open(encoding='utf-8') as f:return list(csv.DictReader(f))
def save(fig,n):
 (HERE/'figures').mkdir(exist_ok=True)
 for ext in ['png','tif','svg']:fig.savefig(HERE/'figures'/f'{n}.{ext}',dpi=600,bbox_inches='tight',facecolor='white')
 plt.close(fig)
plt.rcParams.update({'font.size':9,'font.family':'DejaVu Sans','svg.fonttype':'none'})
rows=read('decision_bounds.csv');fig,ax=plt.subplots(1,2,figsize=(9.4,4.1),sharex=True,sharey=True)
names=['baseline','baseline_capacity','one_extra','two_extra','two_extra_noise005','two_extra_noise010']
labels=['Baseline only','Baseline + total capacity','Add r = 1.25','Add r = 1.25 and 1.5','Two extra states, error ±0.005','Two extra states, error ±0.010']
for a,model in zip(ax,['A','B']):
 for j,name in enumerate(names):
  z=next(x for x in rows if x['model']==model and x['information']==name);lo=float(z['Ca_lower']);hi=float(z['Ca_upper'])
  col={'all_recover':'#12805c','none_recover':'#a04239','opposite_witnesses':'#365f91','unresolved':'#777777'}[z['decision']]
  a.plot([lo,hi],[j,j],lw=5,color=col,solid_capstyle='butt');a.plot([(lo+hi)/2],[j],'o',color=col,ms=4)
 a.axvline(float(rows[0]['target']),ls='--',c='black',lw=1,label='Recovery target')
 a.set_title(f'{model}  Shared observations from structure {model}')
 a.set_yticks(range(len(names)),labels);a.set_xlabel('Outer content bounds (mL/dL)');a.grid(axis='x',alpha=.2)
ax[0].invert_yaxis();fig.tight_layout();save(fig,'Fig10')
alias=read('alias_curves.csv');cap=read('capacity_exclusion.csv')
fig,ax=plt.subplots(1,2,figsize=(9.1,3.7))
rr=[z for z in alias if float(z['exposure_cap'])==16]
ax[0].plot([float(z['r']) for z in rr],[float(z['fractional_s']) for z in rr],label='Homogeneous, β = 1/2',color='#365f91',lw=2)
ax[0].plot([float(z['r']) for z in rr],[float(z['fixed_mixture_s']) for z in rr],label='Finite mixture, β = 1',color='#bb6723',ls='--',lw=1.5)
ax[0].set(xlabel='Flow ratio r',ylabel='Residual exchange fraction',title='A  Different mechanisms share a response')
ax[0].legend(fontsize=8);ax[0].grid(alpha=.2)
cc=[z for z in cap if z['cells']=='8192']
xx=np.arange(len(cc))
for i,z in enumerate(cc):
 ax[1].plot([i,i],[float(z['K_lower']),float(z['K_upper'])],lw=5,c='#365f91')
ax[1].axhline(10,ls='--',color='black',lw=1,label='Illustrative independent upper bound 10')
ax[1].set_xticks(xx,[z['residual_halfwidth'] for z in cc]);ax[1].set(xlabel='Allowed residual halfwidth',ylabel='Minimum total capacity (L/min)',title='B  Independent capacity can exclude a class')
ax[1].legend(fontsize=7,loc='upper right');ax[1].grid(axis='y',alpha=.2)
fig.tight_layout();save(fig,'Fig11')
