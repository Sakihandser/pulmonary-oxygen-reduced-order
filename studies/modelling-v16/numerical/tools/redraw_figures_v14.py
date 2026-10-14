from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures_regenerated'
(OUT/'Figures_v14').mkdir(parents=True,exist_ok=True)
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
# Regenerate submission figures from retained data at native 300 dpi.
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'path'})
DATA=ROOT/'10_diagnostic_revision/results';FIG=OUT/'Figures_v14'
def savefig(fig,n):
 from PIL import Image
 fig.savefig(FIG/f'Fig{n}.png',dpi=300)
 fig.savefig(FIG/f'Fig{n}.svg')
 with Image.open(FIG/f'Fig{n}.png') as im:im.convert('RGB').save(FIG/f'Fig{n}.tif',compression='tiff_lzw',dpi=(300,300))
 plt.close(fig)
fig,ax=plt.subplots(figsize=(7.4,3.2));fig.subplots_adjust(left=.015,right=.985,top=.98,bottom=.03);ax.axis('off');ax.set(xlim=(0,1),ylim=(0,1))
for x,y,w,h,t in [(.015,.52,.29,.39,'Archived five-lobe flows\nk = 1 and 25\nClosure fixes PVR and Q'),(.355,.52,.29,.39,'Four exchange alternatives\nShared blood relation\nExplicit capacity calibration'),(.695,.52,.29,.39,'Flow-weighted content\nCm, Ca and saturation\nUptake and oxygen delivery'),(.355,.04,.29,.31,'Boundary factors\nHemoglobin, inlet feedback\nPost-exchange admixture'),(.695,.04,.29,.31,'Independent observation test\nNot available in this study\nSame-state data required')]:
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.005',facecolor='#edf2f6',edgecolor='#315b78'))
 ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=9)
for a,b in [((.305,.715),(.355,.715)),((.645,.715),(.695,.715)),((.5,.35),(.5,.52)),((.84,.52),(.84,.35))]:ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'->','color':'#315b78'})
ax.text(.018,.2,'Hydraulic feedback to recruitment\nand distension remains a separate\nunvalidated research prototype.',fontsize=8.5,va='center');savefig(fig,1)
common=read(DATA/'common_boundary_comparison.csv');x=np.arange(4)
fig,axs=plt.subplots(1,2,figsize=(7.4,3.4),layout='constrained')
labels=['Saturation\nrelaxation','Content\nrelaxation','Endpoint\nmatched D','Gas\nmapped D']
for k,dx,marker in [(1,-.07,'o'),(25,.07,'s')]:axs[0].scatter(x+dx,[float(r[f'k{k}_mixed_end_ml_dl']) for r in common],marker=marker,label=f'k = {k}')
axs[0].set(ylabel='Mixed terminal content (mL/dL)',title='A  Exchange alternatives');axs[0].legend(frameon=False)
vals=[float(r['loss_ml_dl']) for r in common];axs[1].bar(x,vals,color=['#315b78','#548db3','#a9633c','#46856e'])
for i,v in enumerate(vals):axs[1].text(i,v+.04,f'{v:.6f}',ha='center',fontsize=8)
axs[1].set(ylabel='k = 1 minus k = 25 (mL/dL)',ylim=(0,2.12),title='B  Common boundary contrast')
for a in axs:a.set_xticks(x,labels,fontsize=8)
savefig(fig,2)
prof=read(DATA/'equivalent_conductance.csv');fig,ax=plt.subplots(figsize=(7.4,2.6),layout='constrained')
for key,label,ls in [('saturation_equivalent_d','Saturation relaxation','-'),('content_equivalent_d','Content relaxation','-'),('endpoint_matched_d','Endpoint matched D','--'),('gas_mapped_d','Gas mapped D',':')]:ax.plot([float(r['po2']) for r in prof],[float(r[key]) for r in prof],ls,label=label)
ax.set(xlabel='Capillary PO2 (mmHg)',ylabel='Equivalent D (mL/min/mmHg)');ax.legend(frameon=False,fontsize=9);savefig(fig,3)
rows=read(DATA/'ablation_matrix.csv');path=read(DATA/'sequential_decomposition.csv')
fig,axs=plt.subplots(1,3,figsize=(7.5,3.4),layout='constrained')
axs[0].plot(range(7),[float(r['redistribution_content_loss']) for r in path],'o-',color='#315b78')
axs[0].set(xticks=range(7),xticklabels=['Legacy','Content','Match D','Gas D','Hb','Inlet','Mix'],ylabel='Content contrast (mL/dL)',title='A  Table 2 path');axs[0].tick_params(axis='x',labelsize=8,rotation=60)
for ax,mode,col,title in [(axs[1],'conductance_matched','#a9633c','B  Endpoint matched D'),(axs[2],'conductance_conditional','#46856e','C  Gas mapped D')]:
 losses=[]
 for fb,mix in [(0,0),(0,1),(1,0),(1,1)]:
  pair=[next(r for r in rows if r['k']==str(k) and r['mode']==mode and r['allocation']=='fixed' and r['hb']=='15.0' and r['feedback']==str(fb) and r['admixture']==str(mix)) for k in [1,25]]
  losses.append(float(pair[0]['ca_ml_dl'])-float(pair[1]['ca_ml_dl']))
 display=np.array(losses)*(1000 if mode=='conductance_conditional' else 1)
 ax.plot(range(4),display,'o-',color=col);ax.set(xticks=range(4),xticklabels=['Fixed\n0%','Fixed\n3%','Feedback\n0%','Feedback\n3%'],ylabel='Content contrast (0.001 mL/dL)' if mode=='conductance_conditional' else 'Content contrast (mL/dL)',title=title,ylim=(0,max(display)*1.15));ax.tick_params(axis='x',labelsize=8)
savefig(fig,4)

