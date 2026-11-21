"""Final main Figures 1 and 2: matched comparison and evidence provenance."""
from pathlib import Path
import sys,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];F=R/'figures';F.mkdir(exist_ok=True)
sys.path.insert(0,str(R.parent/'selection_extension/code'))
import figures as original
def read(n):
 with (R/'results'/n).open(encoding='utf-8') as f:return list(csv.DictReader(f))
def save(fig,n):
 fig.savefig(F/f'Fig{n}.png',dpi=180,bbox_inches='tight')
 fig.savefig(F/f'Fig{n}.tif',dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'})
 fig.savefig(F/f'Fig{n}.svg',bbox_inches='tight');plt.close(fig)
def fig1():
 fig,ax=plt.subplots(figsize=(7,6.4));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
 original.box(ax,.03,.84,.94,.12,'37 published quantitative records from 8 studies\nPrimary response set: 11 changes and 6 baseline calibrations\nNo new patient data in the matched-scaling analysis',fontsize=9)
 original.box(ax,.03,.60,.45,.17,'Primary nonhyperoxic strata\nEarly Zurich: 4 changes, thermodilution\nRecent Zurich: 3 changes, direct Fick\nHistorical: 4 changes in 2 people\n3 source families; measured Pv','#DDF0EB',fontsize=8)
 original.box(ax,.54,.60,.43,.17,'Separate evidence strata\n1 pure-oxygen condition\n3 exercise checks using Sv inversion\n6 COPD/IPF changes\n7 cross-sectional rows',fontsize=8.5)
 for x in [.25,.75]:ax.annotate('',(x,.77),(x,.84),arrowprops={'arrowstyle':'->','color':'#516578'})
 original.box(ax,.03,.34,.94,.18,'Calibrate once at the permitted baseline\nCompare 4 exchange maps × 2 capacity-scaling rules\nKeep the original 5-law and 2-baseline comparison\nCurrent Q cancels in every proportional-capacity control',fontsize=10)
 ax.annotate('',(.25,.52),(.25,.60),arrowprops={'arrowstyle':'->','color':'#516578'})
 original.box(ax,.03,.07,.94,.18,'Evaluate conditional predictions and explicit model scenarios\nEqual-family RMSE; modern and historical strata; 108 assumptions\nState-dependent curve shifts and compensation-boundary stress\nFick-Q dependence and aggregate-data limits remain explicit','#E8EEF6',fontsize=9.5)
 ax.annotate('',(.5,.25),(.5,.34),arrowprops={'arrowstyle':'->','color':'#516578'})
 fig.canvas.draw();renderer=fig.canvas.get_renderer()
 for patch,label in fig._label_boxes:
  b=patch.get_window_extent(renderer);t=label.get_window_extent(renderer)
  assert b.x0+4<=t.x0 and b.x1-4>=t.x1 and b.y0+4<=t.y0 and b.y1-4>=t.y1,label.get_text()
 save(fig,1)
def fig2():
 ms=read('matched_metrics.csv');pp=read('matched_predictions.csv')
 fig=plt.figure(figsize=(7,8.1),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1.6])
 laws=['saturation','content','pressure','constant_D'];labels=['Saturation','Content','Pressure','Conductance']
 for j,scope in enumerate(['air_hypoxia','modern']):
  ax=fig.add_subplot(gs[0,j]);yy=np.arange(4)
  for beta,off,color,label in [(0,-.17,'#007F7A','Capacity ∝ Q'),(1,.17,'#DA8A2E','Fixed capacity')]:
   vals=[float(next(r for r in ms if r['scope']==scope and r['model']==f'{law}_beta{beta}')['equal_family_RMSE']) for law in laws]
   ax.barh(yy+off,vals,height=.31,color=color,label=label)
   for i,val in enumerate(vals):ax.text(val+.12,i+off,f'{val:.2f}',fontsize=7,va='center')
  ax.set_yticks(yy,labels);ax.invert_yaxis();ax.set_xlim(0,18);ax.set_xlabel('Equal-family Pa RMSE (mmHg)');ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
  ax.set_title(('A  All 11 nonhyperoxic changes' if j==0 else 'B  Seven modern changes'),loc='left',fontweight='bold',fontsize=9)
  if j==0:ax.legend(loc='upper right',fontsize=7,frameon=False)
 ax=fig.add_subplot(gs[1,:]);rows=[r for r in pp if r['model']=='content_beta0' and float(r['FiO2'])<.3]
 mods=['saturation_beta0','content_beta0','saturation_beta1','content_beta1'];mat=[];names=[]
 for row in rows:
  unit=row['unit'].replace('M01_visit','M01 v').replace('Haefliger2025','Häfliger').replace('Lichtblau2021','Lichtblau').replace('Carta2022','Carta')
  names.append(unit+' / '+row['condition'].replace('acetazolamide','ACZ').replace('placebo','PBO'))
  mat.append([float(next(r for r in pp if (r['unit'],r['condition'],r['model'])==(row['unit'],row['condition'],model))['error']) for model in mods])
 im=ax.imshow(mat,cmap='RdBu_r',vmin=-22,vmax=22,aspect='auto');ax.set_yticks(range(len(names)),names)
 ax.set_xticks(range(4),['Saturation\ncapacity ∝ Q','Content\ncapacity ∝ Q','Saturation\nfixed capacity','Content\nfixed capacity'])
 for i,line in enumerate(mat):
  for j,val in enumerate(line):ax.text(j,i,f'{val:+.1f}',ha='center',va='center',fontsize=8,color='white' if abs(val)>14 else '#162B3B')
 ax.set_title('C  Saturation and content residuals in every primary condition',loc='left',fontweight='bold',fontsize=9)
 fig.colorbar(im,ax=ax,shrink=.7,label='Prediction − observation (mmHg)')
 save(fig,2)
if __name__=='__main__':fig1();fig2();print('Exported final Figures 1 and 2')
