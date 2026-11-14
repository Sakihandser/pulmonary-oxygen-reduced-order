"""Five publication figures exported from frozen result CSVs."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.colors import ListedColormap,BoundaryNorm
R=Path(__file__).resolve().parents[1];F=R/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','savefig.facecolor':'white'})
C={'equilibrium_admixture':'#007F7A','content':'#DA8A2E','saturation':'#7669AD','constant_D':'#C95454','pressure':'#407AA1','baseline_Pa':'#A6ADB6','baseline_Aa':'#6C747D','flow_adaptive_nonhyperoxia':'#233E55'}
L={'equilibrium_admixture':'Effective admixture','content':'Content relaxation','saturation':'Saturation relaxation','constant_D':'Fixed conductance','pressure':'Pressure relaxation','baseline_Pa':'Baseline Pa held','baseline_Aa':'Baseline A-a held','flow_adaptive_nonhyperoxia':'Fitted flow exponent (held out)'}
def read(name):
 with (R/'results'/name).open(encoding='utf-8') as f:return list(csv.DictReader(f))
def save(fig,n):
 if n==1:
  fig.canvas.draw();renderer=fig.canvas.get_renderer()
  for patch,label in fig._label_boxes:
   bounds=patch.get_window_extent(renderer);text_bounds=label.get_window_extent(renderer)
   assert bounds.x0+4<=text_bounds.x0 and bounds.x1-4>=text_bounds.x1 and bounds.y0+4<=text_bounds.y0 and bounds.y1-4>=text_bounds.y1,('Flowchart text overflow',label.get_text())
 fig.savefig(F/f'Fig{n}.png',dpi=180,bbox_inches='tight');fig.savefig(F/f'Fig{n}.tif',dpi=300,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'});fig.savefig(F/f'Fig{n}.svg',bbox_inches='tight');plt.close(fig)
def box(ax,x,y,w,h,text,color='#ECF2F5',fontsize=9):
 patch=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.012',facecolor=color,edgecolor='#98A9B6',linewidth=.7);ax.add_patch(patch)
 label=ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=fontsize)
 if not hasattr(ax.figure,'_label_boxes'):ax.figure._label_boxes=[]
 ax.figure._label_boxes.append((patch,label))
def fig1():
 fig,ax=plt.subplots(figsize=(7,6.4));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
 box(ax,.03,.84,.94,.12,'37 published quantitative records from 8 studies\n7 earlier individual rows + 7 states in 2 additional patients\n+ 23 aggregate records',fontsize=9)
 box(ax,.03,.61,.45,.16,'Primary response comparison\n11 nonhyperoxic changes\n3 cohort families; measured Pv\n6 separate baseline calibrations','#DDF0EB',fontsize=8.5)
 box(ax,.54,.61,.43,.16,'Separate evidence strata\n1 pure-oxygen condition\n3 exercise checks using Sv inversion\n6 COPD/IPF changes\n7 cross-sectional rows',fontsize=8.5)
 ax.annotate('',(.25,.77),(.25,.84),arrowprops={'arrowstyle':'->','color':'#516578'});ax.annotate('',(.75,.77),(.75,.84),arrowprops={'arrowstyle':'->','color':'#516578'})
 box(ax,.03,.35,.94,.18,'Fit each candidate once to the permitted baseline\nFreeze the coefficient and predict changed conditions\nFive exchange models + two baseline comparators\nShared flow exponent: leave an entire cohort family out',fontsize=10)
 ax.annotate('',(.25,.53),(.25,.61),arrowprops={'arrowstyle':'->','color':'#516578'})
 box(ax,.03,.08,.94,.18,'Score and explain\nEqual-family RMSE + MAE; 108 input scenarios; group deletion\nShow residuals, capacity scaling and flow–admixture compensation\nFick-derived Q does not independently validate the closed balance','#E8EEF6',fontsize=10)
 ax.annotate('',(.5,.26),(.5,.35),arrowprops={'arrowstyle':'->','color':'#516578'})
 save(fig,1)
def fig2():
 ms=read('response_metrics.csv');pp=read('response_predictions.csv');summary=json.loads((R/'results/mechanism_summary.json').read_text())
 models=['equilibrium_admixture','flow_adaptive_nonhyperoxia','saturation','content','pressure','constant_D','baseline_Pa','baseline_Aa']
 vals=[summary['nonhyperoxia_nested_score']['equal_family_RMSE'] if m=='flow_adaptive_nonhyperoxia' else float(next(r for r in ms if r['scope']=='primary_PH_nonhyperoxia' and r['model']==m)['equal_family_RMSE']) for m in models]
 fig,(ax,bx)=plt.subplots(2,1,figsize=(7,8.3),gridspec_kw={'height_ratios':[1,1.65]},layout='constrained')
 ax.barh(np.arange(len(models)),vals,color=[C[m] for m in models],height=.68);ax.invert_yaxis();ax.set_yticks(range(len(models)),[L[m] for m in models]);ax.set_xlabel('Equal-family Pa RMSE (mmHg)');ax.set_xlim(0,20)
 for i,val in enumerate(vals):ax.text(val+.2,i,f'{val:.2f}',va='center',fontsize=8)
 ax.set_title('A  Baseline-calibrated nonhyperoxic response prediction',loc='left',fontweight='bold');ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
 mods=['equilibrium_admixture','saturation','content','pressure','constant_D']
 rr=[r for r in pp if r['primary_PH']=='1' and float(r['FiO2'])<.3 and r['model']=='equilibrium_admixture']
 labels=[];mat=[]
 short={'Carta2022':'Carta','Lichtblau2021':'Lichtblau 2021','Haefliger2025':'Häfliger','Melot1983':'Mélot'}
 for r in rr:
  labels.append((r['unit'].replace('M01_visit','M01 v') if r['study']=='Melot1983' else short[r['study']])+' / '+r['condition'].replace('acetazolamide','ACZ').replace('placebo','PBO'))
  mat.append([float(next(x for x in pp if (x['unit'],x['condition'],x['model'])==(r['unit'],r['condition'],m))['error']) for m in mods])
 im=bx.imshow(mat,cmap='RdBu_r',vmin=-25,vmax=25,aspect='auto');bx.set_yticks(range(len(rr)),labels);bx.set_xticks(range(len(mods)),['Admixture','Saturation','Content','Pressure','Fixed D']);bx.tick_params(axis='x',rotation=20)
 for i,line in enumerate(mat):
  for j,value in enumerate(line):bx.text(j,i,f'{value:+.1f}',ha='center',va='center',fontsize=8,color='white' if abs(value)>17 else '#162B3B')
 bx.set_title('B  Every scored change remains visible',loc='left',fontweight='bold');fig.colorbar(im,ax=bx,shrink=.6,label='Prediction − observation (mmHg)');save(fig,2)
def fig3():
 fig,axs=plt.subplots(2,2,figsize=(7,6.4),layout='constrained');ms=read('response_metrics.csv')
 mods=['equilibrium_admixture','content','pressure','constant_D','saturation']
 ax=axs[0,0];x=np.arange(5)
 for offset,scope,col,label in [(-.18,'primary_PH_nonhyperoxia','#007F7A','Air / hypoxia'),(.18,'primary_PH_all','#D98A2D','Including pure oxygen')]:
  vals=[float(next(r for r in ms if r['scope']==scope and r['model']==m)['equal_family_RMSE']) for m in mods];ax.bar(x+offset,vals,width=.34,color=col,label=label)
 ax.set_xticks(x,['Mix','Content','Pressure','D','Sat']);ax.set_ylabel('Pa RMSE (mmHg)');ax.legend(fontsize=7);ax.set_title('A  Regime composition changes ranking',loc='left',fontweight='bold')
 ax=axs[0,1];p=read('flow_exponent_profile.csv');n=read('nonhyperoxia_exponent_profile.csv')
 ax.plot([float(r['beta']) for r in p],[float(r['all_primary_RMSE']) for r in p],color='#D98A2D',label='Including pure oxygen');ax.plot([float(r['beta']) for r in n],[float(r['RMSE']) for r in n],color='#007F7A',label='Air / hypoxia');ax.axvline(0,c='#007F7A',ls=':',lw=1);ax.axvline(1,c='#D98A2D',ls=':',lw=1);ax.set(xlabel='Flow exponent β',ylabel='All-data fit RMSE (mmHg)');ax.legend(fontsize=7);ax.set_title('B  Fit profiles are not validation',loc='left',fontweight='bold')
 ax=axs[1,0];f=read('nonhyperoxia_exponent_folds.csv');labels=['Early Zurich','Modern Zurich','Mélot'];ax.scatter([float(r['beta']) for r in f],range(3),color='#007F7A',s=50);ax.axvline(0,c='#007F7A',ls=':');ax.axvline(1,c='#D98A2D',ls='--');ax.set_yticks(range(3),labels);ax.set(xlim=(-.3,1.15),xlabel='β fitted without the labelled family');ax.invert_yaxis();ax.set_title('C  Fits with one family omitted',loc='left',fontweight='bold')
 ax=axs[1,1];q=np.linspace(.6,2.,100);s0=.15
 for beta,col,label in [(0,'#007F7A','β = 0: constant admixture'),(1,'#D98A2D','β = 1: fixed content capacity')]:ax.plot(q,np.exp(np.log(s0)*(1/q)**beta),color=col,label=label,lw=2)
 ax.set(xlabel='Flow / baseline flow',ylabel='Effective residual fraction s',ylim=(0,.42));ax.legend(fontsize=7,loc='upper left');ax.set_title('D  Flow scaling separates the models',loc='left',fontweight='bold');save(fig,3)
def fig4():
 rr=read('response_decomposition.csv');fig,(ax,bx)=plt.subplots(2,1,figsize=(7,7.4),layout='constrained',gridspec_kw={'height_ratios':[1.45,1]})
 labels=[(r['unit'].replace('M01_visit','M01 v') if r['study']=='Melot1983' else {'Carta2022':'Carta','Lichtblau2021':'Lichtblau','Haefliger2025':'Häfliger'}[r['study']])+' / '+r['condition'].replace('acetazolamide','ACZ').replace('placebo','PBO') for r in rr]
 pos=np.zeros(len(rr));neg=np.zeros(len(rr))
 for field,col,label in [('alveolar_component','#4F8BC2','Alveolar input'),('venous_component','#007F7A','Venous input'),('Hb_component','#D89B3C','Hb'),('unexplained_residual','#A3A8AF','Unexplained residual')]:
  vals=np.array([float(r[field]) for r in rr]);left=np.where(vals>=0,pos,neg);ax.barh(range(len(rr)),vals,left=left,color=col,label=label,height=.7);pos+=np.maximum(vals,0);neg+=np.minimum(vals,0)
 ax.scatter([float(r['observed_delta_Pa']) for r in rr],range(len(rr)),c='black',marker='|',s=90,label='Observed net change',zorder=5);ax.axvline(0,c='#708090',lw=.6);ax.set_yticks(range(len(rr)),labels);ax.invert_yaxis();ax.set_xlabel('Change from permitted baseline Pa (mmHg)');ax.legend(ncol=1,fontsize=7,loc='upper right');ax.set_title('A  Frozen-model attribution with an explicit residual',loc='left',fontweight='bold')
 changes=[r for r in read('published_hemodynamic_changes.csv') if r['study']!='Lichtblau2025']
 for study,col,mark in [('Haefliger2025','#007F7A','s'),('Melot1983','#7669AD','o')]:
  z=[r for r in changes if r['study']==study];bx.scatter([float(r['delta_PVR']) for r in z],[float(r['delta_Pa']) for r in z],label='Häfliger adjusted means' if study.startswith('H') else 'Mélot individual states',c=col,marker=mark,s=45)
 for r in changes:
  name=r['condition'] if r['study'].startswith('H') else r['unit'].replace('M01_visit','M01v')+' '+r['condition'];bx.annotate(name,(float(r['delta_PVR']),float(r['delta_Pa'])),xytext=(4,4),textcoords='offset points',fontsize=7)
 bx.axhline(0,c='#999999',lw=.7);bx.axvline(0,c='#999999',lw=.7);bx.set(xlabel='Published PVR change (WU)',ylabel='Observed Pa change (mmHg)',ylim=(-23,22),xlim=(-14,2));bx.legend(fontsize=7,loc='lower left');bx.set_title('B  Published paired changes without a causal regression',loc='left',fontweight='bold');save(fig,4)
def fig5():
 fig,(ax,bx)=plt.subplots(1,2,figsize=(7,4.1),layout='constrained');b=read('compensation_boundaries.csv');s=np.array([float(r['s_after']) for r in b]);qa=np.array([float(r['Q_ratio_Pa_preserved']) for r in b]);qd=np.array([float(r['Q_ratio_DO2_preserved']) for r in b]);qf=np.array([float(r['Q_ratio_feasible']) for r in b])
 ax.fill_between(s,.6,qd,color='#F3D8D5');ax.fill_between(s,qd,qa,where=qa>=qd,color='#F5DDA5');ax.fill_between(s,np.maximum(qa,qd),2.6,color='#CEE8DF');ax.fill_between(s,qa,qd,where=qa<qd,color='#D8E5F3')
 ax.plot(s,qa,c='#007F7A',lw=2,label='Arterial oxygen preserved');ax.plot(s,qd,c='#A75A26',lw=2,label='Oxygen delivery preserved');ax.scatter([.20,.25],[1.2,1.2],c='black',s=24,zorder=5)
 ax.text(.258,1.35,'Delivery rises\narterial oxygen falls',fontsize=8,ha='center');ax.set(xlim=(.05,.35),ylim=(.6,2.6),xlabel='Effective admixture after change',ylabel='Flow / baseline flow');ax.legend(fontsize=7,loc='upper left');ax.set_title('A  Two different compensation thresholds',loc='left',fontweight='bold')
 # Q ratio = driving-pressure ratio / PVR ratio; fixed s increase 0.15 -> 0.20.
 u=np.linspace(.6,1.2,100);ra=1.416666666666666;rd=1.0198204755381168
 bx.fill_between(u,.3,u/ra,color='#CEE8DF');bx.fill_between(u,u/ra,u/rd,color='#F5DDA5');bx.fill_between(u,u/rd,1.25,color='#F3D8D5');bx.plot(u,u/ra,c='#007F7A',lw=2);bx.plot(u,u/rd,c='#A75A26',lw=2)
 bx.scatter([1,.8],[1/ra,.8/ra],c='#007F7A',s=25);bx.annotate('29.4% PVR reduction',(1,1/ra),xytext=(.76,.45),arrowprops={'arrowstyle':'-','lw':.7},fontsize=7);bx.annotate('43.5% reduction',(.8,.8/ra),xytext=(.62,.36),arrowprops={'arrowstyle':'-','lw':.7},fontsize=7)
 bx.set(xlim=(.6,1.2),ylim=(.3,1.25),xlabel='Driving pressure / baseline',ylabel='PVR / baseline PVR');bx.set_title('B  Pressure unloading shifts the requirement',loc='left',fontweight='bold')
 fig.supxlabel('Chosen reference: Q = 5 L/min, M = 250 mL/min, Hb = 14 g/dL, PA = 100 mmHg, s = 0.15\nPanel B: s increases to 0.20; demand and alveolar content fixed. These are model scenarios.',fontsize=8)
 save(fig,5)
def main():
 for f in [fig1,fig2,fig3,fig4,fig5]:f()
 print('Exported 5 figures as PNG, TIFF and editable SVG')
if __name__=='__main__':main()
