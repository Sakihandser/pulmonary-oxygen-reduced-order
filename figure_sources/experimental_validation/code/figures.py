from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from validate import ROOT,OUT,LAWS
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
def read(name):
    with (OUT/name).open(encoding='utf-8') as f:return list(csv.DictReader(f))
labels=['Saturation','Content','Constant D','Pressure','Equilibrium + mixing']
colors=['#386cb0','#66a61e','#d95f02','#7570b3','#1b9e77']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
def save(fig,name):
    for ext in ['png','tif','svg']:fig.savefig(FIG/(name+'.'+ext),dpi=300,bbox_inches='tight')
    plt.close(fig)

def main():
    pp=read('predictions.csv');rr=[r for r in pp if r['protocol']=='baseline_frozen']
    keys=list(dict.fromkeys((r['study'],r['condition']) for r in rr))
    standard=[k for k in keys if k[1]!='hyperoxia'];hyper=[k for k in keys if k[1]=='hyperoxia']
    fig,axs=plt.subplots(1,2,figsize=(10,5.6),gridspec_kw={'width_ratios':[4,1]},layout='constrained')
    for ax,kk,title,limit in zip(axs,[standard,hyper],['A  Air / hypoxia / exercise','B  Pure oxygen'],[30,350]):
        data=np.array([[float(next(r['error_mmHg'] for r in rr if (r['study'],r['condition'],r['law'])==(s,c,l))) for l in LAWS] for s,c in kk])
        im=ax.imshow(data.T,cmap='RdBu_r',vmin=-limit,vmax=limit,aspect='auto')
        ax.set_xticks(range(len(kk)),[s.replace('2022','22').replace('2021','21').replace('2010','10').replace('2011','11')+'\n'+c.replace('_',' ') for s,c in kk],rotation=65,ha='right',fontsize=7)
        ax.set_yticks(range(5),labels if len(kk)>1 else ['']*5);ax.set_title(title,loc='left',fontsize=11)
        for i in range(5):
            for j in range(len(kk)):ax.text(j,i,f'{data[j,i]:.1f}',ha='center',va='center',fontsize=7,color='white' if abs(data[j,i])>limit*.55 else 'black')
        cb=fig.colorbar(im,ax=ax,orientation='horizontal',pad=.25,shrink=.9);cb.set_label('Predicted minus observed PaO2 (mmHg)',fontsize=8)
    save(fig,'Fig8')
    fig,axs=plt.subplots(1,2,figsize=(9,4.1),layout='constrained')
    for l,label,col in zip(LAWS,labels,colors):
        r=[r for r in pp if r['protocol']=='Carta_frozen_external' and r['law']==l]
        axs[0].scatter([float(x['observed_Pa']) for x in r],[float(x['predicted_Pa']) for x in r],label=label,color=col,s=25,alpha=.85)
    axs[0].plot([35,130],[35,130],ls='--',color='gray');axs[0].set(xlim=(35,130),ylim=(35,130),xlabel='Observed PaO2 (mmHg)',ylabel='Predicted PaO2 (mmHg)',title='A  Frozen transfer to seven patients')
    axs[0].legend(fontsize=7,frameon=False)
    mm=read('metrics.csv');x=np.arange(5)
    for offset,protocol,label in [(-.18,'leave_one_patient_out','Leave one patient out'),(.18,'Carta_frozen_external','Independent cohort transfer')]:
        values=[float(next(r['RMSE_mmHg'] for r in mm if r['study']=='Dantzker1979' and r['protocol']==protocol and r['law']==l)) for l in LAWS]
        axs[1].bar(x+offset,values,width=.36,label=label,color='#457b9d' if offset<0 else '#a8dadc')
    axs[1].set_xticks(x,labels,rotation=35,ha='right',fontsize=8);axs[1].set(ylabel='RMSE (mmHg)',title='B  Errors on held-out patients');axs[1].legend(frameon=False,fontsize=7)
    save(fig,'Fig9')
    rr=read('network_scenarios.csv');sh=read('network_admixture_sensitivity.csv')
    fig,axs=plt.subplots(2,2,figsize=(9,6.8),layout='constrained')
    pa0=min(float(r['PA']) for r in rr)
    for s in sorted(set(r['effective_admixture'] for r in sh),key=float):
        r=[r for r in sh if r['effective_admixture']==s]
        axs[0,0].plot([float(x['PVR_WU']) for x in r],[float(x['Sa_percent']) for x in r],'o-',label=f'Effective mixing {100*float(s):.1f}%')
    axs[0,0].set(title='A  Mixing changes the oxygen response',ylabel='SaO2 (%)');axs[0,0].legend(fontsize=7,frameon=False)
    for hb in [10,14]:
        r=[r for r in rr if r['law']=='equilibrium_admixture' and r['contrast']=='1' and float(r['Hb'])==hb and float(r['VO2'])==250 and float(r['PA'])==pa0 and r['feasible']=='True']
        axs[0,1].plot([float(x['PVR_WU']) for x in r],[float(x['DO2_ml_min']) for x in r],'o-',label=f'Hb {hb} g/dL')
    axs[0,1].set(title='B  Delivery depends strongly on flow and Hb',ylabel='Systemic O2 delivery (mL/min)');axs[0,1].legend(frameon=False,fontsize=8)
    for pa in [pa0,150.]:
        r=[r for r in rr if r['law']=='equilibrium_admixture' and r['contrast']=='1' and float(r['Hb'])==14 and float(r['VO2'])==250 and float(r['PA'])==pa]
        axs[1,0].plot([float(x['PVR_WU']) for x in r],[float(x['Sa_percent']) for x in r],'o-',label=f'PA {pa:.1f} mmHg')
    axs[1,0].set(title='C  Oxygen boundary is a scenario input',ylabel='SaO2 (%)');axs[1,0].legend(frameon=False,fontsize=8)
    for law,label,col in zip(LAWS,labels,colors):
        r=[r for r in rr if r['law']==law and r['contrast']=='1' and float(r['Hb'])==14 and float(r['VO2'])==250 and float(r['PA'])==pa0 and r['feasible']=='True']
        axs[1,1].plot([float(x['PVR_WU']) for x in r],[float(x['Sa_percent']) for x in r],'o-',label=label,color=col)
    axs[1,1].set(title='D  Competing closures remain different',ylabel='SaO2 (%)');axs[1,1].legend(frameon=False,fontsize=7)
    for ax in axs.flat:ax.set_xlabel('Pulmonary vascular resistance (WU)');ax.grid(alpha=.15)
    save(fig,'Fig10')
    controls=read('primary_flow_controls.csv');counter=read('primary_counterfactuals.csv')[1:]
    fig,axs=plt.subplots(2,2,figsize=(9.6,6.6),layout='constrained')
    for closure,label,color,style in [('inherited_cardiac_closure','Inherited Q(PVR)','#1b9e77','o-'),
                                      ('flow_clamped','Q held at 3.73 L/min','#7570b3','s--')]:
        rows=[r for r in controls if r['closure']==closure]
        for ax,key in [(axs[0,0],'Sa_percent'),(axs[0,1],'DO2_ml_min')]:
            ax.plot([float(r['PVR_WU']) for r in rows],[float(r[key]) for r in rows],style,label=label,color=color)
    axs[0,0].set(title='A  Pressure alone does not change SaO2',ylabel='SaO2 (%)',xlabel='PVR (WU)')
    axs[0,1].set(title='B  Flow closure drives the delivery decline',ylabel='O2 delivery (mL/min)',xlabel='PVR (WU)')
    for ax in axs[0]:ax.legend(fontsize=8,frameon=False)
    labels_cf=['Q reduced 25%','Hb 14 to 10 g/dL','PA 100.4 to 150 mmHg','Mixing 0.144 to 0.30','Demand 250 to 500 mL/min']
    for ax,key,title,xlabel in [(axs[1,0],'delta_Sa_pp','C  One factor changed at a time','Change in SaO2 (percentage points)'),
                                 (axs[1,1],'delta_DO2_ml_min','D  Different effects on oxygen delivery','Change in O2 delivery (mL/min)')]:
        vals=[float(r[key]) for r in counter];y=np.arange(len(vals))
        ax.barh(y,vals,color=['#1b9e77' if v>=0 else '#457b9d' for v in vals],height=.6)
        ax.set_yticks(y,labels_cf,fontsize=8);ax.invert_yaxis();ax.axvline(0,color='gray',lw=.6)
        ax.set(title=title,xlabel=xlabel);ax.margins(x=.22)
        for i,v in enumerate(vals):ax.annotate(f'{v:+.1f}',(v,i),xytext=(3 if v>=0 else -3,0),textcoords='offset points',ha='left' if v>=0 else 'right',va='center',fontsize=8)
    for ax in axs.flat:ax.grid(alpha=.15,axis='x')
    save(fig,'Fig11')

if __name__=='__main__':main()
