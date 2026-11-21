"""Figure 12: empirical ranking variability and observation-design witnesses."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def read(n):
    with (ROOT/'results'/n).open(encoding='utf-8') as f:return list(csv.DictReader(f))
def main():
    b=read('review_bootstrap.csv');w=read('review_nonidentifiability_witness.csv')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'oxygen-review-design'})
    fig,ax=plt.subplots(1,3,figsize=(9.5,3.5),layout='constrained')
    groups=['all_pulmonary_vascular','historical_primary_PH','recurrent_emboli']
    labels=['All seven','Idiopathic PH (4)','Emboli (3)']
    for i,g in enumerate(groups):
        r=next(x for x in b if x['subset']==g and x['law']=='saturation')
        v=float(r['admixture_minus_law_RMSE_mmHg']);lo=float(r['empirical_delta_q025']);hi=float(r['empirical_delta_q975'])
        ax[0].plot([lo,hi],[i,i],color='#287A95',lw=2);ax[0].plot(v,i,'o',color='#173B57')
    ax[0].axvline(0,color='gray',ls='--',lw=1);ax[0].set_yticks(range(3),labels);ax[0].invert_yaxis()
    ax[0].set_xlabel('RMSE difference (mmHg)\nAdmixture minus saturation');ax[0].set_title('A  Empirical resampling',loc='left',fontsize=10.5)
    same=[x for x in w if x['family']=='same_arterial_and_venous']
    q=[float(x['Q']) for x in same]
    ax[1].plot(q,[float(x['Sa_percent']) for x in same],'o-',label='Arterial');ax[1].plot(q,[float(x['Sv_percent']) for x in same],'s-',label='Mixed venous')
    ax[1].set_ylim(0,100);ax[1].set_xlabel('Q (L/min)\nM/Q fixed');ax[1].set_ylabel('Saturation (%)');ax[1].legend(frameon=False,loc='lower right');ax[1].set_title('B  Identical oxygen',loc='left',fontsize=10.5)
    ax[2].plot(q,[float(x['DO2_ml_min']) for x in same],'o-',color='#A75032')
    ax[2].set_xlabel('Q (L/min)\nM/Q fixed');ax[2].set_ylabel('Oxygen delivery (mL/min)');ax[2].set_title('C  Different delivery',loc='left',fontsize=10.5)
    out=ROOT/'figures';out.mkdir(exist_ok=True)
    for ext in ['png','tif','svg']:
        kwargs={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tif' else {}
        fig.savefig(out/f'Fig12.{ext}',dpi=600,**kwargs)
    plt.close(fig)
if __name__=='__main__':main()
