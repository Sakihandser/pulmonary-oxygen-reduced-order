from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(__file__).resolve().parents[1]
def read(n):
    with (R/'results'/n).open(encoding='utf-8') as f:return list(csv.DictReader(f))

def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'pulmonary-adequacy'})
    fig,axs=plt.subplots(2,2,figsize=(10,7.4),layout='constrained')
    grid=read('adequacy_grid.csv')
    for s,col in [(.05,'#2864a6'),(.1,'#208772'),(.2,'#ca6434'),(.4,'#994b9c')]:
        rr=[x for x in grid if float(x['s0'])==s and float(x['beta'])==1]
        axs[0,0].plot([float(x['r']) for x in rr],[float(x['radius']) for x in rr],label=f'Baseline s = {s:.2f}',color=col)
    axs[0,0].axhline(.2,color='#333333',ls='--',lw=1,label='Illustrative tolerance 0.2')
    axs[0,0].set(xlabel='Flow ratio Q/Q0',ylabel='Minimax error bound (mL/dL)',title='A  Baseline-only precision limit')
    axs[0,0].legend(fontsize=8,loc='upper left')
    examples=read('adequacy_examples.csv')
    wanted=['moderate_residual','weak_flow_dependence','narrow_support','proportional_capacity']
    v=[next(x for x in examples if x['case']==n) for n in wanted];y=np.arange(4)
    axs[0,1].barh(y-.17,[float(x['width']) for x in v],height=.32,label='Endpoint worst error',color='#b8cbe1')
    axs[0,1].barh(y+.17,[float(x['radius']) for x in v],height=.32,label='Minimax worst error',color='#2864a6')
    axs[0,1].set(yticks=y,yticklabels=['Fixed capacity','Beta = 0.25','Support [0.15, 0.25]','Proportional capacity'],xlabel='Error bound at Q/Q0 = 2 (mL/dL)',title='B  When a reduced predictor suffices')
    axs[0,1].axvline(.2,color='#333333',ls='--',lw=1);axs[0,1].legend(fontsize=8)
    demo=read('second_state_bounds.csv')
    for ax,name,label in [(axs[1,0],'finite_A','C  Finite model A'),(axs[1,1],'matched_B','D  Matched model B')]:
        rr=[x for x in demo if x['model']==name and x['method']=='baseline_only'];r=np.array([float(x['r']) for x in rr])
        ax.fill_between(r,[float(x['low']) for x in rr],[float(x['high']) for x in rr],color='#d9e0e7',label='Baseline-only sharp interval')
        rr=[x for x in demo if x['model']==name and x['method']=='two_exact']
        ax.fill_between(r,[float(x['low']) for x in rr],[float(x['high']) for x in rr],color='#77bca8',alpha=.85,label='Two exact states: sharp interval')
        rr=[x for x in demo if x['model']==name and x['method']=='two_noisy']
        ax.plot(r,[float(x['low']) for x in rr],color='#cf843a',ls=':',lw=1.4,label='Second s +/- 0.005: outer bound')
        ax.plot(r,[float(x['high']) for x in rr],color='#cf843a',ls=':',lw=1.4)
        ax.plot(r,[float(x['actual']) for x in rr],color='#153c46',lw=1.6,label='Specified finite structure')
        ax.axvline(2,color='#777777',ls='--',lw=.9)
        ax.set(xlabel='Flow ratio Q/Q0',ylabel='Arterial content (mL/dL)',title=label,ylim=(16,18.5))
        ax.legend(fontsize=7,loc='lower left')
    dest=R/'figures';dest.mkdir(exist_ok=True)
    for ext in ['png','tif','svg']:
        kwargs={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tif' else {}
        fig.savefig(dest/f'Fig8.{ext}',dpi=300,**kwargs)
    plt.close(fig)
if __name__=='__main__':main()
