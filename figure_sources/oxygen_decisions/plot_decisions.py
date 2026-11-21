from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
def read(path):
    with path.open(encoding="utf-8",newline="") as stream:
        return list(csv.DictReader(stream))

def plot():
    fig, axes = plt.subplots(1,3,figsize=(10.2,4.25),layout="constrained",width_ratios=[1.05,1,1])
    rows=[r for r in read(HERE/'observed.csv') if float(r['residual_halfwidth'])==.01
          and float(r['relative_flow_ratio_halfwidth'])==.10]
    hill = [r for r in rows if r['conversion']=='hill']
    sat = {(r['unit'],r['condition']):r for r in rows if r['conversion']=='reported_saturation'}
    labels=['Carta hypoxia','Häfliger H','Häfliger NS','Häfliger HS','ACZ air','Placebo hypoxia','ACZ hypoxia','M01 PGE1','M01 hydralazine','M01 nifedipine','M02 nifedipine']
    for i,r in enumerate(hill):
        axes[0].plot(float(r['residual_gap']),i,'o',color='#175676',label='Hill' if i==0 else None)
        pair = sat.get((r['unit'],r['condition']))
        if pair:
            axes[0].plot(float(pair['residual_gap']),i,'+',color='#bb5429',label='Reported saturation' if i==0 else None)
    axes[0].set_yticks(range(len(hill)),labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel('Minimum residual discrepancy')
    axes[0].legend(fontsize=7)
    data=[r for r in read(HERE/'margin_grid.csv') if r['model']=='A' and float(r['content_halfwidth'])==.003]
    for shift in [0,.001,.005]:
        line=[r for r in data if float(r['target_up'])==shift]
        axes[1].plot([100*float(r['demand_up_fraction']) for r in line],[float(r['lower_margin']) for r in line],marker='o',label=f'ΔT = {shift:g}')
    axes[1].axhline(0,color='black',lw=.8)
    axes[1].set(xlabel='Demand increase (%)',ylabel='Lower Ca − target (mL/dL)')
    axes[1].legend(fontsize=7)
    for r in read(HERE/'redistribution.csv'):
        axes[2].plot(float(r['shift'])*100,float(r['Ca'])-float(r['target']),'o',color='#175676')
    axes[2].axhline(0,color='black',lw=.8)
    axes[2].set(xlabel='Flow transferred (percentage points)',ylabel='Ca − target (mL/dL)')
    for ax,title in zip(axes,['A   Published PH pairs','B   Joint error budget (A)','C   Redistribution (A)']):
        ax.set_title(title,fontsize=10)
        ax.tick_params(labelsize=8)
    fig.savefig(HERE/'Fig12.svg')
    fig.savefig(HERE/'Fig12.tif',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)

if __name__ == "__main__":
    plot()
