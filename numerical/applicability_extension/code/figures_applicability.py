from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
R=Path(__file__).resolve().parents[1];F=R/'figures';F.mkdir(exist_ok=True)
with (R/'results/pair_compatibility.csv').open(encoding='utf-8') as f:rows=[r for r in csv.DictReader(f) if int(r['primary']) and float(r['FiO2'])<.3]
labels={'Carta2022':'Carta','Haefliger2025':'Haefliger','Lichtblau2021':'Lichtblau','M01_visit1':'M01 visit 1','M01_visit2':'M01 visit 2','M02':'M02'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
fig,(ax,bx)=plt.subplots(1,2,figsize=(11,6.1),gridspec_kw={'width_ratios':[1.6,1]},layout='constrained')
for i,r in enumerate(rows):
 s0=float(r['s0']);lo=float(r['s_lower'])-s0;hi=float(r['s_upper'])-s0;obs=float(r['s1'])-s0
 ax.plot([lo,hi],[i,i],lw=5,color='#a9c0d0',solid_capstyle='round')
 ax.scatter([obs],[i],marker='D',s=35,color='#207d6c' if int(r['compatible']) else '#ab4d38',zorder=3)
 bx.barh(i,float(r['Pa_halfwidth_mmHg']),color='#607f94',height=.55)
ax.set_yticks(range(len(rows)),[labels[r['unit']]+' / '+r['condition'].replace('acetazolamide','ACZ').replace('placebo','PL') for r in rows]);ax.invert_yaxis();bx.invert_yaxis()
ax.axvline(0,color='.5',lw=.8);ax.set_xlabel('Residual change from baseline');ax.set_title('(a) Allowed range and observed nominal change',loc='left',fontsize=11)
bx.set_yticks(range(len(rows)),['']*len(rows));bx.set_xlabel('Required common Pa halfwidth (mmHg)');bx.set_title('(b) Margin to compatibility',loc='left',fontsize=11);bx.set_xlim(0,3.1)
for aa in [ax,bx]:
 aa.set_ylim(len(rows)-.4,-.6);aa.spines[['top','right']].set_visible(False);aa.grid(axis='x',alpha=.2)
fig.savefig(F/'Fig9.png',dpi=300);fig.savefig(F/'Fig9.svg')
with Image.open(F/'Fig9.png') as im:im.convert('RGB').save(F/'Fig9.tif',compression='tiff_lzw',dpi=(300,300))
plt.close(fig)
