from pathlib import Path
import subprocess,sys,shutil
R=Path(__file__).resolve().parent
for rel in ['reproduce_figures.py','experimental_validation/code/figures.py','experimental_validation/code/review_figures.py','selection_extension/code/figures.py','review_extension/code/figures_review.py','coupling_extension/code/figures_coupling.py','heterogeneity_extension/code/figures_heterogeneity.py','adequacy_extension/code/figures_adequacy.py','applicability_extension/code/figures_applicability.py','certainty_extension/code/figures_certainty.py']:
    subprocess.run([sys.executable,str(R/rel)],cwd=R,check=True)
D=R/'raw_submission_figures';D.mkdir(exist_ok=True)
for i in range(1,8):
    src=R/('review_extension' if i<=2 else ('heterogeneity_extension' if i==7 else ('coupling_extension' if i==6 else 'selection_extension')))/'figures'/f'Fig{i}.tif'
    shutil.copy2(src,D/f'Fig{i}.tif')
for i in range(1,13):
    src=R/('figures' if i<=7 else 'experimental_validation/figures')/f'Fig{i}.tif'
    shutil.copy2(src,D/f'FigB{i}.tif')
for n in ['S1','S2']:shutil.copy2(R/'figures'/f'Fig{n}.tif',D/f'Fig{n}.tif')
shutil.copy2(R/'adequacy_extension/figures/Fig8.tif',D/'Fig8.tif')
shutil.copy2(R/'applicability_extension/figures/Fig9.tif',D/'Fig9.tif')
shutil.copy2(R/'certainty_extension/figures/Fig10.tif',D/'Fig10.tif')
shutil.copy2(R/'certainty_extension/figures/Fig11.tif',D/'Fig11.tif')
from PIL import Image
C=R/'submission_figures';C.mkdir(exist_ok=True)
for path in D.glob('*.tif'):
    with Image.open(path) as im:
        im.load();copy=im.copy();dpi=im.info.get('dpi',(300,300))
    copy.save(C/path.name,compression='tiff_lzw',dpi=dpi)


subprocess.run([sys.executable,str(R/'oxygen_decisions/plot_decisions.py')],check=True)
shutil.copy2(R/'oxygen_decisions/Fig12.tif',C/'Fig12.tif')
subprocess.run([sys.executable,str(R/'oxygen_decisions/external_direction.py')],check=True)
shutil.copy2(R/'oxygen_decisions/Fig13.tif',C/'Fig13.tif')
print('Reproduced all 27 figures')
