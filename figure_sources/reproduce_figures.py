from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
for name in ['plot_main_figures.py','plot_sensitivity_figures.py','plot_calibration_maps.py']:
    subprocess.run([sys.executable,str(root/'tools'/name)],check=True)
