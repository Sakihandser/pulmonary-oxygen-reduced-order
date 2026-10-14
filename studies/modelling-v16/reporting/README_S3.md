# Reporting supplement for Modelling v16

The CSV tables and provenance JSON summarize existing numerical results. No new
simulation or independent verification was performed for v16. Allocation gain is
flow-proportional minus fixed content; unmet demand is 250 minus stored systemic
use. Limiting compartments apply the assigned extraction ceilings to stored
content and are not measured organ deficits.

To export the four article figures at 600 dpi, install the dependencies recorded
in S2 and run `python tools/export_figures_v16.py` from this extracted directory.
This script reads local CSV files and writes Figures_v16/; it does not run models.
Arial is used for the article layout; font substitution can change pixels.
The separately supplied TIFF files are the submission figures.

Public source repository:
https://github.com/Sakihandser/pulmonary-oxygen-reduced-order
The fixed manuscript snapshot is specified in the article Data Availability Statement.
