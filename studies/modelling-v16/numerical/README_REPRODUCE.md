# Modelling v16 numerical supplement

This English distribution supplies the article inputs, complete 128-condition matrix, observation analyses, production code, independent endpoint functions, previously recorded mathematical checks and source provenance. Numerical results derive from the frozen v14/v15 study. Only paths, labels and documentation were translated; numerical cells and model arithmetic were preserved.

From the extracted supplement root:

```
python -m pip install -r requirements-lock.txt
python reproduce.py
```

The runner checks the current English-distribution hashes before recomputing the production matrix, reconstructing hydraulics and independently checking endpoints. It needs no network access, participant records or third-party article full text. Run in a disposable copy because reproduction changes outputs. Use stored numerical values and recorded integration tolerances for comparisons; output hashes containing runtime metadata are not numerical agreement tests. The earlier reproduction was performed in the recorded environment, not a newly isolated dependency environment. No model, test or numerical verification was rerun during this language revision.

The relative English phase directories must remain intact. The current entry point is 11_diagnostic_reproduction/code/reproduce_diagnostic.py. The historical figure script tools/redraw_figures_v14.py reads 10_diagnostic_revision/results and exports its historical figure layout. For the current 600 dpi article figures, use the portable script and CSV files in S3; the four submission TIFF files are supplied separately. Arial is used and font substitution may change pixels.

Public repository: https://github.com/Sakihandser/pulmonary-oxygen-reduced-order. The manuscript identifies the fixed publication snapshot. release_sha256.json covers this translated distribution; translation_provenance.json retains source hashes and explains the revision. Old run manifests are historical records, not current distribution hash lists. Independent physiological oxygenation and PAH validation remain absent.
