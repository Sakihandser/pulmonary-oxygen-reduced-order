# Diagnostic article v15 offline reproducibility

S2 Data supplies the exact article inputs and every row of the 128-condition
diagnostic matrix, observation analyses, current figure files and source data.
Run from the extracted root:

```
python -m pip install -r requirements-lock.txt
python reproduce.py
```

The runner verifies all archived file hashes before recomputing the production
matrix, reconstructing hydraulics and independently checking every endpoint.
No network, participant records or third-party article text is required.
Observation results retain the original phase directories for provenance.
Figure source SVG/PNG/TIFF files are supplied. To redraw the final article figures,
run python tools/redraw_figures_v14.py after installing the locked dependencies.
Use Arial for the publication layout; font substitution may change pixels on
other machines. The numerical runner does not typeset manuscript DOCX.

源代码目录与输入路径保持相对结构。解压后从根目录运行 reproduce.py。
哈希校验在重算前执行；重算改变输出后，再次验证归档应重新解压。
保存的精确值和独立积分容差用于核对，勿以含运行时间的输出哈希比较结果。
本机同一环境下独立解压重跑；未另建全新依赖环境。


The v15 manuscript targets Modelling. All numerical inputs, algorithms, outputs,
figure pixels and the 128-condition matrix are unchanged from the frozen v14
numerical release. Legacy script and figure directory names retain v14 for
provenance. The future research agenda introduces no new validated outcomes.
