# 输出与校准信息字典

| 字段或常数 | 单位 | 定义与来源 |
| --- | --- | --- |
| k | 无量纲 | 上游动脉乘子对比；1 与 25 是归档场景，不是 PAH 严重程度 |
| allocation fixed | 无量纲 | CT 空气容积权重；不等于实测毛细血管血容量份额 |
| allocation flow | 无量纲 | 当下流量份额分配容量；规定的基准，不是测得募集规律 |
| hb | g/dL | 15 或 13.8；后者为来源组 Hb 均值，非与网络匹配的受试者数据 |
| cv_ml_dl | mL/dL | 固定 Pv=40 经 Hill 转换，或由系统反馈求得入口含量 |
| mixed_end_ml_dl | mL/dL | Cm，末端交换含量按区域流量混合 |
| ca_ml_dl | mL/dL | Ca=(1−s)Cm+sCv，s=0 或 0.03；规定的末端混合 |
| sa_fraction | 比例 0–1 | 共同 Hill 关系由 Ca 逆算；不是临床共氧合实测值 |
| do2_ml_min | mL/min | 10QCa，总氧输送率；勿与氧交换导度 D 混淆 |
| uptake_ml_min | mL/min | 10Q(Ca−Cv)，与规定混合条件相容的实际摄氧 |
| systemic_used_ml_min | mL/min | 反馈情景的器官实际利用氧；受提取上限限制，不保证达到250 |
| redistribution_content_loss | mL/dL | 指定条件内 Ca(k=1)−Ca(k=25)；共同无混合表中 Ca=Cm |
| redistribution_sa_loss_pp | 百分点 | 100[Sa(k=1)−Sa(k=25)]；不是百分比相对下降 |
| d_total | mL/min/mmHg | 常导度模型的总 D；松弛模型没有恒定 D，该字段为空 |
| Vc relaxation | mL | 86，旧模型规定容量；非来源气体反演 Vc |
| tau | s | 0.30，规定松弛时间；未由独立多状态氧合数据识别 |
| Dmatch | mL/min/mmHg | 5.47329；一均匀床 Q=3.94644、Pv=40、PA=100、Hb=15 的模型间终点匹配 |
| Dgas Hb15 | mL/min/mmHg | 28.09197；DLCO=16.4、DLNO=48、thetaCO=.45、thetaNO=4.5、alpha=1.97 的条件映射 |
| Dgas Hb13.8 | mL/min/mmHg | 27.67245；thetaO2=1.8×13.8/15 的探索缩放；未重新校正原 DLCO |
| thetaO2 | mL O2/mL blood/min/mmHg | 常氧操作参考1.8；0.6 mmol/mL/min/kPa×22.414/7.50062≈1.793；非低氧通用常数 |
| thetaCO | mL CO/mL blood/min/mmHg | .45，诊断假设；未证明匹配来源 Hb、肺泡氧和校正链 |
| lambda O2/CO | 无量纲 | 1.23，膜导度换算；与 Hb 系数、总氧导度及血侧反应项区分 |

目前完成输入、转换和代码实现核对。反应系数生理适用性、原始 Hb 标准化、来源协方差与区域分配实测约束仍未完成。CSV 小数精度是计算精度，不是实验测量精度。

耦合预研 CSV 的 `venous_fraction_fixture` 只作用于左上肺叶：Rv_LUL=该比例×参考床总阻力。其他肺叶保持基准 Rv=0.2×参考床总阻力。`arterial_scale` 是外层求根确定的共同上游阻力比例，用于恢复规定总 PVR 和 Q；不是拟合得到的 PAH 重塑系数。`oxygen_scale_fixture` 只缩放氧导度，未缩放水力阻力，也不代表真实病理微血管床丢失。
