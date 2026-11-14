# Targeted source expansion and provenance

Search cutoff: 10 November 2026. This was a targeted search for usable paired gas-exchange tables, not an exhaustive systematic review. Searches combined pulmonary hypertension, mixed venous oxygen, arterial oxygen, cardiac output, exercise, hypoxia, sildenafil and nifedipine, and followed references of the existing sources.

Included numerical sources:

| Source | Retrieved evidence | Use and dependency |
| --- | --- | --- |
| Häfliger et al. 2025, DOI 10.3390/jcm14238407 | Publisher full text and Europe PMC XML PMC12693256, Table 2 | Four adjusted condition estimates, mean and SE. Measured Pv. Modern study family BASEC2020-02163. Flow is direct Fick. |
| Lichtblau et al. 2025, DOI 10.1164/rccm.202409-1737RL | Europe PMC XML PMC12005039, Table 1 | Four rest/exercise means. No measured Pv; inverse-Sv extension is secondary. Shares the modern family; participant counts are not added. |
| Mélot et al. 1983, DOI 10.1378/chest.83.2.203 | Author-uploaded full text on ResearchGate, Table 1 p. 204 | Seven duplicate-mean states in two people, three catheterization baselines. M01 has atrial right-to-left shunting. Hb absent. Indexed flow retained. Flow independence not established. |

Mélot author source: https://www.researchgate.net/publication/16377122_Effect_of_nifedipine_on_ventilationperfusion_matching_in_primary_pulmonary_hypertension

Additional sources assessed but not entered into the primary response score:

1. Carta et al. 2021, DOI 10.3389/fmed.2021.776956: cardiac-output comparison from an existing cohort; no complete new paired Pa/Pv table for the stated task.
2. Otulana and Higenbottam 1988, DOI 10.1183/09031936.93.01080732: incomplete accessible paired inputs; abstract alone insufficient for extraction.
3. Dantzker et al. 1984, DOI 10.1164/arrd.1984.130.3.412: complete response table unavailable in the accessed source.
4. Nagata et al. 2022, PMC9301830: relevant mixed-venous/prognostic observations, but individual data require access authorization/request and did not supply this paired prediction dataset.
5. Exercise PAH report PMC10382004: aggregate saturation information without a complete measured Pa/Pv input set.
6. Sildenafil exercise report DOI 10.1183/23120541.01655-2025: overlaps the modern trial family; complete accessible table not obtained.

No participant-level data request, author contact, credentialed access or external submission was performed. Source full texts and retrieval failures remain in the local project audit archive, not in the distributed numerical package. The code package contains attributed selected factual observations, extraction locations and original code.

Native kPa values use 7.500616827 mmHg/kPa. Modern PB=717.7 mmHg and historical PB=760 mmHg are local-standard assumptions. PH2O=47.1 mmHg. The full alveolar gas equation uses measured PaCO2 and assumed R, never held-out arterial Pa or A-a. Missing Hb is varied jointly with R and the blood-curve parameters. Reported pH is retained but no patient-specific Bohr correction is fitted. Modern aggregate SE values are not treated as SD or used to invent individual observations.

All outcome-specific effective admixture estimates in state_diagnostics.csv are reconstructions, not predictive inputs. Direct Fick flow and the same oxygen measurements cannot independently validate a mass balance. Historical original-author counterfactuals are not presented as newly discovered physiology.
