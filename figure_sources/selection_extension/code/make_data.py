"""Attributed additional published numbers; no simulated clinical observations."""
from pathlib import Path
import csv,json,hashlib
R=Path(__file__).resolve().parents[1];D=R/'data';D.mkdir(exist_ok=True)
def write(name,rows):
 with (D/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 rows=[]
 def add(study,unit,cohort,condition,baseline,n,q,pa,pv,paco2,hb,fi,ph,sat='',sv='',pap='',pvr='',vo2='',do2='',aa='',granularity='condition_mean',qmethod='source_method',pv_status='reported',hb_status='reported',doi=''):
  rows.append(dict(study=study,calibration_unit=unit,cohort_group=cohort,condition=condition,baseline=baseline,n_source=n,flow=q,flow_unit='L/min/m2' if study=='Melot1983' else 'L/min',Pa_native=pa,Pv_native=pv,PaCO2_native=paco2,pressure_unit='mmHg' if study=='Melot1983' else 'kPa',Hb_g_dl=hb,Hb_status=hb_status,FiO2=fi,PB_mmHg=760 if study=='Melot1983' else 717.7,PB_status='assumed_local_standard',pH=ph,Sa_percent=sat,Sv_percent=sv,mPAP_mmHg=pap,PVR_native=pvr,PVR_unit='dyn_s_cm5' if study=='Melot1983' else 'WU',M_index_ml_min_m2=vo2,DO2_index_ml_min_m2=do2,Aa_mmHg=aa,granularity=granularity,flow_method=qmethod,Pv_status=pv_status,doi=doi))
 # Haefliger et al. 2025, Table 2, adjusted condition estimates (SE, not SD).
 h=[('N',1,4.7,10.1,4.9,4.1,15.4,.21,7.45,94,67,39,7.10),('H',0,4.9,7.9,4.5,4.0,15.5,.15,7.47,90,63,40,7.16),('NS',0,5.1,10.1,4.9,4.2,15.4,.21,7.44,95,65,36,5.99),('HS',0,4.7,7.6,4.6,4.1,15.4,.15,7.46,90,63,38,6.75)]
 for c,b,q,pa,pv,pc,hb,fi,ph,sa,sv,pap,pvr in h:
  add('Haefliger2025','Haefliger2025','Zurich2021_2023',c,b,22,q,pa,pv,pc,hb,fi,ph,sa,sv,pap,pvr,granularity='adjusted_condition_mean',qmethod='direct_Fick_not_independent_of_oxygen',doi='10.3390/jcm14238407')
 # Lichtblau et al. 2025, Table 1; same trial family, different posture/subset.
 # Pv tension is NOT reported; inverse Sv is a declared secondary input scenario.
 l=[('air_rest',1,4.7,10.2,4.1,15.2,.21,7.45,94,69,39,7.0),('hypoxia_rest',0,4.4,8,4,15.3,.15,7.46,91,65,40,7.6),('air_exercise',0,7.2,8.7,4.4,15.6,.21,7.43,91,45,61,7.5),('hypoxia_exercise',0,7,6.3,4.1,15.7,.15,7.44,83,40,64,8.4)]
 for c,b,q,pa,pc,hb,fi,ph,sa,sv,pap,pvr in l:
  add('Lichtblau2025','Lichtblau2025','Zurich2021_2023',c,b,23,q,pa,'',pc,hb,fi,ph,sa,sv,pap,pvr,qmethod='direct_Fick_not_independent_of_oxygen',pv_status='unreported_inverse_Sv_secondary_only',doi='10.1164/rccm.202409-1737RL')
 # Melot et al. 1983 Table 1, duplicate determinations in two patients.
 # First five columns are patient 1 (two catheterizations); final two are patient 2.
 # Table 1 has no Hb: never reconstruct Hb from outcome-derived oxygen delivery.
 m=[('M01_visit1','baseline',1,1.9,46,27,27,7.46,74,1756,125,292,70),('M01_visit1','PGE1',0,2.8,56,35,29,7.45,53,787,131,465,58),('M01_visit1','hydralazine',0,3.3,63,37,25,7.46,72,963,136,568,56),('M01_visit2','baseline',1,2.0,42,25,26,7.48,75,1733,133,320,77),('M01_visit2','nifedipine',0,2.8,49,30,27,7.48,62,1002,147,470,69),('M02','baseline',1,2.1,81,31,25,7.51,53,992,133,382,42),('M02','nifedipine',0,3.4,85,37,22,7.55,54,641,130,615,41)]
 for unit,c,b,q,pa,pv,pc,ph,pap,pvr,vo,do,aa in m:
  add('Melot1983',unit,'Melot1983',c,b,1,q,pa,pv,pc,15,.21,ph,pap=pap,pvr=pvr,vo2=vo,do2=do,aa=aa,granularity='individual_duplicate_mean',qmethod='method_referenced_not_independence_established',hb_status='assumed_15_sensitivity_12_18',doi='10.1378/chest.83.2.203')
 write('additional_observations.csv',rows)
 sources=[dict(study='Haefliger2025',doi='10.3390/jcm14238407',url='https://www.mdpi.com/2077-0383/14/23/8407',location='Table 2',records=4,source_people=22,role='new_direct_Pv_condition_check',overlap='Same BASEC2020-02163 trial family as Lichtblau2025; supine vs semi-supine; do not pool participants'),dict(study='Lichtblau2025',doi='10.1164/rccm.202409-1737RL',url='https://pmc.ncbi.nlm.nih.gov/articles/PMC12005039/',location='Table 1',records=4,source_people=23,role='secondary_exercise_check_with_Sv_inverse',overlap='Same trial family as Haefliger2025; omit from primary Pa score because Pv not reported'),dict(study='Melot1983',doi='10.1378/chest.83.2.203',url='https://doi.org/10.1378/chest.83.2.203',location='Table 1 p204; author-uploaded full text',records=7,source_people=2,role='new_individual_perturbation_stress_test',overlap='M01 visits same patient; all grouped as one source; M01 has atrial right-to-left shunt')]
 write('additional_sources.csv',sources)
 (D/'input_sha256.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in D.glob('*.csv')},indent=2),encoding='utf-8')
 print('Additional source records',len(rows),'direct Pv',sum(x['Pv_status']=='reported' for x in rows))
if __name__=='__main__':main()
