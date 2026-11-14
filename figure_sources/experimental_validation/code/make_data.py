"""Transcribe published observations; no simulated observations or digitized points.

PDF values were checked visually against original tables. XML source tables are
retained in the project evidence directory, not redistributed with this supplement.
"""
from pathlib import Path
import csv, json
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'; DATA.mkdir(parents=True,exist_ok=True)
rows=[]
def add(study,condition,n,q,pv,pa,paco,hb,**kw):
    r=dict(study=study,condition=condition,n=n,flow=q,flow_unit='L/min',Pv_mmHg=pv,
           Pa_mmHg=pa,PaCO2_mmHg=paco,Hb_g_dl=hb,Hb_status='reported',PA_mmHg='',
           FiO2=.21,PB_mmHg=760.,PH2O_mmHg=47.1,PA_source='alveolar_equation_assumed_R',
           Pa_SD_mmHg='',mPAP_mmHg='',PAWP_mmHg='',PVR_WU='',Sa_percent='',Sv_percent='',
           measured_shunt_fraction='',low_vq_fraction='',Ca_reported='',Cv_reported='',
           baseline=0,population='',table='',doi='',notes='')
    r.update(kw);rows.append(r)
kpa=7.500616827
for i,(cond,q,pv,pa,pco,hb,sd,mpap,wedge,pvr,sa,sv,ca,cv) in enumerate([
 ('air',5.4,4.6,9.4,4.5,14.,1.8,36,12,5.,92,65,17.6,12.3),
 ('hypoxia',5.4,4.3,7.,4.2,14.1,1.2,37,11,5.3,87,62,16.5,11.8),
 ('hyperoxia',5.1,6.,60.5,4.4,14.1,12.9,32,11,4.5,98,78,19.9,14.8)]):
    add('Carta2022',cond,149,q,pv*kpa,pa*kpa,pco*kpa,hb,FiO2=[.21,.15,1.][i],PB_mmHg=717.7,
        baseline=int(i==0),Pa_SD_mmHg=sd*kpa,mPAP_mmHg=mpap,PAWP_mmHg=wedge,PVR_WU=pvr,
        Sa_percent=sa,Sv_percent=sv,Ca_reported=ca,Cv_reported=cv,population='PAH_79_CTEPH_70',
        table='Table 2; methods PB and R',doi='10.3389/fmed.2022.791423',
        notes='Available-case aggregate means; not 149 paired complete records. Hyperoxic shunt derived from outcome is not used as an independent predictor.')
for i,(cond,ci,pv,pa,pco,hb,sd,mpap,pvr,sa,sv) in enumerate([
 ('placebo_air',3.,4.6,8.6,4.4,14.4,1.6,39,5.3,91,63),
 ('acetazolamide_air',2.9,4.7,9.,4.4,14.4,1.8,40,5.5,91,64),
 ('placebo_hypoxia',3.,4.3,7.5,3.9,14.,1.6,39,5.2,88,62),
 ('acetazolamide_hypoxia',3.,4.3,7.3,4.1,14.1,1.4,39,5.,87,61)]):
    add('Lichtblau2021',cond,24,ci,pv*kpa,pa*kpa,pco*kpa,hb,flow_unit='L/min/m2',
        FiO2=.15 if 'hypoxia' in cond else .21,PB_mmHg=717.7,baseline=int(i==0),
        Pa_SD_mmHg=sd*kpa,mPAP_mmHg=mpap,PAWP_mmHg=11,PVR_WU=pvr,Sa_percent=sa,Sv_percent=sv,
        population='PAH_7_CTEPH_17',table='Tables 2 and 3',doi='10.3389/fmed.2021.681473',
        notes='Cardiac index used directly: fitted capacities are indexed, not total. Stable BSA cancels in within-cohort flow ratio. Possible overlap with Zurich Carta cohort; no independent pooled N.')
for i,(q,hb,pa,pco,pv,aa,mpap,sh,low) in enumerate([
 (2.3,14.9,90,22,30,28,50,.003,.003),(4.4,14.6,69,34,34,53,80,.006,.058),
 (5.7,14.7,85,24,45,34,70,.024,.071),(6.4,10.1,46,33,26,61,48,.083,.19),
 (5.1,15.6,85,29,44,42,42,.011,.089),(1.8,12.5,58,28,24,71,65,.061,.061),
 (6.7,13.3,65,24,33,58,49,0.,.17)]):
    add('Dantzker1979','D%02d'%(i+1),1,q,pv,pa,pco,hb,PA_mmHg=pa+aa,
        PA_source='published_Pa_plus_Aa_gradient_reconstructs_original_gas_equation',
        mPAP_mmHg=mpap,measured_shunt_fraction=sh,low_vq_fraction=low,
        population='idiopathic_PH' if i<4 else 'recurrent_pulmonary_embolism',
        table='Table II p1052',doi='10.1172/JCI109542',
        notes='Published individual measurements; subject labels recoded. Raw MIGET distributions and individual P50 unavailable. Low-VQ column includes shunt; do not add twice.')
for i,(cond,q,pv,pa,pco,aa,mpap,wedge,sh) in enumerate([
 ('air_rest',4.9,34,64,42,32,27,7,.039),('sildenafil_rest',5.03,32.9,58,40,39,21,5.5,.038),
 ('air_exercise',9.,26,57,46,36,56,22,.029),('sildenafil_exercise',9.32,25.99,56,44,38,45,17.6,.028)]):
    add('Blanco2010',cond,20,q,pv,pa,pco,15.,PA_mmHg=pa+aa,Hb_status='assumed_15_sensitivity_12_18',
        PA_source='published_Pa_plus_Aa_gradient',baseline=int(i==0),mPAP_mmHg=mpap,
        PAWP_mmHg=wedge,measured_shunt_fraction=sh,population='COPD_PH',
        table='Tables 2 and 3 pp272-273',doi='10.1164/rccm.200907-0988OC',
        notes='Post-treatment values reconstructed as reported baseline mean plus reported mean change. Hb absent from tables. External disease-domain stress test, not PAH validation.')
for i,(cond,q,pv,pa,pco,aa,mpap,wedge,sh) in enumerate([
 ('air_rest',5.5,37,77,39,19,20,7,.012),('NO_rest',5.5,37,75,37,21,16,5,.011),
 ('air_exercise',12.,25,59,37,45,40,9,.016),('NO_exercise',12.,25,58,37,44,28,4,.014)]):
    add('Blanco2011',cond,7,q,pv,pa,pco,15.,PA_mmHg=pa+aa,Hb_status='assumed_15_sensitivity_12_18',
        PA_source='published_Pa_plus_Aa_gradient',baseline=int(i==0),mPAP_mmHg=mpap,
        PAWP_mmHg=wedge,measured_shunt_fraction=sh,population='IPF_not_all_resting_PH',
        table='Tables 3 and 4; author thesis pp92 and 94 (PDF 1-based)',doi='10.1152/japplphysiol.01104.2010',
        notes='Author-hosted typeset proof; Hb unavailable. MIGET predicted-minus-measured Pa: air rest -0.4, air exercise 19.6, NO rest 6.0, NO exercise 21.3 mmHg; these are model-derived, not measured PA.')
with (DATA/'observations.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
support=[
 dict(study='Stubbs2021',doi='10.1183/23120541.00449-2021',n=16,conditions=2,
      values='Rest/exercise Q 3.9/6.3 L/min; mPAP 50/77 mmHg; SpO2 94/89%; SvO2 62/23%; PvO2 4.3/2.4 kPa',
      role='Physiological support only',reason='Missing arterial PO2, Hb, PaCO2/PA. SpO2 is not substituted for arterial blood gas.',table='Table 1'),
 dict(study='Beretta2019',doi='10.1038/s41598-019-52679-4',n=18,conditions=4,
      values='Sea-level rest/work and altitude rest/work SpO2 99/98/88/81%; PA 104/103/52/55 mmHg; Hb 14 g/dL',
      role='Altitude boundary support only',reason='No measured mixed-venous input; exercise output partly assumed; table does not supply paired complete model inputs.',table='Table 1 and methods')]
with (DATA/'supporting_observations.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(support[0]));w.writeheader();w.writerows(support)
print(f'{len(rows)} quantitative rows in 5 studies; 6 supporting conditions in 2 studies')
