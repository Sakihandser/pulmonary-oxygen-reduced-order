"""Extract published table cells without inventing individual participant data.

Snapshots are kept locally for provenance. Downloads are optional; all further
analysis is offline. 2025 abstract summaries are separate from tabulated means.
"""
from pathlib import Path
from html.parser import HTMLParser
import re
import json
import hashlib
import csv

ROOT = Path(__file__).resolve().parents[1]


class TableParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows, self.row, self.buf = [], None, None

    def handle_starttag(self, tag, attrs):
        if tag == 'tr':
            self.row = []
        if tag in ('td', 'th') and self.row is not None:
            self.buf = ''

    def handle_data(self, data):
        if self.buf is not None:
            self.buf += data

    def handle_endtag(self, tag):
        if tag in ('td', 'th') and self.buf is not None:
            self.row.append(' '.join(self.buf.split()))
            self.buf = None
        if tag == 'tr' and self.row is not None:
            self.rows.append(self.row)
            self.row = None


def numeric(cell):
    values = re.findall(r'\d+(?:\.\d+)?', cell)
    if '±' not in cell or len(values) != 2:
        raise ValueError(f'Expected mean and spread: {cell}')
    return tuple(map(float, values))


def extract_2021():
    parser = TableParser()
    parser.feed((ROOT / 'sources/Fermoyle_2021_public.html').read_text(encoding='utf-8'))
    fields = {'Cardiac output, L/min': ('q', 'L/min'),
              'PAPm, mm Hg': ('pap', 'mmHg'), 'PCWPm, mm Hg': ('pw', 'mmHg'),
              'DLCO, mL/min per mm Hg': ('dlco', 'mL/min/mmHg'),
              'DLNO, mL/min per mm Hg': ('dlno', 'mL/min/mmHg'),
              'Dm, mL/min per mm Hg': ('dm', 'mL/min/mmHg'),
              'Vc, mL': ('vc', 'mL'), 'PaO2': ('pao2', 'mmHg'),
              'PaCO2': ('paco2', 'mmHg'), 'SaO2': ('sao2', 'percent'),
              'AaDO2, mm Hg': ('aado2', 'mmHg'), 'RER': ('rer', 'dimensionless'),
              'PVR, mm Hg/L per min': ('pvr', 'WU')}
    records = [dict(study='Fermoyle2021', group=g, stage=s, cohort_n=10 if g=='control' else 20,
                    source_doi='10.1161/JAHA.120.019950', source_location='Table 2',
                    summary_type='mean', spread_type='SEM', q_method='direct_Fick')
               for s in ('rest', '20W', 'peak') for g in ('control', 'HFpEF')]
    ledger=[]
    for row in parser.rows:
        if not row or row[0] not in fields or len(row) < 7:
            continue
        name,unit=fields[row[0]]
        for rec,cell in zip(records,row[1:7]):
            value,spread=numeric(cell)
            rec[name],rec[name+'_spread']=value,spread
            ledger.append(dict(study=rec['study'],group=rec['group'],stage=rec['stage'],
                               variable=name,unit=unit,cell=cell,value=value,spread=spread,
                               spread_type='SEM',location='Table 2'))
    assert all(all(k in r for k in ('q','pap','pw','dm','vc','dlco')) for r in records)
    return records,ledger


def extract_2016():
    parser=TableParser()
    parser.feed((ROOT/'sources/Olson_2016_public.html').read_text(encoding='utf-8'))
    mapping={'Q (L/min)':('q','L/min'), 'DLCO (mL/mmHg/min)':('dlco','mL/min/mmHg'),
             'DM (mL/mmHg/min)':('dm','mL/min/mmHg'), 'VC (mL)':('vc','mL')}
    records=[];ledger=[];stage=None
    for row in parser.rows:
        if not row:continue
        if row[0] in ['Rest','Matched Workload (20 W)','Peak Exercise']:
            stage={'Rest':'rest','Matched Workload (20 W)':'20W','Peak Exercise':'peak'}[row[0]]
            records.extend([dict(study='Olson2016',group=g,stage=stage,cohort_n=26 if g=='control' else 20,
                                 source_doi='10.1016/j.jchf.2016.03.001',source_location='Table 3',
                                 summary_type='mean',spread_type='SD',q_method='echo_SV_times_HR')
                            for g in ('control','HFpEF')])
        if stage and row[0] in mapping:
            name,unit=mapping[row[0]]
            for rec,cell in zip(records[-2:],row[1:3]):
                val,spread=numeric(cell);rec[name],rec[name+'_spread']=val,spread
                ledger.append(dict(study=rec['study'],group=rec['group'],stage=stage,
                                   variable=name,unit=unit,cell=cell,value=val,spread=spread,
                                   spread_type='SD',location='Table 3'))
    assert len(records)==6 and all('dm' in r for r in records)
    return records,ledger


def write_csv(path,records):
    fields=list(dict.fromkeys(k for r in records for k in r))
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(records)


def main():
    data=ROOT/'data';data.mkdir(exist_ok=True)
    r1,l1=extract_2021();r2,l2=extract_2016()
    write_csv(data/'published_groups.csv',r1+r2)
    write_csv(data/'source_cells.csv',l1+l2)
    # DLCO, Q, RVSP have different sample sizes. Do not manufacture paired rows.
    # Vc/Dm figures use medians/IQR; abstract numbers differ from plotted medians.
    rows=[]
    for variable,unit,con,cuff,sc,su,n,location,kind in [
        ('q','L/min',12.1,13.0,1.9,2.3,14,'Table 2','mean_SD'),
        ('rvsp','mmHg',38.9,44.7,5.5,6.1,6,'Table 2','mean_SD'),
        ('dlco','mL/min/mmHg',38,41,6,6,15,'Abstract','mean_SD'),
        ('dm','mL/min/mmHg',90,138,39,55,15,'Abstract','abstract_type_unresolved'),
        ('vc','mL',98,85,27,18,15,'Abstract; units checked against Figure 1','abstract_type_unresolved')]:
        for stage,val,spread in [('CON',con,sc),('CUFF',cuff,su)]:
            rows.append(dict(study='DSouza2025',stage=stage,variable=variable,
                             value=val,spread=spread,n=n,unit=unit,summary_type=kind,
                             source_doi='10.1152/ajplung.00358.2024',source_location=location))
    write_csv(data/'pressure_challenge_summaries.csv',rows)
    manifest=[]
    for p in sorted((ROOT/'sources').iterdir()):
        if p.is_file():manifest.append(dict(path=p.name,bytes=p.stat().st_size,
                    sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    (data/'source_snapshot_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('Extracted 12 group-stage records and 10 pressure-challenge summaries; no individual data.')


if __name__=='__main__':main()
