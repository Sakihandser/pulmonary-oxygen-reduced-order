"""Read the attributed aggregate transcriptions supplied with this study."""
from pathlib import Path
import csv
ROOT=Path(__file__).resolve().parents[1]
def write_csv(path,rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
def main():
    for name,expected in [('published_groups.csv',12),('source_cells.csv',102)]:
        with (ROOT/'data'/name).open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
        assert len(rows)==expected,(name,len(rows))
    assert (ROOT/'data/pressure_challenge_summaries.csv').exists()
