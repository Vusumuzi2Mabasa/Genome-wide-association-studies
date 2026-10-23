#!/usr/bin/env python3
"""Export compact, synthetic-only evidence from a successful workflow run."""
import argparse,json,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('--out',type=Path,default=Path('results'));a=p.parse_args()
v=json.loads((a.run/'validation.json').read_text())
if v.get('status')!='PASS' or v.get('stage')!='all' or not v.get('inputs_unchanged'):raise ValueError('A complete successful run is required')
a.out.mkdir(parents=True,exist_ok=True)
for name in ['validation.json','versions.json']:
 shutil.copy2(a.run/name,a.out/name)
for folder in ['summary','genotype','reads']:
 dest=a.out/folder;dest.mkdir(exist_ok=True)
 for f in (a.run/folder).glob('*.tsv'):shutil.copy2(f,dest/f.name)
for name in ['assoc.PHENO.glm.logistic.hybrid','pcs.eigenvec','pcs.eigenval','ld.prune.in','hwe_0.hardy','hwe_1.hardy','qc_stats.afreq']:
 shutil.copy2(a.run/'genotype'/name,a.out/'genotype'/name)
print('Exported synthetic evidence to',a.out)
