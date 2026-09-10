#!/usr/bin/env python3
"""Rebuild data/lumi_nominal_vs_calibrated.csv from the root extract.

The file lists, per (dataset, eps_y, seed), the .ref run that defines the
"nominal" (512x512x25, n_m=1e5) and "calibrated" (n_y^req / n_m^req locus)
points. Its luminosity columns were historically copied from the .ref files,
i.e. per-bunch-crossing m^-2. This script keeps the row set and every non-lumi
column exactly as they are, and re-sources `lumi_ee` / `lumi_fine` from
data/lumi_extracted.csv (which converts to cm^-2 s^-1 per row from the .ref's
own f_rep and n_b, see extract_lumi_data.py) plus the raw `*_m2` columns.

    python3 rebuild_lumi_nominal_vs_calibrated.py
"""
import sys
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
TGT  = HERE / 'data' / 'lumi_nominal_vs_calibrated.csv'
ROOT = HERE / 'data' / 'lumi_extracted.csv'
NEW  = ['lumi_ee', 'lumi_fine', 'lumi_ee_m2', 'lumi_fine_m2', 'f_rep', 'n_b']

tgt  = pd.read_csv(TGT)
root = pd.read_csv(ROOT, usecols=['source_dir', 'filename'] + NEW)
if 'lumi_ee_m2' not in root.columns:
    sys.exit('lumi_extracted.csv predates the unit fix - run extract_lumi_data.py first')

keep = [c for c in tgt.columns if c not in NEW]
out  = tgt[keep].merge(root, on=['source_dir', 'filename'], how='left', validate='one_to_one')
miss = out['lumi_ee_m2'].isna().sum()
if miss:
    sys.exit(f'{miss} rows of {TGT.name} have no match in {ROOT.name}; aborting')

# original column order, new raw/rate-factor columns appended
order = [c for c in tgt.columns if c in out.columns] + [c for c in NEW if c not in tgt.columns]
out = out[order].sort_values(['dataset', 'eps_y_nm', 'seed']).reset_index(drop=True)
out.to_csv(TGT, index=False)

s = (out.groupby(['dataset', 'eps_y_nm'])
        .agg(n=('seed', 'size'), L=('lumi_ee', 'mean'), L_m2=('lumi_ee_m2', 'mean'))
        .reset_index())
s['L [1e34 cm-2 s-1]'] = s.pop('L') / 1e34
s['L_m2 [1e34 m-2/xing]'] = s.pop('L_m2') / 1e34
print(f'wrote {TGT}  ({len(out)} rows)\n')
print(s.round(4).to_string(index=False))
