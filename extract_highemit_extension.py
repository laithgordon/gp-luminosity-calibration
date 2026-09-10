#!/usr/bin/env python3
"""Extract the eps_y = 20-100 nm extension runs into a CSV matching the form of
`lumi_nominal_vs_calibrated.csv` (one row per dataset/eps_y/seed).

Two datasets are emitted:

  tuned   -- grid + n_m from the conservative calibration formulas evaluated at
             each eps_y:  n_y^req = C_y D_y^q (C_y=34.70, q=0.312) rounded UP to
             a power of two; n_m = max(NM_FLOOR, C_m D_y^(s-1/2) n_x n_y n_z)
             with C_m=2.305e-7, s=3.8, NM_FLOOR=5e4.
  frozen  -- grid + n_m held at the eps_y = 20 nm values of those same formulas.

NOTE: the formulas SATURATE over 20-100 nm -- n_y^req falls only 90.4 -> 69.9
(all rounding up to 128) and the locus n_m collapses to 1.6e3-2.4e4, far under
the floor -- so both datasets resolve to the SAME configuration
(512 x 128 x 64, n_m = 50000) at every emittance. The two dataset labels
therefore point at the same underlying runs; `ny_req_raw` records the differing
provenance (per-eps_y value for `tuned`, the frozen 20 nm value for `frozen`).

Reads lumi_ee straight from the .ref files: these runs postdate the last
`lumi_extracted.csv` scrape.

UNITS: GP++ prints lumi_ee/lumi_fine as the per-bunch-crossing luminosity in
m^-2. As in extract_lumi_data.py, `lumi_ee`/`lumi_fine` here are converted to a
rate in cm^-2 s^-1 with the f_rep and n_b echoed in the same .ref
(x 1e-4 * n_b * f_rep = 1.596 for C3-250); the raw values are kept in
`lumi_ee_m2`/`lumi_fine_m2`.
"""
from __future__ import annotations
import json, re, contextlib, io
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')

ANA  = Path(__file__).resolve().parent
REPO = ANA.parent
GH   = REPO / 'GitHub_Analysis'
OUT  = GH / 'data' / 'lumi_tuned_vs_frozen_highemit.csv'

# D_y from the notebook, so the formula inputs match the calibration exactly.
nb = json.load(open(GH / 'GP_CALIBRATION.ipynb')); ns = {'__name__': '__main__'}
with contextlib.redirect_stdout(io.StringIO()):
    for i in [2, 3, 4, 8]:
        exec(compile(''.join(nb['cells'][i]['source']), f'<c{i}>', 'exec'), ns)
D_y_full = ns['D_y_full']

C_M_CONS, S_CONS = 2.305e-07, 3.8
C_Y, Q           = 34.70, 0.312
NM_FLOOR         = 50_000
NX, NZ           = 512, 64
EPS              = [20.0, 40.0, 60.0, 80.0, 100.0]
SEEDS = [100805, 105469, 130471, 135135, 143455, 145992, 153768, 159302,
         167029, 167815, 527875, 733579, 360177, 286038, 479728]
FIXED = dict(emitt_x=0.9, beta_x=12.0, beta_y=0.12, sigma_z=100.0, offset_y=0.0,
             particles=0.624, integration_method=2, n_t=6,
             cut_x=20.0, cut_y=20.0, cut_z=3.5)

def tuned_params(eps_nm):
    Dy  = float(D_y_full(eps_nm))
    nyr = C_Y * Dy ** Q
    ny  = int(2 ** np.ceil(np.log2(nyr)))
    nm  = int(max(NM_FLOOR, C_M_CONS * Dy ** (S_CONS - 0.5) * NX * ny * NZ))
    return ny, nm, nyr

RE_EE   = re.compile(r'\blumi_ee\b\s*=\s*([0-9.eE+-]+)')
RE_FINE = re.compile(r'\blumi_fine\b\s*=\s*([0-9.eE+-]+)')
RE_FREP = re.compile(r'\bf_rep\s*=\s*([0-9.eE+-]+)')
RE_NB   = re.compile(r'\bn_b\s*=\s*([0-9]+)')

def ref_path(eps_tag, ny, nm, seed):
    return (REPO / 'output_nm' / 'C3_250' /
            f'{NX}_{ny}_{NZ}_{nm}_test_C3_250_particles_0.624_emittx_0.9'
            f'_emitty_{eps_tag}_betax_12.0_betay_0.12_sigmaz_100_offsety_0'
            f'_integration_method_2_nt6_seed_{seed}.ref')

def eps_tag(eps_nm):
    """Filename tag as written by the submit scripts (verbatim from the scan)."""
    return '0.020' if eps_nm == 20.0 else f'{eps_nm/1000.:g}'

rows, missing = [], []
ny20, nm20, nyr20 = tuned_params(20.0)

for eps in EPS:
    ny_t, nm_t, nyr_t = tuned_params(eps)
    for ds, ny, nm, nyr in (('tuned',  ny_t, nm_t, nyr_t),
                            ('frozen', ny20, nm20, nyr20)):
        for sd in SEEDS:
            p = ref_path(eps_tag(eps), ny, nm, sd)
            if not p.exists() or p.stat().st_size == 0:
                missing.append((ds, eps, sd, p.name)); continue
            txt = p.read_text(errors='ignore')
            m_ee, m_fi = RE_EE.search(txt), RE_FINE.search(txt)
            m_fr, m_nb = RE_FREP.search(txt), RE_NB.search(txt)
            if not (m_ee and m_fr and m_nb):
                missing.append((ds, eps, sd, p.name)); continue
            f_rep, n_b = float(m_fr.group(1)), int(m_nb.group(1))
            conv = 1e-4 * n_b * f_rep                     # m^-2/crossing -> cm^-2 s^-1
            ee_m2 = float(m_ee.group(1))
            fi_m2 = float(m_fi.group(1)) if m_fi else np.nan
            rows.append(dict(dataset=ds, eps_y_nm=eps, n_x=NX, n_y=ny, n_z=NZ,
                             n_m=nm, ny_req_raw=round(nyr, 1), seed=sd,
                             lumi_ee=ee_m2 * conv, lumi_fine=fi_m2 * conv,
                             **FIXED, source_dir='output_nm', filename=p.name,
                             lumi_ee_m2=ee_m2, lumi_fine_m2=fi_m2,
                             f_rep=f_rep, n_b=n_b))

df = pd.DataFrame(rows)[
    ['dataset','eps_y_nm','n_x','n_y','n_z','n_m','ny_req_raw','seed','lumi_ee',
     'lumi_fine','emitt_x','beta_x','beta_y','sigma_z','offset_y','particles',
     'integration_method','n_t','cut_x','cut_y','cut_z','source_dir','filename',
     'lumi_ee_m2','lumi_fine_m2','f_rep','n_b']
].sort_values(['dataset','eps_y_nm','seed'])
df.to_csv(OUT, index=False)
print(f'wrote {OUT}  ({len(df)} rows)')
if missing:
    print(f'\nMISSING {len(missing)}:')
    for m in missing[:10]: print('  ', m)

s = (df.groupby(['dataset','eps_y_nm'])
       .agg(seeds=('seed','nunique'), n_y=('n_y','first'), n_m=('n_m','first'),
            ny_req=('ny_req_raw','first'),
            mean_L=('lumi_ee','mean'), std_L=('lumi_ee','std')).reset_index())
s['mean_L'] /= 1e34; s['std_L'] /= 1e34
print('\nmean lumi_ee [1e34 cm-2 s-1]:')
print(s.round(4).to_string(index=False))
