#!/usr/bin/env python3
"""Render the GP++ BIB statistics table (mean +- sample STD over seeds)."""
import os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import DATA_DIR  # noqa: E402

D = str(DATA_DIR)
df = pd.read_csv(os.path.join(D, 'bib_stats_cache.csv'))

MET = ['n_pairs','n_reach','reach_pct','frac_of_total_pct','reach_frac_of_total_pct',
       'pt_median_MeV','pt_mean_MeV','pt_p90_MeV','pt_p99_MeV','pt_max_MeV',
       'pt_median_reach_MeV','E_median_GeV','theta_median_mrad','frac_pt_gt_20MeV']
g = (df.groupby(['eps_y_nm','dataset','channel'])
       .agg(seeds=('seed','nunique'), n_x=('n_x','first'), n_y=('n_y','first'),
            n_z=('n_z','first'), n_m=('n_m','first'),
            **{m:(m,'mean') for m in MET},
            **{m+'_sd':(m,'std') for m in MET}).reset_index())
g.to_csv(os.path.join(D,'bib_stats_summary.csv'), index=False)

def c(r,k,fmt,sd=True):
    if not len(r) or not np.isfinite(r[k].iloc[0]): return '-'
    v=r[k].iloc[0]; s=r[k+'_sd'].iloc[0]
    return f'{fmt.format(v)} ± {fmt.format(s)}' if (sd and np.isfinite(s)) else fmt.format(v)

L = ['# GP++ beam-induced background — pair yield, composition and detector reach','',
     'Per bunch crossing, mean ± sample STD over 10 canonical seeds. Reach test',
     'delegated to `reachability_analysis.compute_boundary` (SiD-o2-v04: B = 5 T,',
     'r_det = 14 mm, z_max = 76 mm); a pair particle "reaches" iff its helix from the',
     'IP crosses the barrel wall before exiting the endcap.','',
     'Production channels: **BW** Breit–Wheeler (real γγ), **BH** Bethe–Heitler',
     '(virtual–real), **LL** Landau–Lifshitz (virtual–virtual). `share` is that',
     "channel's fraction of all pairs produced; `reach share` its fraction of all",
     'pairs that reach the detector.','',
     '- **nominal** — pre-calibration grid 512×512×25, n_m = 10⁵ at every ε_y',
     '- **tuned** — conservative-formula grid, n_x = 512, n_z = 64, n_y and n_m per ε_y','']

# ---- headline table: totals ----
L += ['## 1. Totals (all channels)','',
      '| ε_y [nm] | set | n_y | n_m | total pairs | reach | reach % | med p_T [MeV] | med p_T of reaching [MeV] | med θ [mrad] |',
      '|---|---|---|---|---|---|---|---|---|---|']
for eps in sorted(g.eps_y_nm.unique()):
    for ds in ['nominal','tuned']:
        r = g[(g.eps_y_nm==eps)&(g.dataset==ds)&(g.channel=='ALL')]
        if not len(r): continue
        cfg = df[(df.eps_y_nm==eps)&(df.dataset==ds)].iloc[0]
        L.append(f'| {eps:g} | {ds} | {cfg.n_y} | {cfg.n_m:,} | '
                 f'{c(r,"n_pairs","{:,.0f}")} | {c(r,"n_reach","{:.1f}")} | '
                 f'{c(r,"reach_pct","{:.3f}")} | {c(r,"pt_median_MeV","{:.3f}")} | '
                 f'{c(r,"pt_median_reach_MeV","{:.1f}")} | '
                 f'{c(r,"theta_median_mrad","{:.1f}",False)} |')
L.append('')

# ---- per-channel tables ----
L += ['## 2. Per production channel','']
for eps in sorted(g.eps_y_nm.unique()):
    L += [f'### ε_y = {eps:g} nm','']
    for ds in ['nominal','tuned']:
        sub = g[(g.eps_y_nm==eps)&(g.dataset==ds)]
        if not len(sub): 
            L += [f'_**{ds}**: no pair dumps available._','']; continue
        cfg = df[(df.eps_y_nm==eps)&(df.dataset==ds)].iloc[0]
        L += [f'**{ds}** — n_x={cfg.n_x}, n_y={cfg.n_y}, n_z={cfg.n_z}, n_m={cfg.n_m:,}','',
              '| channel | pairs | share % | reach | reach share % | reach % | '
              'med p_T [MeV] | mean p_T | p99 p_T | med p_T reaching | med E [GeV] | med θ [mrad] |',
              '|---|---|---|---|---|---|---|---|---|---|---|---|']
        for ch in ['BW','BH','LL','ALL']:
            r = sub[sub.channel==ch]
            if not len(r): continue
            row=[ch, c(r,'n_pairs','{:,.0f}'), c(r,'frac_of_total_pct','{:.2f}',False),
                 c(r,'n_reach','{:.1f}'), c(r,'reach_frac_of_total_pct','{:.1f}',False),
                 c(r,'reach_pct','{:.3f}'), c(r,'pt_median_MeV','{:.3f}'),
                 c(r,'pt_mean_MeV','{:.2f}',False), c(r,'pt_p99_MeV','{:.1f}',False),
                 c(r,'pt_median_reach_MeV','{:.1f}'), c(r,'E_median_GeV','{:.3f}',False),
                 c(r,'theta_median_mrad','{:.1f}',False)]
            if ch=='ALL': row=[f'**{x}**' for x in row]
            L.append('| '+' | '.join(row)+' |')
        L.append('')

out=os.path.join(D,'bib_stats_table.md'); open(out,'w').write('\n'.join(L)+'\n')
print('\n'.join(L[:40])); print(f'\n... wrote {out}')
