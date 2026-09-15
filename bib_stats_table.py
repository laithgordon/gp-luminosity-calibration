#!/usr/bin/env python3
"""GP++ beam-induced-background statistics at eps_y = 2, 4, 8, 20 nm (tab:bib_yields).

Per production channel {BW, BH, LL} and for all channels combined, per seed:
counts, detector-reaching counts and fraction, and p_T / energy / polar-angle
statistics. Reachability physics is delegated to
reachability_analysis.compute_boundary (SiD-o2-v04: B=5 T, r_det=14 mm,
z_max=76 mm); this script only selects runs and aggregates.

Datasets:
  nominal  - pre-calibration grid 512x512x25, n_m = 1e5
  tuned    - pre-campaign tuned grid, n_x=512, n_z=64: n_y = 34.70*D_y^0.312
             rounded up to a power of two, n_m at or above the n_m locus.
             This is NOT the published conservative locus (which has twice
             the n_y); see README.md, 'Which rule set each configuration'.
"""
from __future__ import annotations
import argparse, glob, hashlib, os, re, sys
import numpy as np, pandas as pd

ANA = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ANA)
from paths import DATA_DIR, RAW_ROOT  # noqa: E402  (set GP_RAW_ROOT, see paths.py)
from pairs_reachability import load_pairs_dat, reaches_detector  # noqa: E402

CHAN = {0: 'BW', 1: 'BH', 2: 'LL'}
GEO = dict(mag_field=5.0, detector_radius=14.0, z_max=76.0, charge=0.3)
CACHE = str(DATA_DIR / 'bib_stats_cache.csv')

# (dataset, eps_y_nm, n_x, n_y, n_z, n_m)
CONFIGS = [
    # nominal: pre-calibration grid, identical at every eps_y
    ('nominal', 2.0,  512, 512, 25, 100000),
    ('nominal', 4.0,  512, 512, 25, 100000),
    ('nominal', 8.0,  512, 512, 25, 100000),
    ('nominal', 20.0, 512, 512, 25, 100000),
    # tuned: pre-campaign tuned grid (not the conservative locus), n_y and n_m per eps_y
    ('tuned',   2.0,  512, 256, 64, 2235661),
    ('tuned',   4.0,  512, 128, 64, 353910),
    ('tuned',   8.0,  512, 128, 64, 250000),
    ('tuned',   20.0, 512, 128, 64, 80000),
    # eps_y = 1 nm tuned (n_y=256, n_m=7.05e6) is omitted: every seed exceeded
    # the 10-day wall without completing (see README.md, "Which rule set each configuration").
]

def _sha256(path, chunk=1 << 22):
    """Content digest of a file, read in chunks (pair dumps reach ~50 MB)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()

def find_dumps(eps_nm, nx, ny, nz, nm, cap):
    """seed -> pair-dump path, globbing both output dirs and parsing emitt_y
    (the submit scripts write it verbatim, so 20 nm appears as both 0.02/0.020)."""
    # Collect EVERY match per seed rather than first-wins. Both globs pin
    # `_offsety_0_`, so beam-offset scan files (offset_y = 0.313 sigma_y) are not
    # matched, and an ambiguity that survives that is raised rather than resolved
    # arbitrarily.
    cand: dict[int, set] = {}
    for src in ('output_nm', 'output'):
        d = os.path.join(str(RAW_ROOT), src, 'C3_250')
        for pat in (f'{nx}_{ny}_{nz}_{nm}_testC3_250_pairs_*emittx_0.9_emitty_*'
                    f'_offsety_0_*.dat',
                    f'{nx}_{ny}_{nz}_testC3_250_pairs_*emittx_0.9_emitty_*'
                    f'_offsety_0_*.dat'):
            for p in glob.glob(os.path.join(d, pat)):
                b = os.path.basename(p)
                m = re.search(r'_emitty_([0-9.eE+-]+)_', b)
                s = re.search(r'_seed_(\d+)\.dat$', b)
                if not m or not s:
                    continue
                if not np.isclose(float(m.group(1)), eps_nm / 1000.0):
                    continue
                # belt and braces: the glob pins offset 0, re-assert it here so a
                # future pattern edit cannot silently reintroduce the defect.
                o = re.search(r'_offsety_([0-9.eE+-]+)_', b)
                if not o or float(o.group(1)) != 0.0:
                    continue
                if os.path.getsize(p) == 0:
                    continue
                cand.setdefault(int(s.group(1)), set()).add(p)

    out = {}
    for sd, paths in cand.items():
        if len(paths) > 1:
            # GP++ wrote some runs under two filenames (with and without the n_m
            # prefix, in output/ and output_nm/). Two candidates are accepted as
            # the same run written twice ONLY if they are byte-identical; equal
            # size is not taken as evidence. Any difference in content means two
            # different runs claim this (config, seed), and that is refused.
            digests = {_sha256(p) for p in paths}
            if len(digests) > 1:
                raise RuntimeError(
                    'ambiguous pair dump for seed %d at eps_y=%g nm, grid '
                    '(%s,%s,%s,%s): %d candidates with differing content:\n  %s'
                    % (sd, eps_nm, nx, ny, nz, nm, len(paths),
                       '\n  '.join(sorted(paths))))
        # One match, or byte-identical copies of a single dump: every candidate
        # gives identical statistics, so which path is recorded cannot change a
        # result. sorted() only makes the recorded path deterministic.
        out[sd] = sorted(paths)[0]
    return {k: out[k] for k in sorted(out)[:cap]}

def stats(path):
    d = load_pairs_dat(str(path))
    if d['pt'].size == 0:
        return []
    mask = reaches_detector(d['pt'], d['theta'], **GEO)
    th = np.minimum(d['theta'], np.pi - d['theta'])       # folded polar angle
    rows = []
    def row(name, sel):
        n, nr = int(sel.sum()), int((sel & mask).sum())
        pt, E = d['pt'][sel], d['E'][sel]
        ptr = d['pt'][sel & mask]
        q = lambda a, p: float(np.percentile(a, p)) if a.size else np.nan
        return dict(channel=name, n_pairs=n, n_reach=nr,
                    frac_of_total_pct=100.0 * n / d['pt'].size,
                    reach_frac_of_total_pct=(100.0 * nr / int(mask.sum())
                                             if mask.sum() else np.nan),
                    reach_pct=100.0 * nr / n if n else np.nan,
                    pt_median_MeV=1e3 * q(pt, 50), pt_mean_MeV=1e3 * float(pt.mean()) if n else np.nan,
                    pt_p90_MeV=1e3 * q(pt, 90), pt_p99_MeV=1e3 * q(pt, 99),
                    pt_max_MeV=1e3 * float(pt.max()) if n else np.nan,
                    pt_median_reach_MeV=1e3 * q(ptr, 50),
                    E_median_GeV=q(E, 50), E_mean_GeV=float(E.mean()) if n else np.nan,
                    theta_median_mrad=1e3 * q(th[sel], 50),
                    frac_pt_gt_20MeV=100.0 * float((pt > 0.020).mean()) if n else np.nan)
    for code, nm_ in CHAN.items():
        rows.append(row(nm_, d['process'] == code))
    rows.append(row('ALL', np.ones(d['pt'].size, dtype=bool)))
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed-cap', type=int, default=10)
    ap.add_argument('--refresh', action='store_true')
    a = ap.parse_args()

    cache = pd.read_csv(CACHE) if (os.path.exists(CACHE) and not a.refresh) else pd.DataFrame()
    done = set() if cache.empty else set(zip(cache.dataset, cache.eps_y_nm, cache.seed))
    recs = []
    for ds, eps, nx, ny, nz, nm in CONFIGS:
        dumps = find_dumps(eps, nx, ny, nz, nm, a.seed_cap)
        if not dumps:
            print(f'  {ds:<8} {eps:>5g} nm : NO pair dumps '
                  f'(n_x={nx} n_y={ny} n_z={nz} n_m={nm})')
            continue
        todo = [s for s in dumps if (ds, eps, s) not in done]
        print(f'  {ds:<8} {eps:>5g} nm : {len(dumps)} dumps, {len(todo)} to process')
        for s in todo:
            for r in stats(dumps[s]):
                recs.append(dict(dataset=ds, eps_y_nm=eps, seed=s, n_x=nx, n_y=ny,
                                 n_z=nz, n_m=nm, **r))
    df = pd.concat([cache, pd.DataFrame(recs)], ignore_index=True) if recs else cache
    if df.empty:
        print('no data'); return
    df.to_csv(CACHE, index=False)
    print(f'\ncached {len(df)} rows -> {CACHE}')

if __name__ == '__main__':
    main()
