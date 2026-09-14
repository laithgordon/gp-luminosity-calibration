#!/usr/bin/env python3
"""Luminosity vs n_m for n_x = 256, 512, 1024, 2048 at eps_y = 20 nm.

Every point is the single common seed 100805, so all curves share the same
initial-particle realisation and differences between them are not seed noise.
Fixed: n_y = 512, n_z = 64, n_t = 6, FFT solver, C3-250 beam.

Converged value per n_x: Richardson extrapolation in n_m through the three
largest points (1e6, 5e6, 1e7), L = L_inf + C * n_m^-p, solved exactly.
A series counts as converged only if its last step (5e6 -> 1e7) moves L by
less than CONV_TOL; unconverged series are plotted but excluded from the
quoted spread.

Repository-only figure (not in the manuscript).

Input is data/L_vs_nm_by_nx_20nm.csv, one luminosity per (n_x, n_m). These runs
were made after the data/lumi_extracted.csv snapshot, so the table is rebuilt
straight from the raw .ref files with --refresh, using the extractor's own
parsing and unit conversion (needs GP_RAW_ROOT, see paths.py).

    python3 nx_nm_convergence.py              # plot from the committed table
    python3 nx_nm_convergence.py --refresh    # rebuild the table from raw .ref first

Writes plots/L_vs_nm_by_nx_20nm.{png,pdf}.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
TABLE = os.path.join(HERE, "data", "L_vs_nm_by_nx_20nm.csv")
SEED = 100805
NM_MIN = 1e5
NX_SERIES = [256, 512, 1024, 2048]
CONV_TOL = 0.001          # relative change over the last n_m doubling
FIXED = dict(emitt_y=0.020, emitt_x=0.9, n_y=512, n_z=64, n_t=6,
             integration_method=2, beta_x=12.0, beta_y=0.12,
             sigma_z=100, offset_y=0.0, particles=0.624)

def table_from_raw() -> pd.DataFrame:
    """Parse every seed-SEED .ref on the (n_y, n_z) = (512, 64) grid exactly as
    extract_lumi_data.py does, and reduce to one luminosity per (n_x, n_m)."""
    from extract_lumi_data import SOURCE_DIRS, build_row, parse_content, parse_filename
    rows = []
    for key, src in SOURCE_DIRS.items():
        for path in sorted(glob.glob(os.path.join(src, f"*_512_64_*_seed_{SEED}.ref"))):
            fn = os.path.basename(path)
            params = parse_filename(fn)
            metrics = parse_content(path) if params is not None else None
            if metrics is not None:
                rows.append(build_row(key, fn, params, metrics))
    df = pd.DataFrame(rows)
    for c in df.columns:
        if c not in ("source_dir", "filename"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
    m = df.lumi_ee.notna() & df.n_m.notna() & (df.seed == SEED)
    for k, v in FIXED.items():
        m &= df[k] == v
    m &= df.n_x.isin(NX_SERIES) & (df.n_m >= NM_MIN)
    t = (df[m].groupby(["n_x", "n_m"]).lumi_ee.mean() / 1e34).reset_index()
    return t.rename(columns={"lumi_ee": "L"}).sort_values(["n_x", "n_m"])


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--refresh", action="store_true", help="rebuild data/L_vs_nm_by_nx_20nm.csv from raw .ref files")
args = ap.parse_args()
if args.refresh or not os.path.exists(TABLE):
    table_from_raw().to_csv(TABLE, index=False)
    print(f"rebuilt {TABLE} from raw .ref files")
g = pd.read_csv(TABLE).sort_values(["n_x", "n_m"])


def richardson3(n, L):
    """Exact L = Linf + C n^-p through three (n, L) points with n1<n2<n3."""
    (n1, n2, n3), (L1, L2, L3) = n, L
    r = (L1 - L2) / (L2 - L3)
    f = lambda p: (n1**-p - n2**-p) / (n2**-p - n3**-p) - r
    p = brentq(f, 1e-3, 5.0)
    C = (L2 - L3) / (n2**-p - n3**-p)
    return L3 - C * n3**-p, p


rows = []
for nx in NX_SERIES:
    s = g[g.n_x == nx].set_index("n_m").L
    top = s.loc[[1e6, 5e6, 1e7]]
    step = abs(top.iloc[2] - top.iloc[1]) / top.iloc[2]
    try:
        Linf, p = richardson3(top.index.values, top.values)
    except ValueError:
        # No power-law solution: the top three points are flat to within noise
        # (a decaying n^-p tail cannot produce equal successive steps). The
        # series has already converged, so L_inf is its largest-n_m value.
        Linf, p = top.iloc[2], np.nan
    rows.append(dict(n_x=nx, L_1e7=top.iloc[2], last_step_pct=100 * step,
                     p=p, L_inf=Linf, converged=step < CONV_TOL))
res = pd.DataFrame(rows)

# ── figure ────────────────────────────────────────────────────────────────────
colors = {256: "#4C72B0", 512: "#000000", 1024: "#DD8452", 2048: "#C44E52"}
marks = {256: "s", 512: "o", 1024: "^", 2048: "D"}
fig, ax = plt.subplots(figsize=(6.4, 4.4))
conv = res[res.converged]
Lbar = conv.L_inf.mean()
ax.axhline(Lbar, color="0.55", lw=1.0, ls=":", zorder=0)
for nx in NX_SERIES:
    s = g[g.n_x == nx]
    ok = bool(res.loc[res.n_x == nx, "converged"].iloc[0])
    ax.plot(s.n_m, s.L, marker=marks[nx], ms=5.5, lw=1.4, color=colors[nx],
            mfc=colors[nx] if ok else "white", mew=1.3,
            label=rf"$n_x = {nx}$" + ("" if ok else " (not converged)"))
ax.set_xscale("log")
ax.set_xlabel(r"$n_m$  (macroparticles per beam)")
ax.set_ylabel(r"$\mathscr{L}$  [$10^{34}\,$cm$^{-2}\,$s$^{-1}$]")
ax.set_title(r"$\varepsilon_y^* = 20$ nm,  $n_y = 512$,  $n_z = 64$,  "
             rf"seed {SEED}", fontsize=10)
ax.legend(frameon=False, fontsize=9, loc="upper right")
ax.grid(alpha=0.25, lw=0.6)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, "plots", f"L_vs_nm_by_nx_20nm.{ext}"), dpi=200)

# ── report ────────────────────────────────────────────────────────────────────
pd.set_option("display.width", 120)
print(g.pivot(index="n_m", columns="n_x", values="L").round(4).to_string())
print("\nper-series convergence (Richardson through 1e6, 5e6, 1e7):")
print(res.round(4).to_string(index=False))
print(f"\nmean L_inf over converged series: {Lbar:.4f}")
print("pairwise differences between converged series:")
c = conv.set_index("n_x")
for i, a in enumerate(c.index):
    for b in c.index[i + 1:]:
        for col in ("L_1e7", "L_inf"):
            d = 100 * (c.loc[a, col] - c.loc[b, col]) / Lbar
            print(f"  n_x {a:>4} vs {b:>4}  [{col:>5}]  {d:+.3f}%")
for _, r in res[~res.converged].iterrows():
    print(f"  n_x {int(r.n_x)} excluded: last step {r.last_step_pct:.2f}%, "
          f"L(1e7) {100*(r.L_1e7-Lbar)/Lbar:+.2f}% and L_inf "
          f"{100*(r.L_inf-Lbar)/Lbar:+.2f}% vs converged mean")
