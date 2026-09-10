#!/usr/bin/env python3
"""BIB pair yield and detector reachability: nominal vs calibrated GP++ grids.

For each (dataset, eps_y, seed) it reads the GUINEA-PIG++ `pairs.dat` dump and
reports, per production channel {BW, BH, LL} and for all channels combined:

    count, median p_T, and the percentage that reach the detector.

All reachability physics is delegated to `reachability_analysis.compute_boundary`
via `pairs_reachability.reaches_detector` (SiD-o2-v04: B = 5 T, r_det = 14 mm,
z_max = 76 mm) -- this script only selects runs and aggregates.

    nominal    : n_x=512, n_y=512, n_z=25, n_m=1e5   (pre-calibration baseline)
    calibrated : n_x=512, n_z=64, n_y from the n_y calibration rounded up to a
                 simulated value, n_m from the conservative locus.

Usage:  python bib_nominal_vs_calibrated.py [--seed-cap N] [--eps 1 8 20]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ANA = Path(__file__).resolve().parent
REPO = ANA.parent
sys.path.insert(0, str(ANA))
from pairs_reachability import load_pairs_dat, reaches_detector  # noqa: E402

CHAN = {0: "BW", 1: "BH", 2: "LL"}
GEO = dict(mag_field=5.0, detector_radius=14.0, z_max=76.0, charge=0.3)

# Full run identity. Matching on (n_x,n_y,n_z,n_m,eps_y) alone is NOT enough:
# the BIB campaigns swept eps_x over the same grids, so a looser match silently
# picks up eps_x=0.1 runs.
FIXED = dict(emitt_x=0.9, beta_x=12.0, beta_y=0.12, sigma_z=100.0, offset_y=0.0,
             particles=0.624, integration_method=2, n_t=6,
             cut_x=20.0, cut_y=20.0, cut_z=3.5)

CSV_DIR = REPO / "GitHub_Analysis" / "data"


def pair_path(source_dir: str, ref_filename: str) -> Path:
    """`.ref` filename -> the pair-dump path GP++ writes beside it."""
    return REPO / source_dir / "C3_250" / (
        ref_filename.replace("test_C3_250_particles_",
                             "testC3_250_pairs_particles_").replace(".ref", ".dat"))


def find_dumps(full: pd.DataFrame, eps_nm, nx, ny, nz, nm, seeds) -> dict[int, Path]:
    """seed -> existing pair-dump path, searching every source_dir."""
    c = full[(full.n_x == nx) & (full.n_y == ny) & (full.n_z == nz) & (full.n_m == nm)
             & np.isclose(full.emitt_y, eps_nm / 1000.0) & full.seed.isin(seeds)
             & (full.source_dir != "output_no_pairs")]
    for k, v in FIXED.items():
        c = c[np.isclose(pd.to_numeric(c[k], errors="coerce"), v)]
    out: dict[int, Path] = {}
    for r in c.itertuples():
        p = pair_path(r.source_dir, r.filename)
        if p.exists() and p.stat().st_size > 0:
            out.setdefault(int(r.seed), p)
    # Fresh runs are absent from lumi_extracted.csv until it is re-scraped, so
    # also probe the deterministic path the pair-dump submit scripts write.
    for sd in sorted(seeds):
        if sd in out:
            continue
        # Glob rather than format the emitt_y field: the submit scripts write it
        # verbatim from the scan file, so eps_y = 20 nm lands as "0.020", which
        # an f-string "{0.02:g}" would never match.
        pat = (f"{int(nx)}_{int(ny)}_{int(nz)}_{int(nm)}_testC3_250_pairs_particles"
               f"_0.624_emittx_0.9_emitty_*_betax_12.0_betay_0.12"
               f"_sigmaz_100_offsety_0_integration_method_2_nt6_seed_{sd}.dat")
        for src in ("output_nm", "output"):
            for q in (REPO / src / "C3_250").glob(pat):
                m = re.search(r"_emitty_([0-9.eE+-]+)_betax", q.name)
                if not m or not np.isclose(float(m.group(1)), eps_nm / 1000.0):
                    continue
                if q.stat().st_size > 0:
                    out[sd] = q
                    break
            if sd in out:
                break
    return out


def analyse(path: Path) -> list[dict]:
    """Per-channel count / median p_T / reach fraction for one pair dump."""
    d = load_pairs_dat(str(path))
    if d["pt"].size == 0:
        return []
    mask = reaches_detector(d["pt"], d["theta"], **GEO)

    def row(name: str, sel: np.ndarray) -> dict:
        n, nr = int(sel.sum()), int((sel & mask).sum())
        return dict(
            channel=name, n_pairs=n, n_reach=nr,
            reach_pct=100.0 * nr / n if n else np.nan,
            median_pt_MeV=1e3 * float(np.median(d["pt"][sel])) if n else np.nan,
            median_pt_reach_MeV=(1e3 * float(np.median(d["pt"][sel & mask]))
                                 if nr else np.nan))

    rows = [row(name, d["process"] == code) for code, name in CHAN.items()]
    rows.append(row("ALL", np.ones(d["pt"].size, dtype=bool)))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed-cap", type=int, default=1)
    ap.add_argument("--eps", type=float, nargs="+", default=[1.0, 8.0, 20.0])
    ap.add_argument("--out", default=str(CSV_DIR / "bib_nominal_vs_calibrated.csv"))
    a = ap.parse_args()

    full = pd.read_csv(CSV_DIR / "lumi_extracted.csv")
    for c in ["n_x", "n_y", "n_z", "n_m", "seed", "emitt_y"]:
        full[c] = pd.to_numeric(full[c], errors="coerce")
    full["n_m"] = full["n_m"].fillna(100000)
    for c, v in [("n_t", 6), ("integration_method", 2), ("cut_x", 20.0),
                 ("cut_y", 20.0), ("cut_z", 3.5), ("offset_y", 0.0), ("beta_y", 0.12)]:
        full[c] = pd.to_numeric(full[c], errors="coerce").fillna(v)
    tgt = pd.read_csv(CSV_DIR / "lumi_nominal_vs_calibrated.csv")

    recs = []
    for eps in a.eps:
        for ds in ["nominal", "calibrated"]:
            blk = tgt[(tgt.dataset == ds) & (tgt.eps_y_nm == eps)]
            if not len(blk):
                print(f"  {ds:<10} {eps:>5g} nm : no config row")
                continue
            r0 = blk.iloc[0]
            dumps = find_dumps(full, eps, r0.n_x, r0.n_y, r0.n_z, r0.n_m, set(blk.seed))
            tag = (f"n_x={int(r0.n_x)} n_y={int(r0.n_y)} n_z={int(r0.n_z)} "
                   f"n_m={int(r0.n_m)}")
            if not dumps:
                print(f"  {ds:<10} {eps:>5g} nm : 0 pair dumps  ({tag})  -- PENDING")
                continue
            use = sorted(dumps)[: a.seed_cap]
            print(f"  {ds:<10} {eps:>5g} nm : {len(dumps)} dump(s), analysing "
                  f"{len(use)}  ({tag})")
            for sd in use:
                for row in analyse(dumps[sd]):
                    recs.append(dict(dataset=ds, eps_y_nm=eps, seed=sd,
                                     n_x=int(r0.n_x), n_y=int(r0.n_y),
                                     n_z=int(r0.n_z), n_m=int(r0.n_m), **row))

    if not recs:
        print("no data")
        return
    df = pd.DataFrame(recs).sort_values(["eps_y_nm", "dataset", "seed", "channel"])
    df.to_csv(a.out, index=False)
    print(f"\nwrote {a.out}  ({len(df)} rows)")


if __name__ == "__main__":
    main()
