#!/usr/bin/env python3
"""Per-seed detector-reach counts behind the luminosity-vs-BIB trade-off figure.

Writes data/lumi_bib_per_seed_cache.csv, read by the `lumi_bib_tradeoff_PRL`
cell of GP_CALIBRATION.ipynb. For every run on the two curves of that figure
(eps_x = 900 nm with eps_y swept, and eps_y = 20 nm with eps_x swept, all on
the nominal 512x512x25 grid at n_m = 1e5) it opens the GUINEA-PIG++ pairs.dat
dump written beside the .ref and counts the pair particles that reach the
detector (SiD-o2-v04: B = 5 T, r_det = 14 mm, z_max = 76 mm), via
pairs_reachability.count_reaching_particles.

Run selection is the one the figure uses, applied to data/lumi_extracted.csv:
offset_y = 0, FFT solver, n_t = 6, at most 15 seeds per (eps_x, eps_y) with the
lowest seed numbers kept, and source_dir 'output' preferred when a run exists in
more than one output tree. The selection and counting are those of the notebook that built the cache; the
raw-output paths are configurable.

    python3 bib_reach_per_seed.py [--workers N]
"""
from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import DATA_DIR, raw_dir  # noqa: E402
from pairs_reachability import count_reaching_particles  # noqa: E402

MAX_SEEDS = 15
CURVE_CELLS = set(
    [(900.0, ey) for ey in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 20.0)]
    + [(ex, 20.0) for ex in (100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0)]
)
SRC_RANK = {"output": 0, "output_nm": 1, "output_no_pairs": 2}
COLUMNS = ["filename", "source_dir", "seed", "emitx_nm", "emity_nm", "n_reach", "n_total"]


def select_runs(df: pd.DataFrame) -> pd.DataFrame:
    """The runs the trade-off figure averages over (notebook cell 5)."""
    mask = (
        (df["n_x"] == 512) & (df["n_y"] == 512) & (df["n_z"] == 25)
        & (df["offset_y"] == 0)
        & (df["particles"] == 0.624)
        & (df["beta_x"] == 12.0) & (df["beta_y"] == 0.12)
        & (df["sigma_z"] == 100)
        & ((df["integration_method"] == 2) | df["integration_method"].isna())
        & ((df["n_t"] == 6) | df["n_t"].isna())
        & ((df["n_m"] == 100000) | df["n_m"].isna())
        & (df["emitt_x"].isin([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]))
    )
    sub = df[mask].copy()
    sub["emitx_nm"] = sub["emitt_x"] * 1000
    sub["emity_nm"] = sub["emitt_y"] * 1000
    sub = sub[sub["emity_nm"] <= 20]
    sub["src_rank"] = sub["source_dir"].map(SRC_RANK).fillna(3)
    sub = sub.sort_values(["emitx_nm", "emity_nm", "seed", "src_rank"])
    sub = sub.drop_duplicates(subset=["emitx_nm", "emity_nm", "seed"], keep="first")
    sub = (sub.sort_values(["emitx_nm", "emity_nm", "seed"])
              .groupby(["emitx_nm", "emity_nm"], as_index=False).head(MAX_SEEDS))
    on_curve = sub.apply(lambda r: (float(r["emitx_nm"]), float(r["emity_nm"])) in CURVE_CELLS, axis=1)
    return sub[on_curve]


def pair_path(source_dir: str, ref_filename: str):
    """.ref filename -> the pairs.dat dump GP++ writes beside it."""
    return raw_dir(source_dir) / (ref_filename
                                  .replace("test_C3_250_particles_", "testC3_250_pairs_particles_")
                                  .replace(".ref", ".dat"))


def count_one(rec: dict) -> dict | None:
    p = pair_path(rec["source_dir"], rec["filename"])
    if not (p.exists() and p.stat().st_size > 0):
        return None
    n_reach, n_total, _ = count_reaching_particles(str(p), return_mask=True)
    return dict(filename=rec["filename"], source_dir=rec["source_dir"], seed=int(rec["seed"]),
                emitx_nm=float(rec["emitx_nm"]), emity_nm=float(rec["emity_nm"]),
                n_reach=int(n_reach), n_total=int(n_total))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lumi", default=str(DATA_DIR / "lumi_extracted.csv"))
    ap.add_argument("--out", default=str(DATA_DIR / "lumi_bib_per_seed_cache.csv"))
    ap.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    a = ap.parse_args()

    runs = select_runs(pd.read_csv(a.lumi))
    recs = runs[["filename", "source_dir", "seed", "emitx_nm", "emity_nm"]].to_dict("records")
    print(f"{len(recs)} runs on the trade-off curves; counting detector reach "
          f"with {a.workers} worker(s)")
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        rows = [r for r in ex.map(count_one, recs, chunksize=4) if r is not None]
    missing = len(recs) - len(rows)
    out = pd.DataFrame(rows, columns=COLUMNS).sort_values(["emitx_nm", "emity_nm", "seed"])
    out.to_csv(a.out, index=False)
    print(f"wrote {a.out} ({len(out)} rows; {missing} runs had no pair dump)")


if __name__ == "__main__":
    main()
