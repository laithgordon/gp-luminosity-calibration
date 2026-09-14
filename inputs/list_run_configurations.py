#!/usr/bin/env python3
"""Every run configuration behind the tuning ladders, as a table.

The ladder figures and fitted constants are built from data/lumi_extracted.csv,
which records one row per GUINEA-PIG++ run with every varied input parameter
parsed from its filename. This script reduces it to one row per distinct
configuration, with the seeds that were run, and writes
inputs/run_configurations.csv. Any row, combined with an input file in inputs/
and a seed, reproduces that run through inputs/submit_template_S3DF_SLURM.sh.
The columns absent from a filename (n_m, integration_method, n_t, cuts) are
left empty exactly as in the table; the template defaults then apply
(n_m = 1e5, integration_method = 2, n_t = 6, cuts 20/20/3.5).

    python3 inputs/list_run_configurations.py
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
CONFIG = ["n_x", "n_y", "n_z", "n_m", "particles", "emitt_x", "emitt_y", "beta_x", "beta_y",
          "sigma_z", "offset_y", "integration_method", "n_t", "cut_x", "cut_y", "cut_z", "grids"]


def main() -> None:
    df = pd.read_csv(HERE.parent / "data" / "lumi_extracted.csv")
    g = (df.groupby(CONFIG, dropna=False)
           .agg(n_runs=("seed", "size"), n_seeds=("seed", "nunique"),
                source_dirs=("source_dir", lambda s: ";".join(sorted(set(s)))),
                seeds=("seed", lambda s: ";".join(str(int(x)) for x in sorted(set(s)))))
           .reset_index()
           .sort_values(["emitt_y", "emitt_x", "n_x", "n_y", "n_z", "n_m"], na_position="first"))
    out = HERE / "run_configurations.csv"
    g.to_csv(out, index=False)
    print(f"wrote {out}: {len(g)} configurations, {int(g.n_runs.sum())} runs")


if __name__ == "__main__":
    main()
