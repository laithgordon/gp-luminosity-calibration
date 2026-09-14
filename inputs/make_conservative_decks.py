#!/usr/bin/env python3
"""Scan-file rows for the conservative-locus configurations.

This is the rule the conservative block of data/gp_luminosity_for_wx.csv and
tab:ny_recommendations use, with the frozen constants of constants.py:

    n_y = 2^ceil(log2(50 * D_y^0.402))                      rounded UP
    n_m = ceil(2.305e-7 * D_y^3.300 * n_x * n_y * n_z),      n_x = 512, n_z = 64

with n_t = 6, integration_method = 2 (FFT), offset_y = 0 and the C3-250 beam of
inputs/acc_C3_250_nominal.dat. D_y per eps_y is data/D_y_table.json.

The rows are in the scan-file format of inputs/submit_template_S3DF_SLURM.sh.
For eps_y = 2-20 nm they are byte-identical to the production decks in
inputs/decks/round4_cons_*.txt. At 1 nm the formula gives n_m = 14,096,315;
the published 1 nm point is instead the already-simulated ladder rung
n_m = 14,100,000 (ratio 1.0003, inputs/decks/phase9_batch3.txt).

    python3 inputs/make_conservative_decks.py              # print all rows
    python3 inputs/make_conservative_decks.py --check      # compare with inputs/decks/
"""
import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO))
from constants import N_X_CONS as N_X, N_Z_CONS as N_Z, n_m_cons, n_y_cons  # noqa: E402

EPS_NM = [1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20.0]
BEAM = "0.624 0.9 {emitt_y} 12.0 0.12 100 0"   # particles emitt_x emitt_y beta_x beta_y sigma_z offset_y


def row(eps_nm: float, d_y: float) -> str:
    n_y = n_y_cons(d_y)
    n_m = math.ceil(n_m_cons(d_y, n_y, N_X, N_Z))
    return f"{BEAM.format(emitt_y=eps_nm / 1000)} {N_X} {n_y} {N_Z} {n_m} 2 6"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="compare with inputs/decks/round4_cons_*.txt")
    a = ap.parse_args()
    d_y = {float(k): v for k, v in json.loads((REPO / "data" / "D_y_table.json").read_text()).items()}
    rows = {e: row(e, d_y[e]) for e in EPS_NM}
    if not a.check:
        print("# particles emitt_x emitt_y beta_x beta_y sigma_z offset_y n_x n_y n_z n_m integration_method n_t")
        for e in EPS_NM:
            print(rows[e])
        return 0
    bad = 0
    for e in EPS_NM[1:]:
        deck = HERE / "decks" / f"round4_cons_{e:g}nm.txt"
        body = [l.strip() for l in deck.read_text().splitlines() if l.strip() and not l.startswith("#")]
        ok = body == [rows[e]]
        bad += not ok
        print(f"{e:>5g} nm  {'OK  ' if ok else 'DIFF'}  generated: {rows[e]}" + ("" if ok else f"   deck: {body}"))
    print(f"    1 nm  n/a   generated: {rows[1.0]}   (published point: ladder rung n_m = 14100000)")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
