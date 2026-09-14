#!/usr/bin/env python3
"""GUINEA-PIG++ numbers exported for the WarpX code comparison.

Writes two canonical tables that the WarpX analysis reads:

  data/gp_luminosity_for_wx.csv     luminosity per eps_y, recomputed from the raw
                                    .ref files, in three blocks:
      nominal          512x512x25, n_m = 1e5 (tab:code_comparison nominal rows)
      conservative     runs ON the conservative locus (tab:ny_recommendations;
                       frozen constants in constants.py):
                         n_y = 2^ceil(log2(50 D_y^0.402)),
                         n_m = 2.305e-7 D_y^3.300 * 512 * n_y * 64,
                       qualifying iff deck n_y == n_y_cons and
                       deck n_m / n_m_cons in [0.99, 1.10], photons and pairs off
      frozen_extension 40-100 nm at 512x128x64, n_m = 5e4 (NOT conservative)
  data/gp_requirements_for_wx.csv   the tuning-ladder requirements n_y^req and
                                    n_m^req per eps_y, re-exported from
                                    data/gp_calibration_export.csv (written by
                                    GP_CALIBRATION.ipynb) with its sigma_log
                                    uncertainty column; no refit.

Luminosity is lumi_ee [m^-2 per crossing] * 1e-4 * n_b * f_rep, with n_b and
f_rep read from each .ref (factor 1.596 for C3-250), applied once; mean, sample
std (ddof=1) and standard error over every seed on disk for the deck. D_y for
the nominal and conservative rows is data/D_y_table.json (the values printed in
the paper); for the frozen extension it is computed from each run's echoed beam
parameters, D_y = 2 N r_e sigma_z / (gamma sigma_y (sigma_x + sigma_y)).

deck_path is relative to GP_RAW_ROOT (see paths.py).

    python3 export_for_wx.py
"""
from __future__ import annotations

import glob
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from constants import N_X_CONS as N_X, N_Z_CONS as N_Z, Q_N_PRED, Q_P, n_m_cons, n_y_cons  # noqa: E402  (frozen locus)
from paths import DATA_DIR, RAW_ROOT  # noqa: E402
from refparse import GAMMA, parse_ref  # noqa: E402

R_E = 2.8179403262e-15          # classical electron radius [m]
QUALIFY = (0.99, 1.10)

# Every .ref matching one of these patterns is one seed of that row.
NOMINAL_DECKS = {
    0.5: 'output_nm/C3_250/512_512_25_100000_test_C3_250_particles_0.624_emittx_0.9_emitty_0.0005_betax_12.0_betay_0.12_sigmaz_100_offsety_0_integration_method_2_nt6_seed_*.ref',
    1.0: 'output_nm/C3_250/512_512_25_100000_test_C3_250_particles_0.624_emittx_0.9_emitty_0.001_betax_12.0_betay_0.12_sigmaz_100_offsety_0_integration_method_2_nt6_seed_*.ref',
    2.0: 'output_nm/C3_250/512_512_25_100000_test_C3_250_particles_0.624_emittx_0.9_emitty_0.002_betax_12.0_betay_0.12_sigmaz_100_offsety_0_integration_method_2_nt6_seed_*.ref',
    4.0: 'output_nm/C3_250/512_512_25_100000_test_C3_250_particles_0.624_emittx_0.9_emitty_0.004_betax_12.0_betay_0.12_sigmaz_100_offsety_0_integration_method_2_nt6_seed_*.ref',
    8.0: 'output_nm/C3_250/512_512_25_100000_test_C3_250_particles_0.624_emittx_0.9_emitty_0.008_betax_12.0_betay_0.12_sigmaz_100_offsety_0.000000_integration_method_2_nt6_seed_*.ref',
    12.0: 'output_nm/C3_250/512_512_25_100000_test_C3_250_particles_0.624_emittx_0.9_emitty_0.012_betax_12.0_betay_0.12_sigmaz_100_offsety_0_integration_method_2_nt6_seed_*.ref',
    16.0: 'output_nm/C3_250/512_512_25_100000_test_C3_250_particles_0.624_emittx_0.9_emitty_0.016_betax_12.0_betay_0.12_sigmaz_100_offsety_0_integration_method_2_nt6_seed_*.ref',
    20.0: 'output/C3_250/512_512_25_test_C3_250_particles_0.624_emittx_0.9_emitty_0.020_betax_12.0_betay_0.12_sigmaz_100_offsety_0_seed_*.ref',
}
FROZEN_DECKS = {
    40.0: 'output_nm/C3_250/512_128_64_50000_test_C3_250_particles_0.624_emittx_0.9_emitty_0.04_betax_12.0_betay_0.12_sigmaz_100_offsety_0_integration_method_2_nt6_seed_*.ref',
    60.0: 'output_nm/C3_250/512_128_64_50000_test_C3_250_particles_0.624_emittx_0.9_emitty_0.06_betax_12.0_betay_0.12_sigmaz_100_offsety_0_integration_method_2_nt6_seed_*.ref',
    80.0: 'output_nm/C3_250/512_128_64_50000_test_C3_250_particles_0.624_emittx_0.9_emitty_0.08_betax_12.0_betay_0.12_sigmaz_100_offsety_0_integration_method_2_nt6_seed_*.ref',
    100.0: 'output_nm/C3_250/512_128_64_50000_test_C3_250_particles_0.624_emittx_0.9_emitty_0.1_betax_12.0_betay_0.12_sigmaz_100_offsety_0_integration_method_2_nt6_seed_*.ref',
}
# emitt_y exactly as the conservative-run filenames spell it
EMITTY_STR = {0.5: "0.0005", 1.0: "0.001", 2.0: "0.002", 4.0: "0.004", 8.0: "0.008",
              12.0: "0.012", 16.0: "0.016", 20.0: "0.02"}

LUMI_COLS = ["eps_y_nm", "D_y", "block", "n_x", "n_y", "n_z", "n_t", "n_m", "n_m_cons",
             "ratio_nm", "n_seeds", "L_mean_1e34", "L_std_1e34", "L_sem_1e34",
             "do_photons", "do_pairs", "deck_path"]
REQ_COLS = ["quantity", "eps_y_nm", "D_y", "value", "sigma_log", "n_seeds"]


def fmt(x) -> str:
    if x is None:
        return ""
    if isinstance(x, (float, np.floating)):
        return repr(float(x))
    return str(x)


def rel(path: str) -> str:
    return os.path.relpath(path, RAW_ROOT)


def stats(values: list[float]) -> tuple[float, float | None, float | None]:
    a = np.array(values)
    if len(a) < 2:
        return float(a.mean()), None, None
    sd = float(a.std(ddof=1))
    return float(a.mean()), sd, sd / math.sqrt(len(a))


def deck_block(block: str, eps: float, pattern: str, d_y_table: dict) -> dict:
    files = sorted(glob.glob(str(RAW_ROOT / pattern)))
    if not files:
        raise FileNotFoundError(f"{block} eps_y={eps} nm: no .ref matches {pattern}")
    parsed = [parse_ref(f) for f in files]
    ident = {(d["n_x"], d["n_y"], d["n_z"], d["n_t"], d["n_m_1"], d["integration_method"],
              d["do_photons_1"], d["do_pairs"], d["offset_y_nm"], d["emitt_y"]) for d in parsed}
    if len(ident) != 1:
        raise RuntimeError(f"{block} eps_y={eps} nm: runs matching the deck differ: {ident}")
    d0 = parsed[0]
    if block == "frozen_extension":
        n_part = d0["particles"] * 1e10 if d0["particles"] < 1e3 else d0["particles"]
        sx, sy, sz = d0["sigma_x_nm"] * 1e-9, d0["sigma_y_nm"] * 1e-9, d0["sigma_z_um"] * 1e-6
        d_y = 2 * n_part * R_E * sz / (GAMMA * sy * (sx + sy))
    else:
        d_y = d_y_table[eps]
    mean, sd, se = stats([d["lumi_ee_cm2s"] / 1e34 for d in parsed])
    return dict(eps_y_nm=eps, D_y=d_y, block=block, n_x=d0["n_x"], n_y=d0["n_y"], n_z=d0["n_z"],
                n_t=d0["n_t"], n_m=d0["n_m_1"], n_m_cons=None, ratio_nm=None, n_seeds=len(parsed),
                L_mean_1e34=mean, L_std_1e34=sd, L_sem_1e34=se, do_photons=d0["do_photons_1"],
                do_pairs=d0["do_pairs"], deck_path=pattern)


def conservative_block(d_y_table: dict) -> tuple[list[dict], list[str]]:
    rows, notes = [], []
    for eps in sorted(EMITTY_STR):
        d_y = d_y_table[eps]
        nyc = n_y_cons(d_y)
        nmc = n_m_cons(d_y, nyc)
        cands: dict[tuple, dict] = {}
        pat = str(RAW_ROOT / "output_nm" / "C3_250" /
                  f"512_{nyc}_64_*_test_C3_250_particles_0.624_emittx_0.9_emitty_{EMITTY_STR[eps]}_*seed_*.ref")
        for f in sorted(glob.glob(pat)):
            d = parse_ref(f)
            if d["lumi_ee_m2"] is None or (d["n_x"], d["n_y"], d["n_z"]) != (512, nyc, 64):
                continue
            if abs(d["offset_y_nm"] or 0) > 1e-12 or d["n_t"] != 6 or d["integration_method"] != 2:
                continue
            if not (abs(d["cut_x_mult"] - 20) < 1e-3 and abs(d["cut_y_mult"] - 20) < 1e-3
                    and abs(d["cut_z_mult"] - 3.5) < 1e-6):
                continue
            seed = int(re.search(r"_seed_(\d+)\.ref$", f).group(1))
            cands.setdefault((d["n_m_1"], d["do_photons_1"], d["do_pairs"]), {})[seed] = (
                d["lumi_ee_cm2s"] / 1e34, f)
        ok = [(k, v) for k, v in cands.items()
              if QUALIFY[0] <= k[0] / nmc <= QUALIFY[1] and k[1] == 0 and k[2] == 0]
        if not ok:
            on_disk = ", ".join(f"n_m={k[0]} ratio {k[0] / nmc:.4f} ({len(v)} seeds, photons/pairs {k[1]}/{k[2]})"
                                for k, v in sorted(cands.items())) or "none"
            notes.append(f"# EXCLUDED conservative eps_y={eps:g} nm: no run on the locus "
                         f"(n_y_cons={nyc}, n_m_cons={nmc:.1f}; on disk at n_y_cons: {on_disk})")
            continue
        # most seeds first, then the n_m closest to the locus
        (nm, ph, pr), seeds = max(ok, key=lambda kv: (len(kv[1]), -abs(kv[0][0] / nmc - 1)))
        by_seed = [seeds[s] for s in sorted(seeds)]
        mean, sd, se = stats([v[0] for v in by_seed])
        first = sorted(v[1] for v in by_seed)[0]
        rows.append(dict(eps_y_nm=eps, D_y=d_y, block="conservative", n_x=N_X, n_y=nyc, n_z=N_Z, n_t=6,
                         n_m=nm, n_m_cons=nmc, ratio_nm=nm / nmc, n_seeds=len(by_seed),
                         L_mean_1e34=mean, L_std_1e34=sd, L_sem_1e34=se, do_photons=ph, do_pairs=pr,
                         deck_path=re.sub(r"_seed_\d+\.ref$", "_seed_*.ref", rel(first))))
    return rows, notes


def write_luminosity(path) -> None:
    d_y_table = {float(k): v for k, v in json.loads((DATA_DIR / "D_y_table.json").read_text()).items()}
    nominal = [deck_block("nominal", e, p, d_y_table) for e, p in sorted(NOMINAL_DECKS.items())]
    conservative, notes = conservative_block(d_y_table)
    frozen = [deck_block("frozen_extension", e, p, d_y_table) for e, p in sorted(FROZEN_DECKS.items())]
    hdr = [
        "# gp_luminosity_for_wx.csv -- GUINEA-PIG++ C3-250 luminosity for the WarpX code comparison",
        "# written by export_for_wx.py from the raw .ref files; deck_path is relative to GP_RAW_ROOT",
        "# locus: n_y_cons = 2^ceil(log2(50 * D_y^0.402))",
        "# locus: n_m_cons = 2.305e-7 * D_y^3.300 * n_x * n_y_cons * n_z, n_x = 512, n_z = 64",
        "# locus: qualifies for block 'conservative' iff deck n_y == n_y_cons and deck n_m / n_m_cons in [0.99, 1.10]",
        "# locus: fixed n_t = 6, cut multipliers 20/20/3.5, integration_method = 2, offset_y = 0, do_photons = 0, do_pairs = 0",
        "# D_y: data/D_y_table.json for nominal and conservative rows; frozen_extension from each run's echoed beam parameters",
        "# units: L_* in 1e34 cm^-2 s^-1; L = lumi_ee[m^-2 per crossing] * 1e-4 * n_b * f_rep with n_b=133, f_rep=120 read from each .ref (factor 1.596), applied once",
        "# L_std_1e34 = sample standard deviation (ddof=1); L_sem_1e34 = L_std / sqrt(n_seeds)",
        "# blocks: nominal = tab:code_comparison configuration (512x512x25, n_m=1e5); conservative = on-locus runs above; frozen_extension = 40-100 nm, NOT conservative, never merge",
        "# note: the nominal 20 nm row has do_photons=do_pairs=1 (all other nominal rows 0); runs differing only in those flags agree in lumi_ee to < 0.05%",
        "# note: frozen_extension runs are at n_y=128, n_m=50000 -- the pre-campaign 20 nm grid, not the conservative 20 nm grid",
    ] + notes
    with open(path, "w") as fh:
        fh.write("\n".join(hdr) + "\n" + ",".join(LUMI_COLS) + "\n")
        for r in nominal + conservative + frozen:
            fh.write(",".join(fmt(r[c]) for c in LUMI_COLS) + "\n")
    print(f"wrote {path} ({len(nominal)} nominal, {len(conservative)} conservative, {len(frozen)} frozen_extension rows)")


def write_requirements(path) -> None:
    cal = pd.read_csv(DATA_DIR / "gp_calibration_export.csv")
    con = pd.read_csv(DATA_DIR / "gp_constants_export.csv").set_index("name")
    rows = []
    for fam, qty in (("n_y", "n_y_req"), ("n_m", "n_m_req")):
        for r in cal[cal.family == fam].sort_values("eps_y_nm").itertuples():
            rows.append(dict(quantity=qty, eps_y_nm=float(r.eps_y_nm), D_y=float(r.D_y), value=float(r.value),
                             sigma_log=float(r.sigma_log), n_seeds=None))
    published = [("C_y_fit", 26.2, None, "C_y_fit"), ("q_n_fit", 0.402, 0.152, "q_fit"),
                 ("C_m_fit", 1.4e-7, None, "C_m_fit"), ("s_fit", 3.330, 0.250, "s_fit"),
                 ("C_y_cons", 50.0, None, None), ("C_m_cons", 2.305e-7, None, "C_m_cons"),
                 ("s_cons", 3.515, None, "s_cons")]
    hdr = [
        "# gp_requirements_for_wx.csv -- GUINEA-PIG++ tuning-ladder requirements, re-exported unchanged (no refit)",
        "# written by export_for_wx.py",
        "# source: data/gp_calibration_export.csv (written by GP_CALIBRATION.ipynb) for value, D_y and sigma_log",
        "# CHANGED 2026-09-14: the uncertainty column is now sigma_log, the budget behind the published fits. Earlier versions of this file carried sigma_tot from data/error_budget.csv, an older budget that no longer matches the notebook; do not mix the two versions",
        "# n_y_req rows: value = n_y^req (vertical cells)",
        "# n_m_req rows: value = n_m^req / (n_x n_y n_z), i.e. required macroparticles PER CELL, as exported -- multiply by n_x n_y n_z for an absolute n_m",
        "# sigma_log = standard deviation of ln(value): Monte Carlo parameter noise (seeded), window selection and jackknife terms in quadrature (sigma_log_total in GP_CALIBRATION.ipynb), dimensionless",
        f"# q_p = {Q_P} (constants.py): the envelope-integration exponent shared with the WarpX repository; q_n^pred = q_p + 1/4 = {Q_N_PRED!r}",
        "# n_seeds: not recorded in the ladder exports; left empty rather than inferred",
        "# D_y for n_m_req rows is as exported (rounded to 2 decimals in the source); for n_y_req rows at full precision",
    ]
    for name, pub, pub_err, key in published:
        if key is None or key not in con.index:
            hdr.append(f"# constant {name:<9s} published {pub} -- not a ladder fit output (a chosen normalisation); not in data/gp_constants_export.csv")
            continue
        val = float(con.loc[key, "value"])
        err = None if pd.isna(con.loc[key, "error"]) else float(con.loc[key, "error"])
        dec = len(repr(pub).split(".")[1]) if "." in repr(pub) and "e" not in repr(pub) else None
        ok = (round(val, dec) == pub) if dec is not None else abs(val / pub - 1) < 0.005
        if pub_err is not None and err is not None:
            ok = ok and round(err, 3) == pub_err
        note = ""
        hdr.append(f"# constant {name:<9s} published {pub}{' +- ' + str(pub_err) if pub_err is not None else ''}"
                   f" | stored in export {val!r}{' +- ' + repr(err) if err is not None else ''}"
                   f" | reproduces: {'yes' if ok else 'NO'}{note}")
    with open(path, "w") as fh:
        fh.write("\n".join(hdr) + "\n" + ",".join(REQ_COLS) + "\n")
        for r in rows:
            fh.write(",".join(fmt(r[c]) for c in REQ_COLS) + "\n")
    print(f"wrote {path} ({len(rows)} rows)")


def main() -> None:
    write_luminosity(DATA_DIR / "gp_luminosity_for_wx.csv")
    write_requirements(DATA_DIR / "gp_requirements_for_wx.csv")


if __name__ == "__main__":
    main()
