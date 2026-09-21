#!/usr/bin/env python3
"""Regenerate every GUINEA-PIG++ number the manuscript reports, from committed data only.

    python3 reproduce_paper_numbers.py

Writes paper_numbers.json (full precision, every uncertainty with its definition,
every quantity with its provenance) and paper_numbers.md (tables rounded as the
manuscript rounds). Reads nothing outside this repository, uses no network and no
random numbers, and stops with an error rather than emitting a partial or
substituted number if an input is missing or disagrees with the analysis exports.

Contents: luminosity dataset; derived luminosity quantities (conservative/nominal,
H_D, 20 -> 2 nm gain, seed scatter); n_y^req and n_m^req; fitted constants and the
frozen conservative loci; the n_y^req inputs exported for kappa (kappa itself is
computed in the WarpX repository); the recommendation table; D_y; the IPC background
table; the derived background quantities; the n_x convergence study.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import chi2 as chi2_dist

REPO = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO))

INPUTS = {
    "luminosity_export": "data/gp_luminosity_for_wx.csv",
    "requirements_export": "data/gp_calibration_export.csv",
    "requirements_for_warpx": "data/gp_requirements_for_wx.csv",
    "constants_export": "data/gp_constants_export.csv",
    "bib_cache": "data/bib_stats_cache.csv",
    "bib_run_luminosity": "data/bib_run_luminosity.csv",
    "luminosity_snapshot": "data/lumi_extracted.csv",
    "nx_convergence_table": "data/L_vs_nm_by_nx_20nm.csv",
    "D_y_table": "data/D_y_table.json",
    "input_deck": "inputs/acc_lumi_opt_C3_250_no_pairs.dat",
    "beam_deck": "inputs/acc_C3_250_nominal.dat",
    "production_decks": "inputs/decks/round4_cons_20nm.txt",
    "notebook": "GP_CALIBRATION.ipynb",
    "frozen_locus": "constants.py",
}
PARAMETER_BLOCK = "Jim_pars_Aug2023"      # the block inputs/submit_template_S3DF_SLURM.sh runs
RECOMMENDATION_EPS_NM = [0.5, 1.0, 1.5, 2.0, 2.5, 4.0, 8.0, 12.0, 16.0, 20.0]

SD = "sample standard deviation over seeds, ddof = 1"
SE = "standard error of the mean, sample standard deviation (ddof = 1) / sqrt(n_seeds)"
SE_RATIO = ("standard error of the ratio of two means, first-order propagation of the two "
            "standard errors, treated as independent")
SIG_LN_REQ = ("standard deviation of ln(requirement): parameter noise (Monte Carlo over the "
              "Richardson fit, seeded), window selection and jackknife terms in quadrature, "
              "as computed by sigma_log_total in GP_CALIBRATION.ipynb")
SIG_FIT = "one-standard-deviation parameter uncertainty from the weighted least-squares covariance"


class ReproductionError(RuntimeError):
    pass


def fail(msg: str):
    raise ReproductionError(msg)


def unc(value, definition):
    return {"value": value, "definition": definition}


def fnum(x):
    if x is None:
        return None
    x = float(x)
    if not math.isfinite(x):
        fail(f"non-finite value {x} in output")
    return x


def load_inputs() -> dict[str, Path]:
    paths = {k: REPO / v for k, v in INPUTS.items()}
    missing = [INPUTS[k] for k, p in paths.items() if not p.is_file()]
    if missing:
        fail("missing committed input(s): " + ", ".join(missing))
    return paths


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def close(a, b, rel=1e-12, what=""):
    if not (abs(a - b) <= rel * max(abs(a), abs(b), 1e-300)):
        fail(f"{what}: recomputed {a!r} disagrees with committed export {b!r}")


# ── inputs that are not tables ────────────────────────────────────────────────

def notebook_source(nb_path: Path) -> str:
    nb = json.loads(nb_path.read_text())
    return "\n".join("".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code")


def notebook_literal(src: str, name: str) -> float:
    m = re.findall(rf"^\s*{name}\s*=\s*([0-9.eE+-]+)", src, re.M)
    if len(m) != 1:
        fail(f"GP_CALIBRATION.ipynb: expected exactly one literal assignment of {name}, found {len(m)}")
    return float(m[0])


def _kv(body: str) -> dict:
    body = re.sub(r"#[^\n]*", "", body)
    return {k.strip(): v.strip() for k, v in re.findall(r"([A-Za-z_.0-9]+)\s*=\s*([^;]+);", body)}


def deck_block(deck: Path, kind: str, name: str) -> dict:
    m = re.search(rf"\${kind}::\s*{re.escape(name)}\s*\{{(.*?)\}}", deck.read_text(), re.S)
    if not m:
        fail(f"{deck.name}: ${kind}:: {name} block not found")
    return _kv(m.group(1))


def production_beam_columns(decks_dir: Path) -> set:
    """(particles, emitt_x, beta_x, beta_y, sigma_z, offset_y) of every row of every committed deck."""
    cols = set()
    for f in sorted(decks_dir.glob("*.txt")):
        for line in f.read_text().splitlines():
            t = line.split()
            if not t or line.lstrip().startswith("#"):
                continue
            v = [float(x) for x in t]
            cols.add((v[0], v[1], v[3], v[4], v[5], v[6]))
    return cols


def cut_multiplier(expr: str) -> float:
    m = re.fullmatch(r"([0-9.]+)\*sigma_[xyz]\.1", expr)
    if not m:
        fail(f"cut expression {expr!r} is not of the form k*sigma")
    return float(m.group(1))


# ── statistics ────────────────────────────────────────────────────────────────

def ratio_se(a, sa, b, sb):
    r = a / b
    return r, r * math.sqrt((sa / a) ** 2 + (sb / b) ** 2)


def wls(x, y, s):
    """Weighted least squares y = a + b x with weights 1/s^2 (the notebook's closed form)."""
    w = 1.0 / np.asarray(s) ** 2
    sw, sx, sy = w.sum(), (w * x).sum(), (w * y).sum()
    sxx, sxy = (w * x * x).sum(), (w * x * y).sum()
    det = sw * sxx - sx ** 2
    b = (sw * sxy - sx * sy) / det
    a = (sxx * sy - sx * sxy) / det
    chi2 = float((w * (y - a - b * x) ** 2).sum())
    return dict(intercept=float(a), slope=float(b), var_intercept=float(sxx / det),
                var_slope=float(sw / det), chi2=chi2, ndf=len(x) - 2)


def richardson3(n, L):
    (n1, n2, n3), (L1, L2, L3) = n, L
    r = (L1 - L2) / (L2 - L3)
    f = lambda p: (n1 ** -p - n2 ** -p) / (n2 ** -p - n3 ** -p) - r
    p = brentq(f, 1e-3, 5.0)
    C = (L2 - L3) / (n2 ** -p - n3 ** -p)
    return L3 - C * n3 ** -p, p


# ── the report ────────────────────────────────────────────────────────────────

def main() -> None:
    P = load_inputs()
    out: dict = {"about": {
        "script": "reproduce_paper_numbers.py",
        "simulator": "GUINEA-PIG++, C3-250",
        "luminosity_units": "1e34 cm^-2 s^-1",
        "determinism": "no random numbers are drawn; per-emittance requirement uncertainties are read from the notebook export, whose Monte Carlo term used rng_seed = 23",
        "inputs_sha256": {INPUTS[k]: sha256(p) for k, p in P.items()},
    }}

    nb_src = notebook_source(P["notebook"])
    # The C3_250 beam block of the campaign input file holds placeholders that every scan row
    # overrides (particles, emittances, beta*, sigma_z), so the beam actually run is taken from
    # inputs/acc_C3_250_nominal.dat and checked against every committed production deck row.
    acc = deck_block(P["beam_deck"], "ACCELERATOR", "C3_250")
    par = deck_block(P["input_deck"], "PARAMETERS", PARAMETER_BLOCK)
    deck_cols = production_beam_columns(P["production_decks"].parent)
    import constants as frozen  # noqa: E402  (repository file, listed in INPUTS)

    # physical constants exactly as the notebook defines them (cell defining D_y_full)
    R_E = notebook_literal(nb_src, "R_E")
    gamma_m = re.search(r"^\s*GAMMA\s*=\s*125e9\s*/\s*([0-9.eE+-]+)", nb_src, re.M)
    if not gamma_m:
        fail("GP_CALIBRATION.ipynb: GAMMA = 125e9 / m_e c^2 not found")
    GAMMA = 125e9 / float(gamma_m.group(1))
    # q_p (envelope-integration exponent, shared with the WarpX repository) and
    # q_n^pred = q_p + 1/4 are defined once, in constants.py; the notebook imports them.
    Q_P, Q_N = float(frozen.Q_P), float(frozen.Q_N_PRED)
    if not re.search(r"from constants import Q_P as Q_PINCH, Q_N_PRED as Q_NY", nb_src):
        fail("GP_CALIBRATION.ipynb does not take q_p and q_n^pred from constants.py")
    beam = dict(energy_GeV=float(acc["energy"]), particles_per_bunch=float(acc["particles"]) * 1e10,
                beta_x_mm=float(acc["beta_x"]), beta_y_mm=float(acc["beta_y"]),
                emitt_x_normalised_mm_mrad=float(acc["emitt_x"]), sigma_z_um=float(acc["sigma_z"]),
                n_b=int(float(acc["n_b"])), f_rep_Hz=float(acc["f_rep"]))
    if beam["energy_GeV"] != 125.0:
        fail("input deck beam energy is not 125 GeV")
    expected = (float(acc["particles"]), beam["emitt_x_normalised_mm_mrad"], beam["beta_x_mm"], beam["beta_y_mm"], beam["sigma_z_um"], 0.0)
    # Every deck row must run the nominal beam in every column except eps_x: the eps_x sweep of
    # the luminosity-vs-BIB figure (inputs/decks/pair_dump_lumi_bib_grid.txt) varies eps_x only.
    other = {(c[0], c[2], c[3], c[4], c[5]) for c in deck_cols}
    if other != {(expected[0], expected[2], expected[3], expected[4], expected[5])} or expected not in deck_cols:
        fail(f"production decks run beam columns {sorted(deck_cols)}, inconsistent with the beam of {INPUTS['beam_deck']} {expected}")
    grid_defaults = dict(n_t=int(par["n_t"]), integration_method=int(par["integration_method"]),
                         cut_x_sigma=cut_multiplier(par["cut_x"]), cut_y_sigma=cut_multiplier(par["cut_y"]),
                         cut_z_sigma=cut_multiplier(par["cut_z"]))

    def D_y(eps_nm: float) -> float:
        sx = math.sqrt(beam["emitt_x_normalised_mm_mrad"] * 1e-6 * beam["beta_x_mm"] * 1e-3 / GAMMA)
        sy = math.sqrt(eps_nm * 1e-9 * beam["beta_y_mm"] * 1e-3 / GAMMA)
        return 2.0 * R_E * beam["particles_per_bunch"] * beam["sigma_z_um"] * 1e-6 / (GAMMA * sy * (sx + sy))

    def L_geom_1e34(eps_nm: float) -> float:
        sx = math.sqrt(beam["emitt_x_normalised_mm_mrad"] * 1e-6 * beam["beta_x_mm"] * 1e-3 / GAMMA)
        sy = math.sqrt(eps_nm * 1e-9 * beam["beta_y_mm"] * 1e-3 / GAMMA)
        N = beam["particles_per_bunch"]
        return N * N * beam["n_b"] * beam["f_rep_Hz"] / (4 * math.pi * sx * sy) * 1e-4 / 1e34

    # ── 7. disruption parameter ──
    d_table = {float(k): float(v) for k, v in json.loads(P["D_y_table"].read_text()).items()}
    eps_all = sorted(set(RECOMMENDATION_EPS_NM) | set(d_table))
    out["disruption_parameter"] = {
        "expression": "D_y = 2 r_e N sigma_z / (gamma sigma_y (sigma_x + sigma_y)), sigma_{x,y} = sqrt(eps*_{x,y} beta*_{x,y} / gamma)",
        "parameters": {"r_e_m": R_E, "gamma": GAMMA, "N": beam["particles_per_bunch"],
                       "sigma_z_m": beam["sigma_z_um"] * 1e-6, "eps_x_normalised_m_rad": beam["emitt_x_normalised_mm_mrad"] * 1e-6,
                       "beta_x_m": beam["beta_x_mm"] * 1e-3, "beta_y_m": beam["beta_y_mm"] * 1e-3},
        "provenance": f"expression and r_e, gamma from GP_CALIBRATION.ipynb (D_y_full); beam from {INPUTS['beam_deck']}, identical to the beam columns of every deck in inputs/decks/",
        "rows": [{"eps_y_nm": e, "D_y": D_y(e),
                  "D_y_deck_value": d_table.get(e),
                  "D_y_deck_value_definition": "three-decimal D_y used to write the conservative decks (data/D_y_table.json)" if e in d_table else None}
                 for e in eps_all],
    }
    for e, v in d_table.items():
        if abs(round(D_y(e), 3) - v) > 1.5e-3:
            fail(f"D_y table value at {e} nm ({v}) is inconsistent with the expression ({D_y(e)})")

    # ── 1. luminosity dataset ──
    lum = pd.read_csv(P["luminosity_export"], comment="#")
    excluded_header = [l[1:].strip() for l in P["luminosity_export"].read_text().splitlines() if l.startswith("# EXCLUDED")]
    need = {"nominal", "conservative"}
    if not need <= set(lum.block):
        fail(f"luminosity export lacks block(s) {need - set(lum.block)}")
    rows = []
    for r in lum[lum.block.isin(need)].sort_values(["block", "eps_y_nm"]).itertuples():
        if not (r.L_std_1e34 > 0 and r.n_seeds >= 2):
            fail(f"{r.block} {r.eps_y_nm} nm: fewer than 2 seeds or zero scatter")
        close(r.L_std_1e34 / math.sqrt(r.n_seeds), r.L_sem_1e34, 1e-9, f"{r.block} {r.eps_y_nm} nm standard error")
        rows.append({
            "configuration": r.block, "eps_y_nm": float(r.eps_y_nm),
            "L_mean_1e34": float(r.L_mean_1e34),
            "L_std_1e34": unc(float(r.L_std_1e34), SD),
            "L_sem_1e34": unc(float(r.L_sem_1e34), SE),
            "n_seeds": int(r.n_seeds),
            "grid": {"n_x": int(r.n_x), "n_y": int(r.n_y), "n_z": int(r.n_z), "n_t": int(r.n_t), "n_m": int(r.n_m),
                     "integration_method": grid_defaults["integration_method"],
                     "cut_multipliers_sigma": [grid_defaults["cut_x_sigma"], grid_defaults["cut_y_sigma"], grid_defaults["cut_z_sigma"]]},
            "do_photons": int(r.do_photons), "do_pairs": int(r.do_pairs),
            "runs": r.deck_path,
        })
    out["luminosity_dataset"] = {
        "provenance": f"{INPUTS['luminosity_export']} (written by export_for_wx.py from the raw .ref files); cut multipliers, n_t and solver from {INPUTS['input_deck']} $PARAMETERS:: {PARAMETER_BLOCK}, which the scan rows do not override",
        "selection": "nominal: every seed of the 512x512x25, n_m = 1e5 deck at that eps_y. conservative: runs with deck n_y == n_y_cons and deck n_m / n_m_cons in [0.99, 1.10], photons and pairs off, all seeds of the qualifying deck",
        "excluded": ["frozen_extension block of the export: not one of the GUINEA-PIG++ configurations the manuscript reports"] + excluded_header,
        "rows": rows,
    }
    L = {(r["configuration"], r["eps_y_nm"]): r for r in rows}

    # ── 2. derived luminosity quantities ──
    ratios, hd = [], []
    for (cfg, e), r in sorted(L.items()):
        lg = L_geom_1e34(e)
        hd.append({"configuration": cfg, "eps_y_nm": e, "H_D": r["L_mean_1e34"] / lg, "L_geom_1e34": lg,
                   "relative_seed_scatter": unc(r["L_std_1e34"]["value"] / r["L_mean_1e34"], "sample standard deviation (ddof = 1) divided by the mean luminosity")})
        if cfg == "conservative" and ("nominal", e) in L:
            n = L[("nominal", e)]
            v, s = ratio_se(r["L_mean_1e34"], r["L_sem_1e34"]["value"], n["L_mean_1e34"], n["L_sem_1e34"]["value"])
            ratios.append({"eps_y_nm": e, "conservative_over_nominal": v, "standard_error": unc(s, SE_RATIO)})
    c2, c20 = L[("conservative", 2.0)], L[("conservative", 20.0)]
    g, gs = ratio_se(c2["L_mean_1e34"], c2["L_sem_1e34"]["value"], c20["L_mean_1e34"], c20["L_sem_1e34"]["value"])
    out["derived_luminosity"] = {
        "provenance": "computed from the luminosity_dataset rows of this file",
        "conservative_over_nominal": ratios,
        "enhancement_factor": {
            "definition": "H_D = L / L_geom, L_geom = N^2 n_b f_rep / (4 pi sigma_x sigma_y), no hourglass factor, sigma from normalised emittance and beta*",
            "parameters": beam | {"gamma": GAMMA},
            "provenance": f"beam parameters from {INPUTS['beam_deck']} (the beam of every production deck); gamma as in GP_CALIBRATION.ipynb",
            "rows": hd},
        "conservative_luminosity_gain_2nm_over_20nm": {
            "value": g, "standard_error": unc(gs, SE_RATIO),
            "definition": "L(conservative, 2 nm) / L(conservative, 20 nm), from the conservative runs with photons and pairs off",
            "n_seeds_2nm": c2["n_seeds"], "n_seeds_20nm": c20["n_seeds"]},
    }

    # ── 3. requirement extractions ──
    cal = pd.read_csv(P["requirements_export"])
    ny = cal[cal.family == "n_y"].sort_values("eps_y_nm")
    nm = cal[cal.family == "n_m"].sort_values("eps_y_nm")
    if len(ny) < 3 or len(nm) < 3:
        fail("requirements export has fewer than three n_y or n_m rows")
    for r in ny.itertuples():
        close(D_y(r.eps_y_nm), r.D_y, 1e-9, f"D_y at {r.eps_y_nm} nm")
    out["requirements"] = {
        "provenance": f"{INPUTS['requirements_export']} (written by GP_CALIBRATION.ipynb): the +-5 % crossing of each Richardson-fitted tuning ladder built from {INPUTS['luminosity_snapshot']}",
        "n_y_req": [{"eps_y_nm": float(r.eps_y_nm), "D_y": float(r.D_y), "value": float(r.value),
                     "sigma_ln": unc(float(r.sigma_log), SIG_LN_REQ)} for r in ny.itertuples()],
        "n_m_req_per_cell": [{"eps_y_nm": float(r.eps_y_nm), "D_y": float(r.D_y), "value": float(r.value),
                              "definition": "n_m^req / (n_x n_y n_z), required macroparticles per grid cell",
                              "sigma_ln": unc(float(r.sigma_log), SIG_LN_REQ)} for r in nm.itertuples()],
        "coverage": {"n_y_req_eps_nm": [float(e) for e in ny.eps_y_nm], "n_m_req_eps_nm": [float(e) for e in nm.eps_y_nm]},
        "note_D_y": "n_m rows carry D_y rounded to 2 decimals, as in the notebook table the n_m fit reads",
    }

    # ── 4. fitted constants and frozen loci ──
    con = pd.read_csv(P["constants_export"]).set_index("name")
    cv = lambda n: float(con.loc[n, "value"])
    ce = lambda n: float(con.loc[n, "error"])
    x, y, s = np.log(ny.D_y.values), np.log(ny.value.values), ny.sigma_log.values
    fy = wls(x, y, s)
    close(fy["slope"], cv("q_fit"), what="q_fit"); close(math.sqrt(fy["var_slope"]), ce("q_fit"), what="sigma q_fit")
    close(math.exp(fy["intercept"]), cv("C_y_fit"), what="C_y_fit")
    w = 1 / s ** 2
    lnc = float((w * (y - Q_N * x)).sum() / w.sum())
    chi2_c = float((w * (y - lnc - Q_N * x) ** 2).sum())
    xm, ym, sm = np.log(nm.D_y.values), np.log(nm.value.values), nm.sigma_log.values
    fm = wls(xm, ym, sm)
    close(fm["slope"], cv("b_fit"), what="b_fit"); close(fm["slope"] + Q_P, cv("s_fit"), what="s_fit")
    close(math.sqrt(fm["var_slope"]), ce("s_fit"), what="sigma s_fit")
    close(math.exp(fm["intercept"]), cv("C_m_fit"), what="C_m_fit")
    close(Q_P, cv("q_pred"), what="q_p")
    out["fitted_constants"] = {
        "method": "weighted least squares in ln-ln space, weights 1/sigma_ln^2 from the requirements section; chi2 = sum of weighted squared residuals",
        "provenance": f"recomputed from {INPUTS['requirements_export']}; checked against {INPUTS['constants_export']} to 1e-12",
        "n_y_free_fit": {"law": "n_y^req = C_y D_y^q_n",
                         "q_n": fy["slope"], "sigma_q_n": unc(math.sqrt(fy["var_slope"]), SIG_FIT),
                         "C_y": math.exp(fy["intercept"]), "sigma_ln_C_y": unc(math.sqrt(fy["var_intercept"]), SIG_FIT),
                         "chi2": fy["chi2"], "ndf": fy["ndf"], "p_value": float(chi2_dist.sf(fy["chi2"], fy["ndf"])),
                         "q_n_minus_prediction_in_sigma": (fy["slope"] - Q_N) / math.sqrt(fy["var_slope"]),
                         "q_n_over_sigma": fy["slope"] / math.sqrt(fy["var_slope"])},
        "n_y_constrained_fit": {"law": "n_y^req = C_y D_y^q_n with q_n fixed", "q_n_fixed": Q_N,
                                "q_n_definition": f"q_p + 1/4 with q_p = {Q_P} from constants.py",
                                "C_y": math.exp(lnc), "sigma_ln_C_y": unc(1 / math.sqrt(w.sum()), SIG_FIT),
                                "chi2": chi2_c, "ndf": len(x) - 1, "p_value": float(chi2_dist.sf(chi2_c, len(x) - 1))},
        "n_m_fit": {"law": "n_m^req / (n_x n_y n_z) = C_m D_y^(s - q_p)",
                    "slope_s_minus_q_p": fm["slope"], "q_p": Q_P, "s": fm["slope"] + Q_P,
                    "sigma_s": unc(math.sqrt(fm["var_slope"]), SIG_FIT),
                    "C_m": math.exp(fm["intercept"]), "sigma_C_m": unc(math.exp(fm["intercept"]) * math.sqrt(fm["var_intercept"]), SIG_FIT + ", propagated to C_m = exp(intercept)"),
                    "chi2": fm["chi2"], "ndf": fm["ndf"], "p_value": float(chi2_dist.sf(fm["chi2"], fm["ndf"]))},
        "conservative_loci": {
            "status": "FROZEN: the values used to define the production runs, fixed when those runs were submitted; not recomputed from the fits above",
            "provenance": INPUTS["frozen_locus"],
            "n_y_cons": {"law": "n_y = 2^ceil(log2(C_y_cons D_y^q_n_cons))", "C_y_cons": frozen.C_Y_CONS, "q_n_cons": frozen.Q_N_CONS},
            "n_m_cons": {"law": "n_m = C_m_cons D_y^(s_cons - q_p) n_x n_y n_z", "C_m_cons": frozen.C_M_CONS,
                         "s_cons_minus_q_p": frozen.NM_EXPONENT_CONS,
                         "s_cons": frozen.NM_EXPONENT_CONS + Q_P, "s_cons_definition": f"s_cons_minus_q_p + q_p with q_p = {Q_P} from constants.py (envelope-integration exponent)",
                         "n_x": frozen.N_X_CONS, "n_z": frozen.N_Z_CONS}},
    }

    # ── 5. inputs exported for kappa (kappa is computed in the WarpX repository) ──
    wxr = pd.read_csv(P["requirements_for_warpx"], comment="#")
    ki = wxr[wxr.quantity == "n_y_req"].sort_values("eps_y_nm")
    if not (list(ki.eps_y_nm) == list(ny.eps_y_nm) and np.array_equal(ki.value.values, ny.value.values)
            and np.array_equal(ki.D_y.values, ny.D_y.values) and np.array_equal(ki.sigma_log.values, ny.sigma_log.values)):
        fail(f"{INPUTS['requirements_for_warpx']} n_y_req rows differ from {INPUTS['requirements_export']}")
    out["kappa_inputs_exported_to_warpx"] = {
        "note": "kappa for both simulators is computed in the WarpX repository from this export; this repository does not compute kappa",
        "provenance": f"{INPUTS['requirements_for_warpx']} (written by export_for_wx.py), identical to the n_y rows of {INPUTS['requirements_export']}",
        "rows": [{"eps_y_nm": float(r.eps_y_nm), "D_y": float(r.D_y), "n_y_req": float(r.value),
                  "sigma_ln": unc(float(r.sigma_log), SIG_LN_REQ)} for r in ki.itertuples()],
    }

    # ── 6. recommendation table ──
    band = (float(ny.D_y.min()), float(ny.D_y.max()))
    rec = []
    for e in RECOMMENDATION_EPS_NM:
        d = D_y(e)
        n_y = frozen.n_y_cons(d)
        n_m = frozen.n_m_cons(d, n_y)
        rec.append({"eps_y_nm": e, "D_y": d, "n_y_rec": int(n_y), "n_m_rec": n_m,
                    "outside_fitted_D_y_band": not (band[0] * (1 - 1e-9) <= d <= band[1] * (1 + 1e-9))})
    out["recommendation_table"] = {
        "definition": "n_y_rec from the frozen n_y locus (rounded up to a power of two), n_m_rec from the frozen n_m locus at that n_y, n_x = 512, n_z = 64",
        "provenance": f"{INPUTS['frozen_locus']} with D_y from the disruption_parameter expression",
        "fitted_D_y_band": list(band), "rows": rec,
    }

    # ── 8. IPC background table ──
    bib = pd.read_csv(P["bib_cache"])
    allc = bib[bib.channel == "ALL"]
    ipc_rows = []
    for (ds, e), grp in sorted(allc.groupby(["dataset", "eps_y_nm"])):
        if grp.seed.duplicated().any():
            fail(f"BIB cache: duplicate seed rows for {ds} {e} nm")
        n = len(grp)
        ipc_rows.append({"dataset": ds, "eps_y_nm": float(e),
                         "grid": {"n_x": int(grp.n_x.iloc[0]), "n_y": int(grp.n_y.iloc[0]), "n_z": int(grp.n_z.iloc[0]), "n_m": int(grp.n_m.iloc[0])},
                         "produced_mean": float(grp.n_pairs.mean()), "produced_std": unc(float(grp.n_pairs.std(ddof=1)), SD),
                         "reaching_mean": float(grp.n_reach.mean()), "reaching_std": unc(float(grp.n_reach.std(ddof=1)), SD),
                         "n_seeds": n, "seeds": sorted(int(v) for v in grp.seed)})
    if len(ipc_rows) != 8:
        fail(f"BIB cache yields {len(ipc_rows)} rows, expected 8")
    out["ipc_background_table"] = {
        "provenance": f"{INPUTS['bib_cache']} (bib_stats_table.py, pair dumps resolved by find_dumps), channel ALL, per bunch crossing",
        "reach_definition": "pairs whose helix from the IP reaches the SiD-o2-v04 barrel (B = 5 T, r_det = 14 mm, z_max = 76 mm)",
        "datasets": "nominal: 512x512x25, n_m = 1e5. tuned: the conservative locus, n_x = 512, n_z = 64, n_y and n_m per row",
        "rows": ipc_rows,
    }

    # ── 9. derived background quantities ──
    T = {(r["dataset"], r["eps_y_nm"]): r for r in ipc_rows}
    t2, t20 = T[("tuned", 2.0)], T[("tuned", 20.0)]
    se = lambda r, k: r[k]["value"] / math.sqrt(r["n_seeds"])
    ipc_r, ipc_s = ratio_se(t2["produced_mean"], se(t2, "produced_std"), t20["produced_mean"], se(t20, "produced_std"))
    rch_r, rch_s = ratio_se(t2["reaching_mean"], se(t2, "reaching_std"), t20["reaching_mean"], se(t20, "reaching_std"))
    runl = pd.read_csv(P["bib_run_luminosity"], comment="#")
    bib_lumi = {}
    for key, r in (("2nm", t2), ("20nm", t20)):
        g_ = r["grid"]
        sel = runl[(runl.dataset == r["dataset"]) & np.isclose(runl.eps_y_nm, r["eps_y_nm"])
                   & (runl.n_x == g_["n_x"]) & (runl.n_y == g_["n_y"]) & (runl.n_z == g_["n_z"])
                   & (runl.n_m == g_["n_m"]) & runl.seed.isin(r["seeds"])
                   & runl.n_pairs.notna() & (runl.n_pairs > 0)]
        per_seed = sel.groupby("seed").lumi_ee.agg(["min", "max"])
        if len(per_seed) != r["n_seeds"]:
            fail(f"BIB-run luminosity at {key}: {len(per_seed)} of {r['n_seeds']} seeds in {INPUTS['bib_run_luminosity']}")
        if not np.allclose(per_seed["min"], per_seed["max"], rtol=0, atol=0):
            fail(f"BIB-run luminosity at {key}: two pair-production runs claim one seed and disagree")
        v = per_seed["min"].values / 1e34
        bib_lumi[key] = (float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v))), len(v))
    lr, ls = ratio_se(bib_lumi["2nm"][0], bib_lumi["2nm"][1], bib_lumi["20nm"][0], bib_lumi["20nm"][1])
    def reduction(nr, ns):
        val = 1 - nr / lr
        err = math.sqrt((ns / lr) ** 2 + (nr * ls / lr ** 2) ** 2)
        return 100 * val, 100 * err
    pr, pe = reduction(ipc_r, ipc_s)
    rr, re_ = reduction(rch_r, rch_s)
    out["derived_background"] = {
        "provenance": f"tuned rows of ipc_background_table (this file); BIB-run luminosity from {INPUTS['bib_run_luminosity']}: read from the .ref of each BIB run itself, paired to that run's pair dump, on the same grid and ten seeds as the BIB rows, offset_y = 0, with pair production on (n_pairs recorded in the .ref)",
        "excluded": ["pairs-off runs on the same grid and seeds: the BIB runs have pairs and photons on, so their own luminosity is the one tab:bib_yields normalises by"],
        "produced_ratio_2nm_over_20nm": {"value": ipc_r, "standard_error": unc(ipc_s, SE_RATIO)},
        "reaching_ratio_2nm_over_20nm": {"value": rch_r, "standard_error": unc(rch_s, SE_RATIO)},
        "BIB_RUN_luminosity_ratio_2nm_over_20nm": {
            "value": lr, "standard_error": unc(ls, SE_RATIO),
            "definition": "luminosity of the conservative-locus BIB runs themselves (photons and pairs on), NOT the conservative-run gain of derived_luminosity",
            "L_2nm_1e34": bib_lumi["2nm"][0], "L_20nm_1e34": bib_lumi["20nm"][0], "n_seeds": bib_lumi["2nm"][2]},
        "produced_per_luminosity_reduction_pct": {"value": pr, "uncertainty": unc(pe, "first-order propagation of the standard errors of the produced ratio and the BIB-run luminosity ratio, treated as independent"),
                                                   "definition": "100 (1 - produced ratio / BIB-run luminosity ratio)"},
        "reaching_per_luminosity_reduction_pct": {"value": rr, "uncertainty": unc(re_, "first-order propagation of the standard errors of the reaching ratio and the BIB-run luminosity ratio, treated as independent"),
                                                  "definition": "100 (1 - reaching ratio / BIB-run luminosity ratio)"},
    }

    # ── 10. n_x convergence at 20 nm ──
    tab = pd.read_csv(P["nx_convergence_table"])
    conv = []
    for nx_, grp in tab.groupby("n_x"):
        s_ = grp.set_index("n_m").L
        if not {1e6, 5e6, 1e7} <= set(s_.index):
            fail(f"n_x convergence table: n_x = {nx_} lacks n_m = 1e6, 5e6 or 1e7")
        top = s_.loc[[1e6, 5e6, 1e7]]
        step = abs(top.iloc[2] - top.iloc[1]) / top.iloc[2]
        try:
            Linf, _ = richardson3(top.index.values, top.values)
        except ValueError:
            Linf = float(top.iloc[2])
        conv.append({"n_x": int(nx_), "L_converged_1e34": float(Linf), "L_at_n_m_1e7_1e34": float(top.iloc[2]),
                     "n_seeds": 1, "converged": bool(step < 1e-3)})
    good = [c for c in conv if c["converged"]]
    mean_good = float(np.mean([c["L_converged_1e34"] for c in good]))
    spread = max(abs(a["L_converged_1e34"] - b["L_converged_1e34"]) for a in good for b in good) / mean_good
    out["nx_convergence_20nm"] = {
        "provenance": f"{INPUTS['nx_convergence_table']}; n_y = 512, n_z = 64, eps_y = 20 nm, seed 100805",
        "definition": "L_converged: exact three-point Richardson extrapolation L = L_inf + C n_m^-p through n_m = 1e6, 5e6, 1e7 (the n_m = 1e7 value where the three points admit no decaying solution); converged: relative change from 5e6 to 1e7 below 0.1 %",
        "rows": conv,
        "agreement_between_converged_n_x": {"max_pairwise_relative_difference": spread,
                                            "definition": "largest |L_converged(i) - L_converged(j)| over converged n_x, divided by their mean"},
        "excluded": [f"n_x = {c['n_x']}: not converged at n_m = 1e7" for c in conv if not c["converged"]],
    }

    # ── write ──
    js = json.dumps(out, indent=1, allow_nan=False, ensure_ascii=False) + "\n"
    (REPO / "paper_numbers.json").write_text(js)
    (REPO / "paper_numbers.md").write_text(render_md(out))
    print("wrote paper_numbers.json and paper_numbers.md")


# ── markdown ─────────────────────────────────────────────────────────────────

def sf(x, n=2):
    if x == 0:
        return "0"
    e = math.floor(math.log10(abs(x)))
    m = round(x / 10 ** e, n - 1)
    if abs(m) >= 10:
        m, e = m / 10, e + 1
    return f"{m:.{n-1}f}×10^{e}"


def render_md(o: dict) -> str:
    M = ["# GUINEA-PIG++ numbers reported in the paper",
         "",
         "Generated by `reproduce_paper_numbers.py` from committed data only. Full precision, uncertainty definitions and provenance are in `paper_numbers.json`. Luminosities in 10³⁴ cm⁻² s⁻¹. Uncertainties are sample standard deviations (ddof = 1) unless marked SE (standard error).",
         ""]
    ld = o["luminosity_dataset"]
    M += ["## 1. Luminosity dataset", "", "| configuration | ε_y [nm] | L | std | SE | seeds | n_x×n_y×n_z | n_t | n_m | cuts [σ] |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in ld["rows"]:
        g = r["grid"]
        M.append(f"| {r['configuration']} | {r['eps_y_nm']:g} | {r['L_mean_1e34']:.3f} | {r['L_std_1e34']['value']:.3f} | {r['L_sem_1e34']['value']:.3f} | {r['n_seeds']} | {g['n_x']}×{g['n_y']}×{g['n_z']} | {g['n_t']} | {g['n_m']:,} | {'/'.join(f'{c:g}' for c in g['cut_multipliers_sigma'])} |")
    M += ["", "Excluded: " + "; ".join(ld["excluded"]), ""]
    dl = o["derived_luminosity"]
    M += ["## 2. Derived luminosity quantities", "", "| ε_y [nm] | conservative / nominal | SE |", "|---|---|---|"]
    M += [f"| {r['eps_y_nm']:g} | {r['conservative_over_nominal']:.3f} | {r['standard_error']['value']:.3f} |" for r in dl["conservative_over_nominal"]]
    M += ["", f"Luminosity gain 20 → 2 nm, conservative runs: **{dl['conservative_luminosity_gain_2nm_over_20nm']['value']:.3f} ± {dl['conservative_luminosity_gain_2nm_over_20nm']['standard_error']['value']:.3f}** (SE)", "",
          f"H_D = L / L_geom, {dl['enhancement_factor']['definition'].split(', ', 1)[1]}.", "", "| configuration | ε_y [nm] | L_geom | H_D | std/mean |", "|---|---|---|---|---|"]
    M += [f"| {r['configuration']} | {r['eps_y_nm']:g} | {r['L_geom_1e34']:.3f} | {r['H_D']:.3f} | {100*r['relative_seed_scatter']['value']:.2f} % |" for r in dl["enhancement_factor"]["rows"]]
    rq = o["requirements"]
    M += ["", "## 3. Requirements", "", "| ε_y [nm] | D_y | n_y^req | σ(ln) |", "|---|---|---|---|"]
    M += [f"| {r['eps_y_nm']:g} | {r['D_y']:.1f} | {r['value']:.1f} | {r['sigma_ln']['value']:.3f} |" for r in rq["n_y_req"]]
    M += ["", "| ε_y [nm] | D_y | n_m^req/(n_x n_y n_z) | σ(ln) |", "|---|---|---|---|"]
    M += [f"| {r['eps_y_nm']:g} | {r['D_y']:.1f} | {sf(r['value'],3)} | {r['sigma_ln']['value']:.3f} |" for r in rq["n_m_req_per_cell"]]
    fc = o["fitted_constants"]; f1, f2, f3 = fc["n_y_free_fit"], fc["n_y_constrained_fit"], fc["n_m_fit"]
    lc = fc["conservative_loci"]
    M += ["", "## 4. Fitted constants", "",
          "| fit | constant | value | χ²/ndf |", "|---|---|---|---|",
          f"| n_y free | q_n | {f1['q_n']:.3f} ± {f1['sigma_q_n']['value']:.3f} | {f1['chi2']:.2f}/{f1['ndf']} |",
          f"| n_y free | C_y | {f1['C_y']:.1f} | |",
          f"| n_y, q_n = {f2['q_n_fixed']} | C_y | {f2['C_y']:.2f} | {f2['chi2']:.2f}/{f2['ndf']} |",
          f"| n_m | s | {f3['s']:.3f} ± {f3['sigma_s']['value']:.3f} | {f3['chi2']:.2f}/{f3['ndf']} |",
          f"| n_m | C_m | {sf(f3['C_m'])} | |",
          "", f"q_n free fit: {f1['q_n_minus_prediction_in_sigma']:+.2f} σ from q_n = {f2['q_n_fixed']}, {f1['q_n_over_sigma']:.2f} σ from 0.", "",
          "Conservative loci, **frozen** values used to define the production runs (not refit):", "",
          f"- n_y^cons = 2^⌈log₂({lc['n_y_cons']['C_y_cons']:g}·D_y^{lc['n_y_cons']['q_n_cons']})⌉",
          f"- n_m^cons = {sf(lc['n_m_cons']['C_m_cons'],4)}·D_y^{lc['n_m_cons']['s_cons_minus_q_p']:.3f}·n_x n_y n_z, i.e. s_cons = {lc['n_m_cons']['s_cons']:.3f} with q_p = {f3['q_p']}", ""]
    k = o["kappa_inputs_exported_to_warpx"]
    M += ["## 5. Inputs exported for κ", "", k["note"].capitalize() + ".", "", "| ε_y [nm] | D_y | n_y^req | σ(ln) |", "|---|---|---|---|"]
    M += [f"| {r['eps_y_nm']:g} | {r['D_y']:.1f} | {r['n_y_req']:.1f} | {r['sigma_ln']['value']:.3f} |" for r in k["rows"]]
    M += [""]
    rt = o["recommendation_table"]
    M += ["## 6. Recommendation table", "", "| ε_y [nm] | D_y | n_y^rec | n_m^rec |", "|---|---|---|---|"]
    M += [f"| {r['eps_y_nm']:g} | {r['D_y']:.1f} | {r['n_y_rec']} | {sf(r['n_m_rec'])}{' (outside fitted D_y band)' if r['outside_fitted_D_y_band'] else ''} |" for r in rt["rows"]]
    dp = o["disruption_parameter"]
    M += ["", "## 7. D_y", "", f"`{dp['expression']}`", "", "| ε_y [nm] | D_y | deck value |", "|---|---|---|"]
    M += [f"| {r['eps_y_nm']:g} | {r['D_y']:.3f} | {'' if r['D_y_deck_value'] is None else f'{r['D_y_deck_value']:.3f}'} |" for r in dp["rows"]]
    ib = o["ipc_background_table"]
    M += ["", "## 8. IPC background (per bunch crossing)", "", "| dataset | ε_y [nm] | n_y | n_m | produced | reaching | seeds |", "|---|---|---|---|---|---|---|"]
    M += [f"| {r['dataset']} | {r['eps_y_nm']:g} | {r['grid']['n_y']} | {r['grid']['n_m']:,} | {r['produced_mean']:,.0f} ± {r['produced_std']['value']:,.0f} | {r['reaching_mean']:.1f} ± {r['reaching_std']['value']:.1f} | {r['n_seeds']} |" for r in ib["rows"]]
    M += ["", "Seeds per row: " + "; ".join(f"{r['dataset']} {r['eps_y_nm']:g} nm: {', '.join(map(str, r['seeds']))}" for r in ib["rows"]), ""]
    db = o["derived_background"]
    M += ["## 9. Derived background quantities (tuned BIB runs, 2 nm vs 20 nm)", "",
          f"- Produced particles, ratio 2 nm / 20 nm: {db['produced_ratio_2nm_over_20nm']['value']:.2f} ± {db['produced_ratio_2nm_over_20nm']['standard_error']['value']:.2f} (SE)",
          f"- Reaching particles, ratio 2 nm / 20 nm: {db['reaching_ratio_2nm_over_20nm']['value']:.2f} ± {db['reaching_ratio_2nm_over_20nm']['standard_error']['value']:.2f} (SE)",
          f"- **BIB-run** luminosity ratio 2 nm / 20 nm (not the conservative-run gain of §2): {db['BIB_RUN_luminosity_ratio_2nm_over_20nm']['value']:.3f} ± {db['BIB_RUN_luminosity_ratio_2nm_over_20nm']['standard_error']['value']:.3f} (SE)",
          f"- Produced per luminosity, reduction: {db['produced_per_luminosity_reduction_pct']['value']:.1f} ± {db['produced_per_luminosity_reduction_pct']['uncertainty']['value']:.1f} %",
          f"- Reaching per luminosity, reduction: {db['reaching_per_luminosity_reduction_pct']['value']:.1f} ± {db['reaching_per_luminosity_reduction_pct']['uncertainty']['value']:.1f} %", ""]
    nx = o["nx_convergence_20nm"]
    M += ["## 10. n_x convergence at ε_y = 20 nm", "", "| n_x | L converged | seeds | converged |", "|---|---|---|---|"]
    M += [f"| {r['n_x']} | {r['L_converged_1e34']:.4f} | {r['n_seeds']} | {'yes' if r['converged'] else 'no'} |" for r in nx["rows"]]
    M += ["", f"Agreement between converged n_x: {100*nx['agreement_between_converged_n_x']['max_pairwise_relative_difference']:.2f} %", ""]
    return "\n".join(M)


if __name__ == "__main__":
    try:
        main()
    except ReproductionError as exc:
        sys.exit(f"reproduce_paper_numbers.py: FAILED: {exc}")
