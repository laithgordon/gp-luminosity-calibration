"""Raw GUINEA-PIG++ .ref parser, used by export_for_wx.py.

Reads ONLY the raw .ref text. No cached CSV value is ever used as a number;
the cached CSVs are consulted for file paths and seed membership only.
"""
from __future__ import annotations
import os, re

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import RAW_ROOT  # noqa: E402  (set GP_RAW_ROOT, see paths.py)

BASE = str(RAW_ROOT)
GAMMA = 125e9 / 0.510998950e6   # C3-250 beam Lorentz factor

def refpath(source_dir: str, filename: str) -> str:
    return os.path.join(BASE, source_dir, "C3_250", filename)

_num = r"([-+0-9.eEnN][^\s;]*)"

def _f(txt, pat, cast=float):
    m = re.search(pat, txt)
    if not m:
        return None
    try:
        return cast(m.group(1))
    except ValueError:
        return None

def parse_ref(path: str) -> dict:
    """Return every field the export needs, straight from the .ref text."""
    with open(path, "r", errors="replace") as fh:
        txt = fh.read()
    st = os.stat(path)
    d = {
        "path": path,
        "mtime": st.st_mtime,
        "size": st.st_size,
        # --- luminosity, raw per-crossing m^-2 as GP++ writes it ---
        "lumi_ee_m2":   _f(txt, r"\blumi_ee=" + _num + r";"),
        "lumi_ee_high_m2": _f(txt, r"\blumi_ee_high=" + _num + r";"),
        "lumi_fine_m2": _f(txt, r"\blumi_fine\s*=\s*" + _num + r"\s*m\*\*"),
        # --- normalisation echoed by the run itself ---
        "f_rep": _f(txt, r"\bf_rep\s*=\s*" + _num),
        "n_b":   _f(txt, r"\bn_b\s*=\s*" + _num),
        # --- physics switches (A3) ---
        "do_photons_1": _f(txt, r"\bdo_photons\.1\s*=\s*" + _num, int),
        "do_photons_2": _f(txt, r"\bdo_photons\.2\s*=\s*" + _num, int),
        "do_pairs":     _f(txt, r"\bdo_pairs\s*=\s*" + _num, int),
        "do_coherent":  _f(txt, r"\bdo_coherent\s*=\s*" + _num, int),
        "do_trident":   _f(txt, r"\bdo_trident\s*=\s*" + _num, int),
        "do_lumi":      _f(txt, r"\bdo_lumi\s*=\s*" + _num, int),
        "do_isr":       _f(txt, r"\bdo_isr\s*=\s*" + _num, int),
        "do_espread":   _f(txt, r"\bdo_espread\s*=\s*" + _num, int),
        "integration_method": _f(txt, r"\bintegration_method\s*=\s*" + _num, int),
        "force_symmetric":    _f(txt, r"\bforce_symmetric\s*=\s*" + _num, int),
        "charge_sign":  _f(txt, r"\bcharge_sign\s*=\s*" + _num),
        # --- grid (A3/B5/B6) ---
        "n_x": _f(txt, r"\bn_x\s*=\s*" + _num, int),
        "n_y": _f(txt, r"\bn_y\s*=\s*" + _num, int),
        "n_z": _f(txt, r"\bn_z\s*=\s*" + _num, int),
        "n_t": _f(txt, r"\bn_t\s*=\s*" + _num, int),
        "n_m_1": _f(txt, r"\bn_m\.1\s*=\s*" + _num, int),
        "n_m_2": _f(txt, r"\bn_m\.2\s*=\s*" + _num, int),
        # --- grid extent, absolute, as written ---
        "cut_x_nm": _f(txt, r"\bcut_x\s*=\s*" + _num + r"\s*nm"),
        "cut_y_nm": _f(txt, r"\bcut_y\s*=\s*" + _num + r"\s*nm"),
        "cut_z_um": _f(txt, r"\bcut_z\s*=\s*" + _num + r"\s*micrometer"),
        # --- beam params echoed in the deck (GP++ prints these directly) ---
        "emitt_x": _f(txt, r"emitt_x\s*:\s*" + _num),
        "emitt_y": _f(txt, r"emitt_y\s*:\s*" + _num),
        "beta_x_um":  _f(txt, r"beta_x\s*:\s*" + _num),
        "beta_y_um":  _f(txt, r"beta_y\s*:\s*" + _num),
        "sigma_x_nm": _f(txt, r"sigma_x\s*:\s*" + _num),
        "sigma_y_nm": _f(txt, r"sigma_y\s*:\s*" + _num),
        "sigma_z_um": _f(txt, r"sigma_z\s*:\s*" + _num),
        "energy_GeV": _f(txt, r"energy\s*:\s*" + _num),
        "particles":  _f(txt, r"particles\s*:\s*" + _num),
        "n_macro_tracked": _f(txt, r"number of tracked macroparticles\s*:\s*" + _num, int),
        # --- beam offsets / waist, as the run's own deck set them (round 2) ---
        "offset_x_nm": _f(txt, r"offset_x\s*:\s*" + _num),
        "offset_y_nm": _f(txt, r"offset_y\s*:\s*" + _num),
        "offset_z_um": _f(txt, r"offset_z\s*:\s*" + _num),
        "waist_x_um":  _f(txt, r"waist_x\s*:\s*" + _num),
        "waist_y_um":  _f(txt, r"waist_y\s*:\s*" + _num),
        # --- grid occupancy ---
        "out_1": _f(txt, r"\bout\.1=" + _num, int),
        "out_2": _f(txt, r"\bout\.2=" + _num, int),
        "miss_1": _f(txt, r"beam 1 : miss\s*=\s*" + _num),
        "miss_2": _f(txt, r"beam 2 : miss\s*=\s*" + _num),
    }
    d["do_photons"] = (d["do_photons_1"], d["do_photons_2"])
    # unit conversion, applied EXACTLY ONCE, using this run's own f_rep / n_b
    for k in ("lumi_ee", "lumi_fine"):
        raw = d[k + "_m2"]
        d[k + "_cm2s"] = (None if raw is None or d["f_rep"] is None or d["n_b"] is None
                          else raw * 1e-4 * d["n_b"] * d["f_rep"])
    # cut multipliers, from the sigmas GP++ itself printed
    for ax in ("x", "y"):
        sig, cut = d["sigma_%s_nm" % ax], d["cut_%s_nm" % ax]
        d["cut_%s_mult" % ax] = (cut / sig) if (sig and cut) else None
    d["cut_z_mult"] = ((d["cut_z_um"] / d["sigma_z_um"])
                       if (d["cut_z_um"] and d["sigma_z_um"]) else None)
    return d
