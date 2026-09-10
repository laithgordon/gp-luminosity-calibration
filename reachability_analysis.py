#!/usr/bin/env python3
"""Utility for visualising the reachability boundary in the pT–theta plane."""

import argparse
import os
from typing import Dict, Optional, Sequence

import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import FuncFormatter

try:
    import mplhep as hep

    HAS_MPLHEP = True
except Exception:
    HAS_MPLHEP = False


def find_theta_for_pt(pt: float, B0: float = 5.0, *, q: float = 0.3,
                      r_det: float = 14.0, z_max: float = 76.0,
                      theta_upper: Optional[float] = None) -> float:
    if pt <= 0.0:
        return np.nan

    denom = q * B0
    if denom == 0.0:
        return np.nan

    R_m = pt / denom  # radius in metres when pt is in GeV/c
    if not np.isfinite(R_m) or R_m <= 0.0:
        return np.nan

    R = R_m * 1000.0  # convert to mm to match detector geometry inputs
    if not np.isfinite(R) or R <= 0.0:
        return np.nan

    if 2.0 * R < r_det:
        return np.nan

    if z_max <= 0.0:
        return np.nan

    ratio = r_det / (2.0 * R)
    if ratio > 1.0 + 1e-9:
        return np.nan
    ratio = min(max(ratio, 0.0), 1.0)
    alpha = np.arcsin(ratio)  # Minimal theta occurs when 2R sin(alpha) = r_det
    tan_theta = (2.0 * R * alpha) / z_max
    if not np.isfinite(tan_theta) or tan_theta <= 0.0:
        return np.nan

    theta = np.arctan(tan_theta)

    if theta_upper is not None and theta > theta_upper:
        return np.nan

    return float(theta)


def compute_boundary(pt_values: np.ndarray, q: float = 0.3, B0: float = 5.0,
                      r_det: float = 14.0, z_max: float = 76.0,
                      theta_upper: Optional[float] = None) -> np.ndarray:
    """Return array of minimal theta solutions for each pT in ``pt_values``."""

    pt_arr = np.asarray(pt_values, dtype=float)
    thetas = [
        find_theta_for_pt(
            float(pt),
            B0=B0,
            q=q,
            r_det=r_det,
            z_max=z_max,
            theta_upper=theta_upper,
        )
        for pt in pt_arr
    ]
    return np.asarray(thetas, dtype=float)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Plot the reachability boundary by solving for the minimal theta at each pT.'
    )
    parser.add_argument('--pt-min', type=float, help='Minimum pT in GeV/c (default enforces 2R >= r_det).')
    parser.add_argument('--pt-max', type=float, default=2.0, help='Maximum pT in GeV/c.')
    parser.add_argument('--pt-samples', type=int, default=500, help='Number of pT samples (log spaced).')
    parser.add_argument('--linear-pt', action='store_true', help='Use linear spacing in pT instead of log.')
    parser.add_argument('--charge', type=float, default=0.3, help='Effective charge factor q (GeV*T*mm).')
    parser.add_argument('--mag-field', type=float, default=5.0, help='Magnetic field B0 in Tesla.')
    parser.add_argument('--detector-radius', type=float, default=14.0, help='Detector radius in millimetres.')
    parser.add_argument('--z-max', type=float, default=76.0, help='Maximum z-extent in millimetres.')
    parser.add_argument('--theta-upper', type=float, help='Optional upper bound on theta in radians for solver.')
    parser.add_argument('--out', default='reachability_boundary.png', help='Primary output plot path (png/jpg/etc).')
    parser.add_argument('--out-pdf', help='Optional PDF output path (defaults to primary path with .pdf).')
    parser.add_argument('--plot-theta-min', type=float, help='Override lower theta limit for plot (rad).')
    parser.add_argument('--plot-theta-max', type=float, help='Override upper theta limit for plot (rad).')
    parser.add_argument('--plot-pt-min', type=float, help='Override lower pT limit for plot (GeV/c).')
    parser.add_argument('--plot-pt-max', type=float, help='Override upper pT limit for plot (GeV/c).')
    return parser.parse_args()


def _publication_rcparams() -> Dict[str, float]:
    """Return rcParams tweaks for a publication-ready appearance."""

    return {
        "font.size": 21,
        "axes.titlesize": 23,
        "axes.labelsize": 23,
        "legend.fontsize": 20,
        "xtick.labelsize": 21,
        "ytick.labelsize": 21,
        "axes.linewidth": 1.2,
    }


def _format_sigfig(value: float, sig: int = 2) -> str:
    """Format ``value`` with ``sig`` significant digits without scientific notation."""

    if not np.isfinite(value):
        return str(value)
    if value == 0:
        return "0"

    abs_value = abs(value)
    exponent = int(np.floor(np.log10(abs_value)))
    decimal_places = sig - 1 - exponent
    rounded = round(value, decimal_places)
    if decimal_places > 0:
        formatted = f"{rounded:.{decimal_places}f}"
        formatted = formatted.rstrip('0').rstrip('.')
        return formatted if formatted else "0"
    return f"{rounded:.0f}"


def plot_boundary(theta_values: np.ndarray, pt_values: np.ndarray, out_path: str,
                  B0: float, rect_info: Dict[str, float],
                  extra_paths: Optional[Sequence[str]] = None,
                  plot_theta_min: Optional[float] = None,
                  plot_theta_max: Optional[float] = None,
                  plot_pt_min: Optional[float] = None,
                  plot_pt_max: Optional[float] = None) -> None:
    if HAS_MPLHEP:
        hep.style.use("CMS")

    plt.rcParams.update(_publication_rcparams())

    fig, ax = plt.subplots(figsize=(10, 7))
    #ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\theta$ [rad]')
    ax.set_ylabel(r'$p_T$ [MeV/$c$]')
    ax.grid(which='both', alpha=0.25, linestyle='--', linewidth=0.8)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda val, _: _format_sigfig(val, sig=2)))

    theta_plot = np.asarray(theta_values, dtype=float)
    pt_plot = np.asarray(pt_values, dtype=float) * 1e3  # convert from GeV to MeV for plotting

    theta_min = float(theta_plot.min())
    theta_max = float(theta_plot.max())
    pt_min = float(pt_plot.min())
    pt_max = float(pt_plot.max())

    theta_lo_default = max(theta_min * 0.9, 1e-6)
    theta_hi_default = theta_max * 1.05
    if plot_theta_min is not None:
        theta_lo = max(float(plot_theta_min), 1e-6)
    else:
        theta_lo = theta_lo_default
    if plot_theta_max is not None:
        theta_hi = max(float(plot_theta_max), theta_lo * (1 + 1e-6))
    else:
        theta_hi = theta_hi_default
    if theta_hi <= theta_lo:
        theta_hi = theta_lo * 1.01

    pt_lo_default = max(pt_min * 0.8, 1e-6)
    pt_hi_default = pt_max * 1.15
    plot_pt_min_mev = plot_pt_min * 1e3 if plot_pt_min is not None else None
    plot_pt_max_mev = plot_pt_max * 1e3 if plot_pt_max is not None else None
    if plot_pt_min_mev is not None:
        pt_lo = max(float(plot_pt_min_mev), 1e-6)
    else:
        pt_lo = pt_lo_default
    if plot_pt_max_mev is not None:
        pt_hi = max(float(plot_pt_max_mev), pt_lo * (1 + 1e-6))
    else:
        pt_hi = pt_hi_default
    if pt_hi <= pt_lo:
        pt_hi = pt_lo * 1.01

    if theta_plot[-1] < theta_hi:
        pt_floor = float(pt_plot.min())
        theta_plot = np.append(theta_plot, theta_hi)
        pt_plot = np.append(pt_plot, pt_floor)

    ax.set_xlim(theta_lo, theta_hi)
    ax.set_ylim(pt_lo, pt_hi)

    boundary_line, = ax.plot(
        theta_plot,
        pt_plot,
        color='#1f77b4',
        linewidth=2.5,
        label='Reachability boundary'
    )

    ax.fill_between(
        theta_plot,
        pt_plot,
        np.full_like(pt_plot, pt_hi),
        color='#1f77b4',
        alpha=0.09,
        edgecolor='none',
        zorder=1
    )

    asym_theta = rect_info['theta_asym']
    asym_pt = rect_info['pt_asym'] * 1e3  # convert from GeV to MeV

    if asym_theta < theta_hi and asym_pt < pt_hi:
        ax.plot(
            [asym_theta, asym_theta],
            [asym_pt, pt_hi],
            color='#d62728',
            linestyle='--',
            linewidth=1.8,
        )
        ax.plot(
            [asym_theta, theta_hi],
            [asym_pt, asym_pt],
            color='#d62728',
            linestyle='-.',
            linewidth=1.8,
        )

        rect = patches.Rectangle(
            (asym_theta, asym_pt),
            theta_hi - asym_theta,
            pt_hi - asym_pt,
            facecolor='#d62728',
            alpha=0.08,
            edgecolor='none',
            zorder=0
        )
        ax.add_patch(rect)

    ax.text(
        0.10,
        0.93,
        rf'$B = {B0:.1f}\,\mathrm{{T}}$',
        transform=ax.transAxes,
        fontsize=21,
        color='#1f77b4',
        bbox=dict(facecolor='white', alpha=0.10, edgecolor='none', boxstyle='round,pad=0.3')
    )

    theta_display = _format_sigfig(asym_theta, sig=2)
    pt_display = _format_sigfig(asym_pt, sig=2)
    asym_text = (
        r'$\theta_{\rm asym} = ' + f'{theta_display}' + r'\,\mathrm{{rad}}$' + '\n'
        r'$p_{T,\rm asym} = ' + f'{pt_display}' + r'\,\mathrm{{MeV}/c}$'
    )
    ax.text(
        0.62,
        0.12,
        asym_text,
        transform=ax.transAxes,
        fontsize=21,
        color='#d62728',
        bbox=dict(facecolor='white', alpha=0.10, edgecolor='none', boxstyle='round,pad=0.3')
    )

    handles = [boundary_line, patches.Patch(facecolor='#d62728', alpha=0.18, edgecolor='none')]
    labels = ['Reachability boundary', 'Rectangular asymptote']
    ax.legend(handles, labels, loc='upper right', frameon=True, fancybox=True)
    ax.tick_params(which='both', direction='out', length=5)
    ax.set_title(r'Detector reach in $p_T$–$\theta$ space', loc='left')
    ax.text(
        0.99,
        1.02,
        'SiD_o2_v04',
        transform=ax.transAxes,
        ha='right',
        va='bottom',
        fontsize=17,
        color='#333333'
    )

    fig.tight_layout()

    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    fig.savefig(out_path, dpi=450, bbox_inches='tight')
    print(f"Saved reachability boundary plot to {out_path}")

    if extra_paths:
        for extra in extra_paths:
            if not extra:
                continue
            if os.path.abspath(extra) == os.path.abspath(out_path):
                continue
            extra_dir = os.path.dirname(extra)
            if extra_dir:
                os.makedirs(extra_dir, exist_ok=True)
            is_raster = extra.lower().endswith(('.png', '.jpg', '.jpeg', '.tif', '.tiff'))
            fig.savefig(
                extra,
                dpi=450 if is_raster else None,
                bbox_inches='tight'
            )
            print(f"Saved additional plot to {extra}")

    plt.close(fig)


def main():
    args = parse_args()

    pt_min_condition = (args.charge * args.mag_field * args.detector_radius) / 2000.0
    provided_pt_min = args.pt_min
    pt_min = max(provided_pt_min if provided_pt_min is not None else pt_min_condition, pt_min_condition)

    if pt_min <= 0.0:
        raise ValueError('Computed pt_min is non-positive; check charge, field, and detector radius inputs.')

    if provided_pt_min is not None and provided_pt_min < pt_min_condition:
        print(
            f"Note: provided pt_min={provided_pt_min:.6g} GeV/c is below the 2R>=r_det threshold; using {pt_min:.6g} GeV/c instead."
        )

    pt_max = max(args.pt_max, pt_min * 1.01)
    samples = max(args.pt_samples, 5)

    if args.linear_pt:
        pt_values = np.linspace(pt_min, pt_max, samples)
    else:
        pt_values = np.logspace(np.log10(pt_min), np.log10(pt_max), samples)

    if pt_values.size > 0:
        pt_values[0] = pt_min
        pt_values[-1] = pt_max

    theta_values = compute_boundary(
        pt_values,
        q=args.charge,
        B0=args.mag_field,
        r_det=args.detector_radius,
        z_max=args.z_max,
        theta_upper=args.theta_upper,
    )

    mask = (
        np.isfinite(theta_values)
        & np.isfinite(pt_values)
        & (theta_values > 0)
        & (pt_values > 0)
    )
    if not np.any(mask):
        raise RuntimeError('No valid boundary points to plot. Check the configuration values.')

    invalid = int(np.size(pt_values) - np.count_nonzero(mask))
    if invalid > 0:
        print(f"Warning: discarded {invalid} pT samples without valid reachability solutions.")

    theta_valid = theta_values[mask]
    pt_valid = pt_values[mask]

    order = np.argsort(theta_valid)
    theta_plot = theta_valid[order]
    pt_plot = pt_valid[order]

    theta_asym = float(theta_plot[0])
    pt_asym = float(pt_plot[-1])
    rect_info = {
        'theta_asym': theta_asym,
        'pt_asym': pt_asym,
    }

    extra_outputs = []
    if args.out_pdf:
        extra_outputs.append(args.out_pdf)
    else:
        root, ext = os.path.splitext(args.out)
        if ext.lower() != '.pdf':
            extra_outputs.append(root + '.pdf')

    # Remove duplicates and the primary output path if inadvertently scheduled twice
    normalized_primary = os.path.abspath(args.out)
    unique_extras = []
    seen = set()
    for path in extra_outputs:
        if not path:
            continue
        norm = os.path.abspath(path)
        if norm == normalized_primary or norm in seen:
            continue
        seen.add(norm)
        unique_extras.append(path)

    plot_boundary(
        theta_plot,
        pt_plot,
        args.out,
        args.mag_field,
        rect_info,
        unique_extras,
        plot_theta_min=args.plot_theta_min,
        plot_theta_max=args.plot_theta_max,
        plot_pt_min=args.plot_pt_min,
        plot_pt_max=args.plot_pt_max,
    )

    print(
        f"Rectangular asymptote (legacy ROI): theta >= {theta_asym:.6g} rad, pT >= {pt_asym:.6g} GeV/c"
    )
    print(
        f"Boundary spans theta in [{theta_plot.min():.6g}, {theta_plot.max():.6g}] rad and pT in [{pt_plot.min():.6g}, {pt_plot.max():.6g}] GeV/c"
    )


if __name__ == '__main__':
    main()
