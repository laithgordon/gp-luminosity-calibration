"""
Count GUINEA-PIG++ pair-background particles that reach the detector.

This is a thin companion to `reachability_analysis.py` (the upstream
boundary/plotting tool from dntounis/Beam_Beam_Backgrounds). All boundary
physics is delegated to that module via `find_theta_for_pt` /
`compute_boundary`; this file only:

    * parses a GUINEA-PIG++ pairs.dat file into (pT, theta) arrays, and
    * counts how many particles sit on the "reaches" side of the boundary.

GUINEA-PIG++ pairs.dat column layout (from src/particlesCPP.h PAIR_PARTICLE):
    E[GeV, signed by charge]  vx  vy  vz  x[nm]  y[nm]  z[nm]  process  slice1  slice2
where (vx, vy, vz) is the unit velocity vector. For ultra-relativistic
pairs |v|~c, so p_i = E * v_i to excellent approximation.
"""

from __future__ import annotations

import argparse

import numpy as np

from reachability_analysis import compute_boundary


def load_pairs_dat(path: str) -> dict[str, np.ndarray]:
    """Parse a GUINEA-PIG++ pairs.dat file.

    Returns a dict with arrays: E (GeV, unsigned), charge (+1 / -1),
    px, py, pz (GeV), pt (GeV), theta (rad, [0, pi]),
    x, y, z (nm, as written by GP++), process.
    """
    # Some pair.dat files have ragged rows (occasional truncation). Only the
    # first 8 columns feed the reachability test. pandas.read_csv with the C
    # engine is ~5-10x faster than np.genfromtxt on these multi-million-row
    # text files; `on_bad_lines="skip"` tolerates truncated rows.
    import pandas as _pd
    _df = _pd.read_csv(path, sep=r"\s+", header=None,
                       usecols=range(8), engine="c",
                       on_bad_lines="skip", dtype=float)
    data = _df.to_numpy()
    if data.size == 0:
        empty = np.empty(0)
        return dict(E=empty, charge=empty, px=empty, py=empty, pz=empty,
                    pt=empty, theta=empty, x=empty, y=empty, z=empty,
                    process=np.empty(0, dtype=int))
    if data.ndim == 1:
        data = data[None, :]

    E_signed = data[:, 0]
    vx, vy, vz = data[:, 1], data[:, 2], data[:, 3]
    x, y, z = data[:, 4], data[:, 5], data[:, 6]
    process = data[:, 7].astype(int)

    E = np.abs(E_signed)
    charge = np.where(E_signed >= 0, 1, -1)

    px = E * vx
    py = E * vy
    pz = E * vz
    pt = np.hypot(px, py)

    v_norm = np.sqrt(vx * vx + vy * vy + vz * vz)
    theta = np.arccos(np.clip(vz / v_norm, -1.0, 1.0))

    return dict(E=E, charge=charge, px=px, py=py, pz=pz,
                pt=pt, theta=theta, x=x, y=y, z=z, process=process)


def _build_boundary_interp(pt_particles: np.ndarray,
                           *, mag_field: float, detector_radius: float,
                           z_max: float, charge: float,
                           pt_samples: int = 2000):
    """Return `(pt_grid, theta_min_grid)` from `compute_boundary` covering the
    pT range of the particle sample. Used for fast interpolation so we don't
    call `find_theta_for_pt` once per particle.
    """
    pt_min_condition = (charge * mag_field * detector_radius) / 2000.0  # GeV
    pt_valid = pt_particles[pt_particles > pt_min_condition]
    if pt_valid.size == 0:
        return None, pt_min_condition

    pt_lo = pt_min_condition * (1.0 + 1e-6)
    pt_hi = float(pt_valid.max()) * 1.01
    if pt_hi <= pt_lo:
        pt_hi = pt_lo * 1.01

    pt_grid = np.logspace(np.log10(pt_lo), np.log10(pt_hi), pt_samples)
    theta_grid = compute_boundary(
        pt_grid,
        q=charge,
        B0=mag_field,
        r_det=detector_radius,
        z_max=z_max,
        theta_upper=None,
    )

    good = np.isfinite(theta_grid) & (theta_grid > 0)
    pt_grid = pt_grid[good]
    theta_grid = theta_grid[good]

    # theta_min is monotonically decreasing in pT (higher pT -> easier reach);
    # np.interp needs the x-axis sorted ascending, which log-grid already is.
    return (pt_grid, theta_grid), pt_min_condition


def reaches_detector(pt: np.ndarray, theta: np.ndarray,
                     *, mag_field: float = 5.0,
                     detector_radius: float = 14.0,
                     z_max: float = 76.0,
                     charge: float = 0.3) -> np.ndarray:
    """Vectorized per-particle reachability test.

    Uses `compute_boundary` from `reachability_analysis` as the sole physics
    source: we build a fine (pT, theta_min) curve and check
    `theta_folded >= theta_min(pT)`, where `theta_folded = min(theta, pi-theta)`
    to handle the z-symmetric barrel.
    """
    pt = np.asarray(pt, dtype=float)
    theta = np.asarray(theta, dtype=float)

    interp, pt_min_condition = _build_boundary_interp(
        pt,
        mag_field=mag_field,
        detector_radius=detector_radius,
        z_max=z_max,
        charge=charge,
    )

    mask = np.zeros(pt.shape, dtype=bool)
    if interp is None:
        return mask

    pt_grid, theta_grid = interp
    # Particles below the 2R >= r_det pT threshold can never reach.
    above_pt_threshold = pt > pt_min_condition

    theta_folded = np.minimum(theta, np.pi - theta)

    # np.interp clamps to edge values outside [pt_grid.min, pt_grid.max];
    # that's the physically correct limit (theta_min plateaus at the endpoints).
    theta_min_at_pt = np.interp(pt, pt_grid, theta_grid)

    mask = above_pt_threshold & (theta_folded >= theta_min_at_pt)
    return mask


def count_reaching_particles(path: str,
                              *, mag_field: float = 5.0,
                              detector_radius: float = 14.0,
                              z_max: float = 76.0,
                              charge: float = 0.3,
                              return_mask: bool = False):
    """Count pair particles from `path` that reach the detector barrel.

    Parameters match the CLI of `reachability_analysis.py` exactly so the
    same geometry definition can be shared.

    Returns
    -------
    int
        Number of particles reaching the detector (default), or
        (n_reach, n_total, mask) when `return_mask=True`.
    """
    pairs = load_pairs_dat(path)
    mask = reaches_detector(
        pairs["pt"], pairs["theta"],
        mag_field=mag_field,
        detector_radius=detector_radius,
        z_max=z_max,
        charge=charge,
    )
    n_reach = int(mask.sum())
    n_total = int(mask.size)
    if return_mask:
        return n_reach, n_total, mask
    return n_reach


def _cli() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("pairs_file", help="Path to a GUINEA-PIG++ pairs.dat file.")
    p.add_argument("--mag-field", type=float, default=5.0)
    p.add_argument("--detector-radius", type=float, default=14.0)
    p.add_argument("--z-max", type=float, default=76.0)
    p.add_argument("--charge", type=float, default=0.3,
                   help="helix constant q in R[m] = p_T[GeV/c] / (q * B[T]) for a unit-charge "
                        "track; q = 0.2998 GeV/(c T m). Not a particle charge "
                        "(default: %(default)s)")
    args = p.parse_args()

    n_reach, n_total, _ = count_reaching_particles(
        args.pairs_file,
        mag_field=args.mag_field,
        detector_radius=args.detector_radius,
        z_max=args.z_max,
        charge=args.charge,
        return_mask=True,
    )
    frac = n_reach / n_total if n_total else 0.0
    print(f"{args.pairs_file}")
    print(f"  geometry: B={args.mag_field} T, r_det={args.detector_radius} mm, "
          f"z_max={args.z_max} mm, q={args.charge}")
    print(f"  {n_reach}/{n_total} particles reach detector ({frac:.3%})")


if __name__ == "__main__":
    _cli()
