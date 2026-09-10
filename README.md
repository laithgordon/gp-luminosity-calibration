# GUINEA-PIG++ Luminosity Calibration Analysis

Analysis code, datasets and paper figures for the GUINEA-PIG++ beam–beam
luminosity calibration study (C³-250 collider configuration). The study
calibrates the macroparticle count and grid resolution required for converged
luminosity as a function of vertical beam emittance, and quantifies the
beam-induced background (BIB) at the calibrated settings.

## Contents

| Path | Description |
|------|-------------|
| `GP_CALIBRATION.ipynb` | Main analysis notebook. Produces every figure in `plots/` and the `gp_*_export.csv` datasets. |
| `extract_lumi_data.py` | Scans GUINEA-PIG++ `.ref` output files and writes `data/lumi_extracted.csv`, the root luminosity table every figure is built from. |
| `rebuild_lumi_nominal_vs_calibrated.py` | Re-sources the luminosity columns of `data/lumi_nominal_vs_calibrated.csv` from `data/lumi_extracted.csv`. |
| `extract_highemit_extension.py` | Builds `data/lumi_tuned_vs_frozen_highemit.csv` (tuned vs frozen parameter sets, ε_y = 20–100 nm). |
| `bib_nominal_vs_calibrated.py` | Per-channel BIB pair yield and detector reach, nominal vs calibrated grid → `data/bib_nominal_vs_calibrated.csv`. |
| `bib_stats_table.py`, `bib_stats_render.py` | Per-channel BIB statistics (count, reach fraction, p_T / energy / angle) → `data/bib_stats_cache.csv` (per seed) and `data/bib_stats_summary.csv` (mean ± STD). |
| `pairs_reachability.py`, `reachability_analysis.py` | Detector-reach test used by the BIB scripts (SiD-o2-v04: B = 5 T, r_det = 14 mm, z_max = 76 mm). `reachability_analysis.py` is from [dntounis/Beam_Beam_Backgrounds](https://github.com/dntounis/Beam_Beam_Backgrounds). |
| `inputs/` | GUINEA-PIG++ input file with the nominal C³-250 beam and grid parameters, a template SLURM submit script, and an example scan file. See `inputs/README.md`. |
| `plots/` | Final rendered figures, PDF + PNG. |
| `data/` | Processed datasets (see below). |

## Figures (`plots/`)

| File | Content |
|------|---------|
| `convergence_scan` | Sweep each grid parameter (`n_m`, `n_z`, `n_x`, `n_y`) about the starting operating point and locate luminosity convergence. |
| `convergence_scan_refined` | Re-sweep each parameter about the converged region. |
| `L_vs_nm` | `L` vs macroparticle count at seven ε_y; Richardson fits with the ±5 % crossing `n_m^req` marked. |
| `nm_req_calibration` | `n_m^req/(n_x n_y n_z)` vs disruption `D_y`: the `n_m` law `C_m·D_y^(s−q_p)` and the conservative locus. |
| `L_vs_ny` | `L` vs vertical resolution along the `n_m` locus, with `n_y^req` marked. |
| `ny_req_calibration` | `n_y^req` vs `D_y`: free power-law fit against the derived exponent `q_n`. |
| `kappa_vs_Dy` | `κ`, the cells per pinched vertical σ, vs `D_y` — the test of the pinch model. |
| `L_vs_n_by_emittance_6panel` | `L` vs `n_x` and `n_z` at three emittances, Richardson fits; held-out points marked. |
| `fig4_pseudoflat_instability_all_params` | Coefficient of variation of `L` across grid settings vs ε_y. |
| `variance_3d_lumi_mesh` | Luminosity variance across the (`n_x`, `n_y`, ε_y) grid. |
| `lumi_bib_tradeoff_PRL` | Luminosity gain vs BIB cost. |

## Datasets (`data/`)

| File | Content |
|------|---------|
| `lumi_extracted.csv` | One row per simulation run: run parameters parsed from the filename, luminosity metrics from the `.ref`. **Units:** `lumi_ee`/`lumi_fine` are rates in cm⁻² s⁻¹, converted per row from GP++'s per-crossing m⁻² value with the `f_rep` and `n_b` echoed in the same `.ref` (× 10⁻⁴·n_b·f_rep = 1.596 for C³-250); the raw per-crossing values are in `lumi_ee_m2`/`lumi_fine_m2`. Do not apply the factor again downstream. |
| `lumi_nominal_vs_calibrated.csv` | Per-seed luminosity at the nominal and calibrated grids, ε_y = 0.5–20 nm. |
| `lumi_tuned_vs_frozen_highemit.csv` | Per-seed luminosity, tuned vs frozen parameter sets, ε_y = 20–100 nm. |
| `lumi_bib_per_seed_cache.csv` | Per-seed luminosity and BIB reach feeding `lumi_bib_tradeoff_PRL`. |
| `bib_nominal_vs_calibrated.csv` | Per-seed, per-channel (BW / BH / LL) pair counts and detector reach, nominal vs calibrated grid. |
| `bib_stats_cache.csv`, `bib_stats_summary.csv` | Per-channel BIB statistics at ε_y = 1, 8, 20 nm: per seed, and mean ± STD. |
| `gp_kappa_export.csv` | `κ` per ε_y with `D_y`, `R(D_y)`, `n_y^req`, `c_y` — for cross-code comparison at fixed `κ`. |
| `gp_calibration_export.csv` | The `n_y` and `n_m` calibration families in one long table. |
| `gp_constants_export.csv` | Every fitted and derived constant with its error. |

## Regenerating

The raw simulation output (`.ref`, `pairs.dat`, `output*/`) is **not** in this
repository — it is large and produced by GUINEA-PIG++ runs. To rebuild from raw
output, point `SOURCE_DIRS` at the top of `extract_lumi_data.py` at your
GUINEA-PIG++ output directories, then:

```bash
python3 extract_lumi_data.py          # -> data/lumi_extracted.csv
jupyter nbconvert --execute GP_CALIBRATION.ipynb --to notebook --inplace
                                      # -> plots/*, data/gp_*_export.csv
python3 bib_stats_table.py && python3 bib_stats_render.py
python3 bib_nominal_vs_calibrated.py  # BIB tables (need pairs.dat dumps)
```

The BIB scripts and `extract_highemit_extension.py` locate the raw output
relative to their own directory's parent; adjust the `REPO` path at the top of
each if the layout differs.

## Requirements

Python 3 with `numpy`, `pandas`, `matplotlib`, `scipy`; Jupyter for the notebook.
