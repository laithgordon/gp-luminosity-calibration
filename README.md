# GUINEA-PIG++ Luminosity Tuning Analysis

Analysis code, data and figures behind the GUINEA-PIG++ results of the C³-250
beam–beam study: the grid resolution and macroparticle count required for a
converged luminosity as a function of vertical emittance, the conservative
settings derived from them, the beam-induced background (BIB) at those settings,
and the luminosity exported for the WarpX code comparison.

## Paper numbers

`reproduce_paper_numbers.py` regenerates every GUINEA-PIG++ number the paper
reports, from the committed data alone: no raw output, no network, no random
numbers. Run it from a fresh clone:

```bash
python3 reproduce_paper_numbers.py   # -> paper_numbers.json, paper_numbers.md
```

`paper_numbers.md` holds the tables rounded as the manuscript rounds them;
`paper_numbers.json` holds the same quantities at full precision, with the
definition of every uncertainty and the data file and selection rule behind
every quantity. Its output corresponds to the numbers in the paper, so a
reviewer can compare the two directly. It stops with an error, rather than
writing a partial file, if an input is missing or a recomputed fit disagrees
with the notebook's exports. What it covers is listed in `REPRODUCTION.md`.

## Reproducing the results

Python 3 with `numpy`, `pandas`, `matplotlib` and `scipy`, plus Jupyter
(`nbconvert`) for the notebook. The published figures were produced with Python
3.13, numpy 2.3, pandas 2.3, matplotlib 3.10.5 and scipy 1.16. GUINEA-PIG++ 1.2.2 produced the
runs.

The raw GUINEA-PIG++ output (`.ref` files and `pairs.dat` dumps, several GB) is
not in git. Steps marked **raw** read it; point `GP_RAW_ROOT` at the directory
that holds `output/`, `output_nm/` and `output_no_pairs/` (each with a `C3_250/`
subdirectory). `paths.py` is the only place that path is resolved; if the
variable is unset, the parent of this repository is used. Every other step runs
from the committed data alone. Run from the repository root, in this order:

```bash
export GP_RAW_ROOT=/path/to/GuineaPig_Feb_2025

# 1. raw   luminosity table read by the notebook (committed snapshot, see below)
python3 extract_lumi_data.py --restrict-to data/lumi_extracted.csv

# 2. raw   BIB statistics per seed, then the tab:bib_yields numbers
python3 bib_stats_table.py --refresh       # -> data/bib_stats_cache.csv
python3 bib_stats_render.py                # -> data/bib_stats_summary.csv, data/bib_stats_table.md

# 3. raw   detector reach per seed for the luminosity-vs-BIB figure
python3 bib_reach_per_seed.py              # -> data/lumi_bib_per_seed_cache.csv

# 4.       every notebook figure and the fitted constants
jupyter nbconvert --to notebook --execute --inplace GP_CALIBRATION.ipynb
                                           # -> plots/*, data/gp_calibration_export.csv,
                                           #    data/gp_constants_export.csv

# 5. raw   GP++ numbers for the WarpX comparison
python3 export_for_wx.py                   # -> data/gp_luminosity_for_wx.csv,
                                           #    data/gp_requirements_for_wx.csv

# 6. raw   the n_x x n_m convergence figure (drop --refresh to plot the committed table)
python3 nx_nm_convergence.py --refresh     # -> data/L_vs_nm_by_nx_20nm.csv, plots/L_vs_nm_by_nx_20nm.*

# 7.       conservative-locus scan rows, checked against the production decks
python3 inputs/make_conservative_decks.py --check
```

Step 5 reads `data/gp_calibration_export.csv` and `data/gp_constants_export.csv`,
so it runs after step 4.

The committed `data/lumi_bib_per_seed_cache.csv` has 593 rows, including runs the
figure does not use. Step 3 rewrites
it with only the 263 runs on the figure's curves. For those 263, `n_reach`, the
only column the figure uses, is identical; `n_total` differs by one particle in
two runs. The figure is unchanged either way.

### The luminosity snapshot

`data/lumi_extracted.csv` is the luminosity table of 2026-09-10 (16,803 runs).
Every notebook figure and every fitted constant was produced from it, and it
regenerates exactly from raw output with `--restrict-to`, as in step 1. The raw
tree has since gained 143 runs: the conservative-locus production runs, the
conservative BIB reruns, and the n_x × n_m convergence runs. Some of them pass
the notebook's ladder selection. A plain `python3 extract_lumi_data.py` includes
them, which changes `L_vs_nm`, `L_vs_ny`, `nm_req_tuning`, `ny_req_tuning` and
`kappa_vs_Dy` and moves the fitted constants (for example q_n,fit 0.402 → 0.378,
C_y,fit 26.2 → 27.9). The published values are those of the snapshot. The
conservative-locus luminosities are read straight from the raw `.ref` files by
`export_for_wx.py`, so they do not depend on the snapshot.

## Figures (`plots/`, PDF and PNG)

**In the manuscript**

| File | Content | Produced by |
|---|---|---|
| `L_vs_nm` | `L` vs `n_m` at seven ε_y with Richardson fits and the ±5 % crossing `n_m^req`. | notebook |
| `nm_req_tuning` | `n_m^req/(n_x n_y n_z)` vs disruption `D_y`: the `n_m` law `C_m·D_y^(s−q_p)` and the conservative locus. | notebook |
| `L_vs_ny` | `L` vs `n_y` along the `n_m` locus, with `n_y^req`. | notebook |
| `ny_req_tuning` | `n_y^req` vs `D_y`: free power-law fit against the derived exponent `q_n`. | notebook |
| `L_vs_n_by_emittance_6panel` | `L` vs `n_x` and `n_z` at three emittances with Richardson fits. | notebook |
| `fig4_pseudoflat_instability_all_params` | Coefficient of variation of `L` across grid settings vs ε_y. | notebook |
| `lumi_bib_tradeoff_PRL` | Luminosity gain vs BIB cost. | notebook |

**Repository only**

| File | Content | Produced by |
|---|---|---|
| `convergence_scan` | Sweep of each grid parameter about the starting operating point. | notebook |
| `convergence_scan_refined` | Re-sweep about the converged region. | notebook |
| `kappa_vs_Dy` | `κ`, cells per pinched vertical σ, vs `D_y`, GUINEA-PIG++ points only. Not the manuscript's κ figure: that one shows both simulators and is produced in the WarpX repository, under the same filename. | notebook |
| `variance_3d_lumi_mesh` | Luminosity variance across (`n_x`, `n_y`, ε_y). | notebook |
| `L_vs_nm_by_nx_20nm` | `L` vs `n_m` at ε_y = 20 nm for n_x = 256, 512, 1024, 2048 (seed 100805). | `nx_nm_convergence.py` |

The manuscript's GUINEA-PIG++ discretization figure
(`guinea_pig_discretization_final.pdf`) is not produced by this repository.

## Published numbers

| Result | Source | Script |
|---|---|---|
| `tab:ny_recommendations` (`n_y^rec`, `n_m^rec` per ε_y) | `constants.py`, `data/D_y_table.json` | `inputs/make_conservative_decks.py` |
| `tab:bib_yields` | `data/bib_stats_cache.csv` | `bib_stats_render.py` |
| Nominal, conservative and 40–100 nm luminosities of the code comparison | raw `.ref` | `export_for_wx.py` |
| Fitted constants `C_y,fit`, `q_n,fit`, `C_m,fit`, `s_fit` | `data/lumi_extracted.csv`, with `q_p` from `constants.py` | notebook |
| Requirements with uncertainties (`sigma_log`), exported to WarpX; κ for both simulators is computed in the WarpX repository from this export | `data/gp_calibration_export.csv` | `export_for_wx.py` |

## Which rule set each configuration

The conservative-locus constants are frozen in `constants.py`: they were fixed
when the conservative runs were submitted and are never recomputed from the
current fit.

| Configuration | Used for | `n_y` | `n_m` |
|---|---|---|---|
| Nominal, 512×512×25 | nominal rows of the code comparison and of `tab:bib_yields` | fixed 512 | fixed 100 000 |
| Pre-campaign tuned grid, n_x = 512, n_z = 64: `n_y` = 256, 128, 128, 128 at 2, 4, 8, 20 nm | tuned rows of `tab:bib_yields` | `34.70·D_y^0.312`, rounded up to a power of two | the `n_m` locus at that `n_y` at 2 and 4 nm (2 235 661; 353 910); 250 000 and 80 000 at 8 and 20 nm, above the locus |
| Frozen extension, 40–100 nm, 512×128×64 | `frozen_extension` block of `gp_luminosity_for_wx.csv` | same pre-campaign rule, saturated at 128 | 50 000 floor |
| Conservative locus, 1–20 nm, n_x = 512, n_z = 64 | `conservative` block of `gp_luminosity_for_wx.csv`; `tab:ny_recommendations` | `50·D_y^0.402`, rounded up (512, 512, 256, 256, 256, 256, 256) | `⌈2.305×10⁻⁷·D_y^3.300·512·n_y·64⌉`; at 1 nm the existing ladder rung 14 100 000 (ratio 1.0003) |

The pre-campaign and conservative rules round the same way and differ only in
normalisation and exponent: the conservative `n_y` is exactly twice the
pre-campaign `n_y` at every ε_y.

## Contents

| Path | Role |
|---|---|
| `paths.py` | Locates the raw output (`GP_RAW_ROOT`) and the repository data. |
| `constants.py` | Frozen conservative-locus constants and locus functions, and `q_p` (envelope-integration exponent, shared with the WarpX repository) with `q_n^pred = q_p + 1/4`. |
| `extract_lumi_data.py` | Raw `.ref` → `data/lumi_extracted.csv`. |
| `GP_CALIBRATION.ipynb` | Tuning ladders, fits, and every figure except `L_vs_nm_by_nx_20nm`. |
| `bib_stats_table.py`, `bib_stats_render.py` | BIB statistics from `pairs.dat` (`find_dumps` refuses ambiguous dumps) → `tab:bib_yields`. |
| `bib_reach_per_seed.py` | Detector reach per seed for `lumi_bib_tradeoff_PRL`. |
| `pairs_reachability.py`, `reachability_analysis.py` | Detector-reach test (SiD-o2-v04: B = 5 T, r_det = 14 mm, z_max = 76 mm). `reachability_analysis.py` is from [dntounis/Beam_Beam_Backgrounds](https://github.com/dntounis/Beam_Beam_Backgrounds). |
| `refparse.py`, `export_for_wx.py` | Raw `.ref` parser and the WarpX export. |
| `nx_nm_convergence.py` | The n_x × n_m convergence figure. |
| `inputs/` | GP++ input decks, templates, the conservative-deck generator and every run configuration; see `inputs/README.md`. |

| Data file | Read by |
|---|---|
| `data/lumi_extracted.csv` | notebook, `bib_reach_per_seed.py` |
| `data/lumi_bib_per_seed_cache.csv` | notebook |
| `data/bib_stats_cache.csv` | `bib_stats_render.py` |
| `data/gp_calibration_export.csv`, `data/gp_constants_export.csv` | `export_for_wx.py` |
| `data/D_y_table.json` | `export_for_wx.py`, `inputs/make_conservative_decks.py` |
| `data/L_vs_nm_by_nx_20nm.csv` | `nx_nm_convergence.py` |
| `data/gp_luminosity_for_wx.csv`, `data/gp_requirements_for_wx.csv` | the WarpX analysis |

**Units.** `lumi_ee` and `lumi_fine` in `data/lumi_extracted.csv` are rates in
cm⁻² s⁻¹, converted per row from GUINEA-PIG++'s per-crossing m⁻² value with the
`f_rep` and `n_b` echoed in the same `.ref` (× 10⁻⁴·n_b·f_rep = 1.596 for
C³-250). The raw per-crossing values are in `lumi_ee_m2` and `lumi_fine_m2`. Do
not apply the factor again.
