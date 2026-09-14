# Reproducing the GUINEA-PIG++ numbers in the paper

`reproduce_paper_numbers.py` regenerates every number the paper reports from
GUINEA-PIG++, using only files committed to this repository. It writes
`paper_numbers.json`, at full precision with a definition on every uncertainty
and the data file and selection rule behind every quantity, and
`paper_numbers.md`, with the tables rounded as the manuscript rounds them.

## How to run

```bash
git clone https://github.com/laithgordon/gp-luminosity-calibration.git
cd gp-luminosity-calibration
python3 reproduce_paper_numbers.py
```

It needs Python 3 with `numpy`, `pandas` and `scipy`. It reads no raw
simulation output and uses no network. It draws no random numbers, so repeated
runs give byte-identical files. If an input is missing, or a fit it recomputes
disagrees with the notebook's exports by more than 10⁻¹², it stops with an error
and writes nothing.

## What it covers

| Section | Content | Source |
|---|---|---|
| 1 | Luminosity per configuration and emittance: mean, std, SE, seeds, full grid | `data/gp_luminosity_for_wx.csv`, cut multipliers from `inputs/acc_lumi_opt_C3_250_no_pairs.dat` |
| 2 | Conservative/nominal ratio, H_D with L_geom and its parameters, 20 → 2 nm conservative gain, seed scatter | section 1 rows |
| 3 | n_y^req and n_m^req with σ(ln) | `data/gp_calibration_export.csv` |
| 4 | Free and exponent-constrained n_y fits and the n_m fit, with χ² and ndf; the frozen conservative loci | recomputed from section 3; loci from `constants.py` |
| 5 | κ per emittance, weighted mean, slope of ln κ against ln D_y | section 3 |
| 6 | Recommendation table | frozen loci and D_y |
| 7 | D_y and its expression | beam of `inputs/acc_C3_250_nominal.dat` |
| 8 | IPC background table, eight rows with their ten seeds | `data/bib_stats_cache.csv` |
| 9 | 2 nm / 20 nm produced and reaching ratios, the BIB-run luminosity ratio, background-per-luminosity reductions | section 8 and the pairs-on BIB runs in `data/lumi_extracted.csv` |
| 10 | n_x convergence at 20 nm | `data/L_vs_nm_by_nx_20nm.csv` |

## What it deliberately does not cover

- **Other simulator.** WarpX quantities, which the WarpX repository reproduces.
- **Beyond the snapshot.** The 40–100 nm frozen extension, which is not a GUINEA-PIG++ configuration the paper reports. Also the conservative BIB reruns and the 143 runs added to the raw tree after the 2026-09-10 luminosity snapshot, none of which the paper uses.
- **Raw output and figures.** Regenerating the committed data from raw output, and regenerating figures. The README documents both.
- **Intermediates.** Per-seed values, intermediate caches and diagnostics.

## Verification

- **Clean clone.** From a clean clone, with the raw-output root pointed at a path that does not exist, the script completes and its output is byte-identical to the committed files.
- **Determinism.** A second run gives byte-identical output.
- **Failure.** Removing an input makes it stop with an error.
- **Internal consistency.** Every derived quantity was recomputed from the rows of the same JSON file and agrees: ratios, the gain, H_D, seed scatter, κ from n_y^req, the background ratios and reductions, and the recommendation table from the frozen loci.
- **Independent recomputation.** The luminosity rows, ratios, gain, IPC table and background quantities agree with an earlier recomputation from the raw `.ref` and `pairs.dat` files.
- **Published constants.** q_n = 0.402 ± 0.152, C_y = 26.2, χ² 1.02/5; C_y = 21.32 at q_n = 0.465, χ² 1.20/6; s = 3.329 ± 0.250, C_m = 1.4×10⁻⁷, χ² 7.49/5; κ = 0.284 ± 0.008 with slope −0.063 ± 0.152. All reproduce.

## What does not match the manuscript

These are reported, not adjusted:

| Quantity | Manuscript | This repository | Reason |
|---|---|---|---|
| Exclusion of q_n = 0 | 2.6σ | 2.65σ (0.40187 / 0.15153) | rounding |
| s_cons | 3.515 | 3.514 | s_cons = 3.300 + q_p; the repository adopts q_p = 0.214, while 3.515 needs 0.215. No run used s_cons, only the exponent 3.300. |
| n_m^rec at 0.5 nm | 4.5×10⁷ | 4.44×10⁷, which rounds to 4.4×10⁷ | the frozen locus gives the lower value |
| n_m^rec at 16 nm | 7.1×10⁴ | 7.03×10⁴, which rounds to 7.0×10⁴ | the frozen locus gives the lower value |

Other notes a reviewer may meet:

- **No reference for H_D.** Nothing committed records H_D values to check against. The script states its definition: L_geom = N² n_b f_rep / (4π σ_x σ_y), without an hourglass factor.
- **Two uncertainty sets.** The published fits use the per-emittance σ(ln n^req) of `data/gp_calibration_export.csv`. `data/gp_requirements_for_wx.csv`, the file exported to WarpX, instead carries an older uncertainty budget, `data/error_budget.csv`. For n_y it differs at 2, 4 and 20 nm, for example 0.551 against 0.534 at 2 nm.
- **Two D_y evaluations.** The committed D_y table was computed from the beam sizes echoed by each run. It agrees with the D_y expression to 2×10⁻⁶, but the two round differently at 8 nm: 34.215 in the table, 34.216 from the expression.
- **A notebook label.** The notebook's printed label for the exponent-constrained n_y fit says q = 0.214, but the fit fixes q_n = 0.465. The values it prints are for q_n = 0.465.
