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

The committed `paper_numbers.json` is reproduced to floating-point last-digit precision: a few values can differ in their final digit between platforms, and no quoted value depends on those digits.

## What it covers

| Section | Content | Source |
|---|---|---|
| 1 | Luminosity per configuration and emittance: mean, std, SE, seeds, full grid | `data/gp_luminosity_for_wx.csv`, cut multipliers from `inputs/acc_lumi_opt_C3_250_no_pairs.dat` |
| 2 | Conservative/nominal ratio, H_D with L_geom and its parameters, 20 → 2 nm conservative gain, seed scatter | section 1 rows |
| 3 | n_y^req and n_m^req with σ(ln) | `data/gp_calibration_export.csv` |
| 4 | Free and exponent-constrained n_y fits and the n_m fit, with χ² and ndf; the frozen conservative loci | recomputed from section 3; loci and q_p from `constants.py` |
| 5 | The inputs exported for κ: n_y^req, σ(ln) and D_y per emittance. κ itself, for both simulators, is computed in the WarpX repository from this export | `data/gp_requirements_for_wx.csv` |
| 6 | Recommendation table | frozen loci and D_y |
| 7 | D_y and its expression | beam of `inputs/acc_C3_250_nominal.dat` |
| 8 | IPC background table, eight rows with their ten seeds | `data/bib_stats_cache.csv` |
| 9 | 2 nm / 20 nm produced and reaching ratios, the BIB-run luminosity ratio, background-per-luminosity reductions | section 8 and the pairs-on BIB runs in `data/lumi_extracted.csv` |
| 10 | n_x convergence at 20 nm | `data/L_vs_nm_by_nx_20nm.csv` |

q_p = 0.2153 is the envelope-integration exponent, defined once in `constants.py`
and shared with the WarpX repository; q_n^pred = q_p + 1/4 = 0.4653 is derived
from it. The notebook, the export and this script all read both from there.

## What it deliberately does not cover

- **κ.** This repository exports κ's inputs and does not compute κ. The WarpX repository computes it for both simulators with a single definition of R(D_y) and a single weighting.
- **Other simulator.** WarpX quantities, which the WarpX repository reproduces.
- **Beyond the snapshot.** The 40–100 nm frozen extension, which is not a GUINEA-PIG++ configuration the paper reports. Also the conservative BIB reruns and the 143 runs added to the raw tree after the 2026-09-10 luminosity snapshot, none of which the paper uses.
- **Raw output and figures.** Regenerating the committed data from raw output, and regenerating figures. The README documents both.
- **Intermediates.** Per-seed values, intermediate caches and diagnostics.

## Verification

- **Clean clone.** From a clean clone, with the raw-output root pointed at a path that does not exist, the script completes and its output matches the committed files to floating-point last-digit precision.
- **Determinism.** A second run gives byte-identical output.
- **Failure.** Removing an input makes it stop with an error.
- **Internal consistency.** Every derived quantity was recomputed from the rows of the same JSON file and agrees: ratios, the gain, H_D, seed scatter, the background ratios and reductions, and the recommendation table from the frozen loci.
- **Independent recomputation.** The luminosity rows, ratios, gain, IPC table and background quantities agree with a recomputation from the raw `.ref` and `pairs.dat` files.
- **Published constants.** Free n_y fit: q_n = 0.402 ± 0.152, C_y = 26.2, χ² 1.02/5. Constrained fit at q_n = 0.4653: C_y = 21.3, χ² 1.20/6. n_m fit: s = 3.330 ± 0.250, C_m = 1.4×10⁻⁷, χ² 7.49/5. s_cons = 3.515. All reproduce.

## Notes

- **No reference for H_D.** Nothing committed records H_D values to check against. The script states its definition: L_geom = N² n_b f_rep / (4π σ_x σ_y), without an hourglass factor.
- **Two D_y evaluations.** The committed D_y table was computed from the beam sizes echoed by each run. It agrees with the D_y expression to 2×10⁻⁶, but the two round differently at 8 nm: 34.215 in the table, 34.216 from the expression.
- **Committed figures predate q_p = 0.2153.** They were made with q_p = 0.214 and q_n = 0.465, and were deliberately not regenerated. Rerunning the notebook now would change the drawn s_cons label in `nm_req_tuning` from 3.51 to 3.52, and the constrained-fit line in `ny_req_tuning`. It would also change the repository-only `kappa_vs_Dy`, whose manuscript version comes from the WarpX repository.
