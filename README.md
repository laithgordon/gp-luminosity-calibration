# GUINEA-PIG++ Luminosity Calibration Analysis

Analysis code and paper-ready figures for the GUINEA-PIG++ beam-beam luminosity
calibration study (C3_250 collider configuration). The study calibrates the
macroparticle count and grid resolution required for converged luminosity as a
function of vertical beam emittance.

## Contents

| Path | Description |
|------|-------------|
| `extract_lumi_data.py` | Scans GUINEA-PIG++ `.ref` output files, parses run parameters from filenames and luminosity metrics from file contents, and writes `data/lumi_extracted.csv`. |
| `GP_CALIBRATION.ipynb` | Main analysis notebook. Produces four calibration figures and saves each to `plots/` as PDF + PNG. |
| `data/lumi_extracted.csv` | Extracted luminosity metrics (one row per simulation run). Processed summary — **not** raw simulation output. |
| `plots/` | Final rendered figures (PDF + PNG). |

## Figures produced

1. **Convergence scan** — sweep each grid parameter (`n_m`, `n_z`, `n_x`, `n_y`) and locate where luminosity `L` converges; fits `L = L∞ − C/nᵖ`.
2. **Refined scan** — re-sweep each parameter around the converged region.
3. **Required macroparticle count** `n_m^req` vs vertical emittance `ε_y`.
4. **Required vertical resolution** `n_y^req` vs vertical emittance `ε_y`.

## Regenerating the data

The raw simulation data (the `.ref` files, `beam*.dat`, and `output*/` directories)
is **not** included in this repository — it is large and produced by GUINEA-PIG++
runs. To rebuild `data/lumi_extracted.csv`, point the `SOURCE_DIRS` paths at the
top of `extract_lumi_data.py` at your GUINEA-PIG++ output directories and run:

```bash
python3 extract_lumi_data.py
```

Then run `GP_CALIBRATION.ipynb` to regenerate the figures.

## Requirements

- Python 3
- `numpy`, `pandas`, `matplotlib`, `scipy` (for the notebook)
- Jupyter (to run the notebook)
