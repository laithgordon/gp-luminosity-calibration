# GUINEA-PIG++ inputs

| File | Purpose |
|------|---------|
| `acc_C3_250_nominal.dat` | One default parameter set: the nominal C³-250 beam and the nominal grid. |
| `acc_lumi_opt_C3_250_no_pairs.dat` | The input file the luminosity campaigns ran with (pairs, photons and dumps off). |
| `acc_lumi_opt_C3_250_do_dump_opt.dat` | The input file the BIB pair-dump campaigns ran with (pairs, photons, hadrons, jets, coherent and dumps on). |
| `submit_template_S3DF_SLURM.sh` | SLURM array template: numeric parameters from each scan-file row, physics/output switches from its `SETTINGS`. |
| `example_scan.txt` | Scan-file format, with the nominal row and a small `n_m` sweep. |
| `decks/` | The production scan files of every published configuration (table below). |
| `make_conservative_decks.py` | Writes the conservative-locus rows from the frozen constants in `constants.py`; `--check` compares them with `decks/round4_cons_*.txt`. |
| `run_configurations.csv` | Every distinct configuration in `data/lumi_extracted.csv` with its seeds: the ladder runs behind the fitted constants and figures. |
| `list_run_configurations.py` | Rebuilds `run_configurations.csv` from `data/lumi_extracted.csv`. |

## Nominal (default) parameters

**Beam (`$ACCELERATOR:: C3_250`)** — 125 GeV e⁻/e⁺, N = 0.624×10¹⁰ per bunch,
ε_x = 0.9 µm·rad, ε_y = 8 nm·rad, β_x/β_y = 12/0.12 mm, σ_z = 100 µm,
energy spread 0.3 %, beam-1 polarisation −0.8, f_rep = 120 Hz × n_b = 133,
zero offset. The ε_y scan in the paper varies `emitt_y` from 0.0005 to 0.1
(0.5–100 nm·rad) with everything else fixed.

**Grid (`$PARAMETERS:: Jim_pars_Aug2023`)** — the pre-calibration operating
point:

| parameter | value | meaning |
|---|---|---|
| `n_x`, `n_y`, `n_z` | 512, 512, 25 | transverse and longitudinal cells |
| `n_m` | 100 000 | macroparticles per beam |
| `n_t` | 6 | push sub-steps per slice crossing |
| `integration_method` | 2 | FFT field solver |
| `cut_x`, `cut_y`, `cut_z` | 20 σ_x, 20 σ_y, 3.5 σ_z | half-widths of the grid |
| `do_pairs` … `do_coherent` | 0 | secondary production off (luminosity runs) |

Three switches are off in the default and are turned on from the submit
script rather than by editing the file:

| `SETTINGS` in the template | alters | used for |
|---|---|---|
| `PAIRS=1` | `do_pairs`, `track_pairs`, `store_pairs`, `do_photons`, `do_hadrons`, `do_jets`, `do_coherent` → 1; writes `pairs.dat` | the BIB tables |
| `DUMP=1` (+ `DUMP_STRIDE`) | `do_dump=1`, `dump_particle=<stride>`, `dump_step=1`; one `b1.N`/`b2.N` particle file per slice stage | the pinched-beam-size analysis |
| `SIZE_LOG=1` | `do_size_log=1`; per-step rms size per slice to `beamsize1/2.dat` | cross-check of the size analysis |

The tuned and conservative settings are not separate input files: they are an
input file with `n_y`, `n_z` and `n_m` overridden per ε_y by a scan-file row.
Which rule set each published configuration is tabulated in the top-level
`README.md`. The production scan files are in `decks/`:

| Deck | Configuration | Used for |
|---|---|---|
| `time_step_study_nt6.txt` | nominal 512×512×25, `n_m` = 10⁵, ε_y = 0.5–16 nm | nominal luminosity |
| `pair_dump_lumi_bib_grid.txt` | nominal grid at ε_y = 20 nm, pairs on | nominal luminosity at 20 nm |
| `round4_cons_{2,4,8,12,16,20}nm.txt` | conservative locus | conservative luminosity |
| `phase9_batch3.txt` | ladder rungs, including 512×512×64, `n_m` = 1.41×10⁷ at 1 nm | conservative luminosity at 1 nm |
| `pair_dump_tuned_highemit.txt` | 512×128×64, `n_m` = 5×10⁴, ε_y = 40–100 nm | frozen extension |
| `pair_dump_2nm.txt`, `pair_dump_4nm.txt`, `pair_dump_calibrated.txt` | pre-campaign tuned grid at 2, 4, 8, 20 nm, pairs on | tuned rows of `tab:bib_yields` |

## Running

```bash
# from the directory holding build/bin/guinea and this inputs/ folder
sbatch inputs/submit_template_S3DF_SLURM.sh
# then
python3 extract_lumi_data.py        # parses the .ref files into data/lumi_extracted.csv
```

Fill in `--account`, `--partition`, `--mail-user`, and the `SETTINGS` block at
the top of the template. The output filename encodes every varied quantity and
is what `extract_lumi_data.py` parses, so keep the `stem` pattern if you change
the script.

Three things the template does that matter:

- **Per-task scratch directory.** GP++ writes fixed-name files (`pairs.dat`,
  `beamsize*.dat`, `b1.N`) into the working directory; concurrent tasks sharing
  one directory silently overwrite each other's output.
- **Idempotent.** A row whose `.ref` already exists is skipped, so a failed or
  extended campaign is resubmitted with the same command.
- **Non-preemptable QoS.** GP++ writes its results only on exit; a preempted
  multi-day run yields nothing.
