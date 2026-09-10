# GUINEA-PIG++ inputs

| File | Purpose |
|------|---------|
| `acc_C3_250_nominal.dat` | The GP++ input file: **one** default parameter set — the nominal C³-250 beam and the nominal grid. |
| `submit_template_S3DF_SLURM.sh` | SLURM array template that alters the default per run: numeric parameters from the scan-file row, physics/output switches from its `SETTINGS`. |
| `example_scan.txt` | Scan-file format, with the nominal row and a small `n_m` sweep. |

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

The **calibrated** ("tuned") settings in the paper are not a separate file: they
are the nominal file with `n_y` and `n_m` overridden per ε_y from the
calibration laws, `n_z = 64`, exactly as the template does for any scan row.

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
