#!/bin/bash
#
# Template SLURM array job for the GUINEA-PIG++ parameter scans in this study.
# One array task per seed; each task walks every row of the scan file, so a
# 15-row array with a 4-row scan file is 60 runs.
#
# Edit the four SBATCH lines marked <...> and the three SETTINGS below.
#
#SBATCH --account=<account>            # e.g. atlas:usatlas
#SBATCH --partition=<partition>        # e.g. roma (CPU)
#SBATCH --qos=normal                   # see note on preemption below
#SBATCH --job-name=GP_scan
#SBATCH --output=GP_scan-%A_%a.txt
#SBATCH --error=GP_scan-%A_%a.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16g              # 48g was needed for n_m >= 7e6
#SBATCH --time=2-00:00:00
#SBATCH --array=0-14                   # one task per seed in SEEDS below
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL
#
# PREEMPTION. On S3DF the account decides the QoS: atlas:default only offers
# `preemptable`, under which two multi-day campaigns in this study were killed
# and lost (GP++ writes its outputs at the very end, so a preempted run yields
# nothing). Use an account that offers qos=normal for anything longer than a
# few hours.

set -u
SUBMIT_DIR=${SLURM_SUBMIT_DIR:-$PWD}
cd "$SUBMIT_DIR"

# ── SETTINGS ────────────────────────────────────────────────────────────────
GUINEA="$SUBMIT_DIR/build/bin/guinea"                  # GP++ executable
ACC_FILE="$SUBMIT_DIR/inputs/acc_C3_250_nominal.dat"   # the ONE default input
SCAN="$SUBMIT_DIR/inputs/example_scan.txt"             # one run per row
OUTDIR="$SUBMIT_DIR/output_nm/C3_250"                  # where .ref files land
#
# Physics / output switches. The default input has all three OFF; set to 1 to
# alter them for every run of this submission (numeric parameters are changed
# per row through the scan file instead).
PAIRS=0      # 1: pair production + photons/hadrons/jets/coherent, writes pairs.dat
             #    (the BIB studies). Roughly 3x the run time at n_m ~ 1e5.
DUMP=0       # 1: dump every macroparticle once per slice stage (b1.N / b2.N),
             #    used for the pinched-beam-size analysis. DUMP_STRIDE keeps only
             #    every k-th particle; without it n_m = 2e7 writes ~400 GB.
DUMP_STRIDE=1
SIZE_LOG=0   # 1: per-step rms beam size per slice (beamsize1.dat / beamsize2.dat)
# ────────────────────────────────────────────────────────────────────────────

# The 15 canonical seeds. Every point in the paper is a mean over exactly
# these; do not substitute.
SEEDS=(100805 105469 130471 135135 143455 145992 153768 159302 167029 167815 527875 733579 360177 286038 479728)
SEED=${SEEDS[$SLURM_ARRAY_TASK_ID]}
collider="C3_250"
mkdir -p "$OUTDIR"

rows=()
while IFS= read -r line; do
    [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
    rows+=("$line")
done < "$SCAN"

# GP++ writes several fixed-name files (pairs.dat, beamsize*.dat, b1.N ...) into
# the working directory, so concurrent tasks MUST NOT share one. Each task gets
# a scratch directory that is removed on exit.
WORK="$SUBMIT_DIR/.gp_work/${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$WORK"; trap 'rm -rf "$WORK"' EXIT

for row in "${rows[@]}"; do
    IFS=' ' read -r -a p <<< "$row"
    particles=${p[0]}; emitt_x=${p[1]}; emitt_y=${p[2]}
    beta_x=${p[3]};    beta_y=${p[4]};  sigma_z=${p[5]}; offset_y=${p[6]}
    n_x=${p[7]}; n_y=${p[8]}; n_z=${p[9]}; n_m=${p[10]}; im=${p[11]}; n_t=${p[12]}

    # Output name carries every varied quantity; extract_lumi_data.py parses it.
    stem="${n_x}_${n_y}_${n_z}_${n_m}_test_${collider}_particles_${particles}_emittx_${emitt_x}_emitty_${emitt_y}_betax_${beta_x}_betay_${beta_y}_sigmaz_${sigma_z}_offsety_${offset_y}_integration_method_${im}_nt${n_t}_seed_${SEED}"
    out_ref="$OUTDIR/${stem}.ref"
    if [ -s "$out_ref" ]; then echo "skip (exists): $stem"; continue; fi   # idempotent
    echo "run: $stem"

    # Per-run copy of the input file with the row's values substituted.
    cfg="$WORK/acc_${stem}.dat"
    cp "$ACC_FILE" "$cfg"
    sed -i "s/particles=[^;]*/particles=$particles/;s/emitt_x=[^;]*/emitt_x=$emitt_x/" "$cfg"
    sed -i "s/emitt_y=[^;]*/emitt_y=$emitt_y/;s/beta_x=[^;]*/beta_x=$beta_x/"       "$cfg"
    sed -i "s/beta_y=[^;]*/beta_y=$beta_y/;s/sigma_z=[^;]*/sigma_z=$sigma_z/"       "$cfg"
    sed -i "s/offset_y=[^;]*/offset_y=$offset_y/;s/rndm_seed=[0-9]\+/rndm_seed=$SEED/g" "$cfg"
    sed -i "s/\(n_x=\)[^;]*;/\1${n_x};/g;s/\(n_y=\)[^;]*;/\1${n_y};/g"             "$cfg"
    sed -i "s/\(n_z=\)[^;]*;/\1${n_z};/g;s/\(n_m=\)[^;]*;/\1${n_m};/g"             "$cfg"
    sed -i "s/\(n_t=[[:space:]]*\)[^;]*;/\1${n_t};/g;s/\(integration_method=\)[^;]*;/\1${im};/g" "$cfg"
    # ── switches (alter the default input's off-by-default settings) ──
    if [ "$PAIRS" = 1 ]; then
        sed -i "s/do_pairs=[^;]*/do_pairs=1/;s/track_pairs=[^;]*/track_pairs=1/;s/store_pairs=[^;]*/store_pairs=1/" "$cfg"
        sed -i "s/do_photons=[^;]*/do_photons=1/;s/do_hadrons=[^;]*/do_hadrons=1/;s/do_jets=[^;]*/do_jets=1/;s/do_coherent=[^;]*/do_coherent=1/" "$cfg"
    fi
    if [ "$DUMP" = 1 ]; then
        sed -i "s/do_dump=[^;]*/do_dump=1/;s/dump_particle=[^;]*/dump_particle=${DUMP_STRIDE}/;s/dump_step=[^;]*/dump_step=1/" "$cfg"
    fi
    [ "$SIZE_LOG" = 1 ] && sed -i "s/do_size_log=[^;]*/do_size_log=1/" "$cfg"

    START=$(date +%s)
    ( cd "$WORK" && "$GUINEA" --acc_file "$cfg" "$collider" "Jim_pars_Aug2023" "$out_ref" )
    echo "TIMING_RESULT_GUINEA: $(( $(date +%s) - START )) seconds"

    # Collect the optional outputs beside the .ref under matching names.
    [ -s "$WORK/pairs.dat" ] && mv "$WORK/pairs.dat" \
        "$OUTDIR/${stem/test_${collider}_particles_/test${collider}_pairs_particles_}.dat"
    if [ "$DUMP" = 1 ]; then
        mkdir -p "$OUTDIR/dumps/$stem" && mv "$WORK"/b1.* "$WORK"/b2.* "$OUTDIR/dumps/$stem/" 2>/dev/null
    fi
    for b in 1 2; do
        [ -s "$WORK/beamsize$b.dat" ] && mv "$WORK/beamsize$b.dat" "$OUTDIR/beamsize${b}_${stem}.dat"
    done
    rm -f "$cfg" "$WORK"/b1.* "$WORK"/b2.* "$WORK"/bp1.* "$WORK"/bp2.*
    echo "  -> $(grep -o 'lumi_ee=[0-9.eE+-]*' "$out_ref" 2>/dev/null | head -1)"
done
echo "SEED $SEED COMPLETE"
