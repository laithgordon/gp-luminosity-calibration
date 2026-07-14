"""
extract_lumi_data.py
--------------------
Scan .ref files in output/C3_250/, output_nm/C3_250/, output_no_pairs/C3_250/,
parse filename parameters and content metrics, save to analysis/data/lumi_extracted.csv.
"""

import os
import re
import glob
import csv
import math

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE = "/fs/ddn/sdf/group/atlas/d/laithg/GuineaPig_Feb_2025"
SOURCE_DIRS = {
    "output":       os.path.join(BASE, "output",         "C3_250"),
    "output_nm":    os.path.join(BASE, "output_nm",      "C3_250"),
    "output_no_pairs": os.path.join(BASE, "output_no_pairs", "C3_250"),
}
OUT_CSV = os.path.join(BASE, "analysis", "data", "lumi_extracted.csv")
BATCH_SIZE = 500

# ── Filename regex patterns ────────────────────────────────────────────────────
# Common tail shared by both patterns (after offsety value):
#   optional: _integration_method_{val}
#   optional: _cutx{val}_cuty{val}_cutz{val}
#   required: _seed_{seed}.ref
_TAIL = (
    r"(?:_integration_method_(?P<integration_method>[0-9.eE+-]+))?"
    r"(?:_nt(?P<n_t>[0-9.eE+-]+))?"
    r"(?:_cutx(?P<cut_x>[0-9.eE+-]+)_cuty(?P<cut_y>[0-9.eE+-]+)_cutz(?P<cut_z>[0-9.eE+-]+))?"
    r"_seed_(?P<seed>\d+)\.ref$"
)

_COMMON_PARAMS = (
    r"particles_(?P<particles>[0-9.eE+-]+)_"
    r"emittx_(?P<emitt_x>[0-9.eE+-]+)_emitty_(?P<emitt_y>[0-9.eE+-]+)_"
    r"betax_(?P<beta_x>[0-9.eE+-]+)_betay_(?P<beta_y>[0-9.eE+-]+)_"
    r"sigmaz_(?P<sigma_z>[0-9.eE+-]+)_offsety_(?P<offset_y>[0-9.eE+-]+)"
)

# Pattern A: no n_m field   {nx}_{ny}_{nz}_test_C3_250_...
PATTERN_NO_NM = re.compile(
    r"^(?P<n_x>\d+)_(?P<n_y>\d+)_(?P<n_z>\d+)_test_C3_250_"
    + _COMMON_PARAMS + _TAIL
)

# Pattern B: with n_m field  {nx}_{ny}_{nz}_{nm}_test_C3_250_...
PATTERN_NM = re.compile(
    r"^(?P<n_x>\d+)_(?P<n_y>\d+)_(?P<n_z>\d+)_(?P<n_m>\d+)_test_C3_250_"
    + _COMMON_PARAMS + _TAIL
)

# ── Content regex patterns ─────────────────────────────────────────────────────
RE_LUMI_EE   = re.compile(r"\blumi_ee\b\s*=\s*([0-9.eE+-]+)")
RE_LUMI_FINE = re.compile(r"\blumi_fine\b\s*=\s*([0-9.eE+-]+)")
RE_PHOT_E1   = re.compile(r"\bphot-e1\b\s*=\s*([0-9.eE+-]+)")
RE_PHOT_E2   = re.compile(r"\bphot-e2\b\s*=\s*([0-9.eE+-]+)")
RE_N_PHOT1   = re.compile(r"\bn_phot1\b\s*=\s*([0-9.eE+-]+)")
RE_N_PHOT2   = re.compile(r"\bn_phot2\b\s*=\s*([0-9.eE+-]+)")
RE_DE1       = re.compile(r"\bde1\b\s*=\s*([0-9.eE+-]+)")
RE_DE2       = re.compile(r"\bde2\b\s*=\s*([0-9.eE+-]+)")
RE_N_PAIRS   = re.compile(r"\bn_pairs\b\s*=\s*([0-9.eE+-]+)")
RE_E_PAIRS   = re.compile(r"\be_pairs\b\s*=\s*([0-9.eE+-]+)")
RE_OUT1      = re.compile(r"out\.1=([0-9]+)")
RE_OUT2      = re.compile(r"out\.2=([0-9]+)")

CONTENT_PATTERNS = [
    ("lumi_ee",   RE_LUMI_EE,   float),
    ("lumi_fine", RE_LUMI_FINE, float),
    ("phot_e1",   RE_PHOT_E1,   float),
    ("phot_e2",   RE_PHOT_E2,   float),
    ("n_phot1",   RE_N_PHOT1,   float),
    ("n_phot2",   RE_N_PHOT2,   float),
    ("de1",       RE_DE1,       float),
    ("de2",       RE_DE2,       float),
    ("n_pairs",   RE_N_PAIRS,   float),
    ("e_pairs",   RE_E_PAIRS,   float),
    ("out_1",     RE_OUT1,      int),
    ("out_2",     RE_OUT2,      int),
]

# ── CSV columns ───────────────────────────────────────────────────────────────
COLUMNS = [
    "source_dir", "filename",
    "n_x", "n_y", "n_z", "n_m",
    "particles",
    "emitt_x", "emitt_y",
    "beta_x", "beta_y",
    "sigma_z", "offset_y",
    "seed",
    "integration_method", "n_t",
    "cut_x", "cut_y", "cut_z",
    "lumi_ee", "lumi_fine",
    "phot_e1", "phot_e2",
    "n_phot1", "n_phot2",
    "de1", "de2",
    "n_pairs", "e_pairs",
    "out_1", "out_2",
]

# ── Helpers ────────────────────────────────────────────────────────────────────

def parse_filename(fn):
    """Return dict of parameters or None if no pattern matches."""
    for pat in (PATTERN_NM, PATTERN_NO_NM):
        m = pat.match(fn)
        if m:
            return m.groupdict()
    return None


def parse_content(path):
    """Read .ref file and return dict of extracted metric values."""
    try:
        with open(path, "r", errors="replace") as fh:
            content = fh.read()
    except OSError:
        return {}
    result = {}
    for col, regex, cast in CONTENT_PATTERNS:
        m = regex.search(content)
        if m:
            try:
                result[col] = cast(m.group(1))
            except (ValueError, OverflowError):
                pass
    return result


def process_files(paths, source_key, test_mode=False):
    """Process a list of .ref paths; return (records, n_parsed, n_failed)."""
    records = []
    n_parsed = 0
    n_failed = 0
    limit = 10 if test_mode else len(paths)
    for path in paths[:limit]:
        fn = os.path.basename(path)
        params = parse_filename(fn)
        if params is None:
            n_failed += 1
            continue
        metrics = parse_content(path)
        row = {col: "" for col in COLUMNS}
        row["source_dir"] = source_key
        row["filename"] = fn
        # filename params
        for k, v in params.items():
            if k in row and v is not None:
                row[k] = v
        # content metrics
        for col, _, _ in CONTENT_PATTERNS:
            if col in metrics:
                row[col] = metrics[col]
        n_parsed += 1
        records.append(row)
    return records, n_parsed, n_failed


def sanity_check(records):
    """Print mean lumi_ee by emitt_y."""
    from collections import defaultdict
    buckets = defaultdict(list)
    for r in records:
        ey = r.get("emitt_y", "")
        lee = r.get("lumi_ee", "")
        if ey != "" and lee != "":
            try:
                buckets[float(ey)].append(float(lee))
            except (ValueError, TypeError):
                pass
    print("\n--- Sanity check: mean lumi_ee by emitt_y ---")
    for ey in sorted(buckets):
        vals = buckets[ey]
        mean = sum(vals) / len(vals)
        print(f"  emitt_y={ey:.4f} mm.mrad  n={len(vals):5d}  mean_lumi_ee={mean:.4e}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    # ── Test run on 10 files first ──
    print("=== TEST RUN (10 files each dir) ===")
    test_ok = True
    for src_key, src_dir in SOURCE_DIRS.items():
        sample = glob.glob(os.path.join(src_dir, "*.ref"))[:10]
        recs, nparsed, nfailed = process_files(sample, src_key, test_mode=True)
        print(f"  {src_key}: {nparsed} parsed, {nfailed} failed")
        if recs:
            r0 = recs[0]
            print(f"    sample: emitt_y={r0['emitt_y']}  lumi_ee={r0['lumi_ee']}")
        if nfailed > 0 and nparsed == 0:
            print(f"  WARNING: all test files failed to parse in {src_key}")
            test_ok = False

    if not test_ok:
        print("Aborting — fix filename patterns before full run.")
        return

    print("\n=== FULL RUN ===")
    all_records = []
    total_parsed = 0
    total_failed = 0

    for src_key, src_dir in SOURCE_DIRS.items():
        all_paths = sorted(glob.glob(os.path.join(src_dir, "*.ref")))
        n_total = len(all_paths)
        print(f"\n{src_key}: {n_total} .ref files")
        src_parsed = 0
        src_failed = 0

        # process in batches
        for batch_start in range(0, n_total, BATCH_SIZE):
            batch = all_paths[batch_start : batch_start + BATCH_SIZE]
            recs, np_, nf = process_files(batch, src_key)
            all_records.extend(recs)
            src_parsed += np_
            src_failed += nf
            if (batch_start // BATCH_SIZE + 1) % 5 == 0 or batch_start + BATCH_SIZE >= n_total:
                print(f"  ... {batch_start + len(batch)}/{n_total} processed")

        print(f"  parsed={src_parsed}  failed={src_failed}")
        total_parsed += src_parsed
        total_failed += src_failed

    print(f"\nTotal: {total_parsed} parsed, {total_failed} failed")

    # ── Write CSV ──
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(all_records)
    print(f"Saved: {OUT_CSV}  ({len(all_records)} rows)")

    # ── Sanity check ──
    sanity_check(all_records)


if __name__ == "__main__":
    main()
