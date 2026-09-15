#!/usr/bin/env python3
"""
extract_lumi_data.py
--------------------
Scan .ref files in $GP_RAW_ROOT/{output,output_nm,output_no_pairs}/C3_250/ (paths.py),
parse run parameters from the filename and luminosity metrics from the content,
and write data/lumi_extracted.csv.

The CSV is REBUILT IN FULL on every run, so parse errors never compound and
deleted .ref files disappear from the output. What is cached is only the *file
reading*: a file whose (mtime, size) is unchanged since the last run is not
re-read from disk, its previously parsed row is reused. Combined with a thread
pool for the reads that remain, a warm run over ~16k files takes seconds rather
than ~20 minutes. Use --full to ignore the cache and re-read everything.

Usage
    python3 extract_lumi_data.py              # incremental; writes data/lumi_extracted.csv
    python3 extract_lumi_data.py --full       # re-read every .ref file
    python3 extract_lumi_data.py --workers 8  # fewer threads
    python3 extract_lumi_data.py --out /tmp/x.csv    # override destination(s)
    python3 extract_lumi_data.py --restrict-to data/lumi_extracted.csv --out /tmp/x.csv
                                              # re-extract exactly the runs listed in the
                                              # committed snapshot, ignoring later runs
    python3 extract_lumi_data.py --test       # 10 files/dir, writes nothing

UNITS (read this before using lumi_ee / lumi_fine)
    GP++ prints `lumi_ee` and `lumi_fine` in the .ref file as the luminosity of
    ONE bunch crossing in m^-2 (see src/resultsCPP.cc: f_rep*n_b is applied only
    to the lumi[j1][j2] matrix, never to lumi_ee). This script converts them to
    a rate in cm^-2 s^-1 using the f_rep and n_b echoed in the SAME .ref file:

        lumi_ee [cm^-2 s^-1] = lumi_ee_m2 [m^-2 / crossing] * 1e-4 * n_b * f_rep

    (= x 1.5960 for the C3-250 deck, n_b = 133, f_rep = 120 Hz). The raw values
    are kept in `lumi_ee_m2` / `lumi_fine_m2`, and `f_rep` / `n_b` are recorded
    per row. A row whose .ref lacks the f_rep/n_b echo gets an EMPTY lumi_ee
    and is listed in the report - it is never silently left unconverted.

Outputs
    data/lumi_extracted.csv             (what GP_CALIBRATION.ipynb reads)
    data/lumi_extract_report.txt        (unparsed filenames, rows missing lumi_ee)
    data/.lumi_extract_cache.json       (read cache; safe to delete any time)
"""

import argparse
import csv
import gzip
import json
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor

# ── Paths ──────────────────────────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import RAW_ROOT  # noqa: E402  (set GP_RAW_ROOT, see paths.py)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

SOURCE_DIRS = {key: os.path.join(str(RAW_ROOT), key, "C3_250")
               for key in ("output", "output_nm", "output_no_pairs")}

# Only this repository's copy is written by default. Pass --out (repeatable)
# to refresh a mirror elsewhere as well.
DEFAULT_OUT_CSVS = [os.path.join(SCRIPT_DIR, "data", "lumi_extracted.csv")]
CACHE_PATH  = os.path.join(SCRIPT_DIR, "data", ".lumi_extract_cache.json")
REPORT_PATH = os.path.join(SCRIPT_DIR, "data", "lumi_extract_report.txt")
CACHE_VERSION = 3          # bump to invalidate every cached row (3: lumi unit conversion)

# ── Filename regex patterns ────────────────────────────────────────────────────
# Common tail shared by both patterns (after the offsety value). Every optional
# group has been seen in production filenames; `grids` is the one whose absence
# silently dropped 280 grid-extent files before CACHE_VERSION 2.
_TAIL = (
    r"(?:_integration_method_(?P<integration_method>[0-9.eE+-]+))?"
    r"(?:_nt(?P<n_t>[0-9.eE+-]+))?"
    r"(?:_cutx(?P<cut_x>[0-9.eE+-]+)_cuty(?P<cut_y>[0-9.eE+-]+)_cutz(?P<cut_z>[0-9.eE+-]+))?"
    r"(?:_grids(?P<grids>[0-9.eE+-]+))?"
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
# Byte patterns: .ref files are ASCII, and skipping the UTF-8 decode of ~180 MB
# is measurably cheaper than decoding and then matching str patterns.
_NUM = rb"([0-9.eE+-]+)"
CONTENT_PATTERNS = [
    # raw per-crossing values exactly as printed (m^-2); converted in build_row
    ("lumi_ee_m2",   re.compile(rb"\blumi_ee\b\s*=\s*"   + _NUM), float),
    ("lumi_fine_m2", re.compile(rb"\blumi_fine\b\s*=\s*" + _NUM), float),
    # collider rate factors echoed by GP++ in the SWITCHES block of every .ref
    ("f_rep",     re.compile(rb"\bf_rep\s*=\s*"     + _NUM), float),
    ("n_b",       re.compile(rb"\bn_b\s*=\s*([0-9]+)"),        int),
    ("phot_e1",   re.compile(rb"\bphot-e1\b\s*=\s*"   + _NUM), float),
    ("phot_e2",   re.compile(rb"\bphot-e2\b\s*=\s*"   + _NUM), float),
    ("n_phot1",   re.compile(rb"\bn_phot1\b\s*=\s*"   + _NUM), float),
    ("n_phot2",   re.compile(rb"\bn_phot2\b\s*=\s*"   + _NUM), float),
    ("de1",       re.compile(rb"\bde1\b\s*=\s*"       + _NUM), float),
    ("de2",       re.compile(rb"\bde2\b\s*=\s*"       + _NUM), float),
    ("n_pairs",   re.compile(rb"\bn_pairs\b\s*=\s*"   + _NUM), float),
    ("e_pairs",   re.compile(rb"\be_pairs\b\s*=\s*"   + _NUM), float),
    ("out_1",     re.compile(rb"out\.1=([0-9]+)"),             int),
    ("out_2",     re.compile(rb"out\.2=([0-9]+)"),             int),
]

# ── CSV columns ───────────────────────────────────────────────────────────────
# Order is load-bearing for downstream analysis: everything through "out_2" is
# identical to the original schema. `grids` is appended LAST so any positional
# reader of the original columns is unaffected.
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
    "grids",
    # raw per-crossing values and rate factors.
    # lumi_ee / lumi_fine above are cm^-2 s^-1 = *_m2 * 1e-4 * n_b * f_rep.
    "lumi_ee_m2", "lumi_fine_m2", "f_rep", "n_b",
]

# m^-2 per crossing -> cm^-2 s^-1
def lumi_rate(L_m2, f_rep, n_b):
    return L_m2 * 1e-4 * n_b * f_rep

# ── Helpers ────────────────────────────────────────────────────────────────────

def parse_filename(fn):
    """Return dict of parameters, or None if no pattern matches."""
    for pat in (PATTERN_NM, PATTERN_NO_NM):
        m = pat.match(fn)
        if m:
            return m.groupdict()
    return None


def parse_content(path):
    """Read one .ref file; return dict of metrics, or None if unreadable."""
    try:
        with open(path, "rb") as fh:
            blob = fh.read()
    except OSError:
        return None
    result = {}
    for col, regex, cast in CONTENT_PATTERNS:
        m = regex.search(blob)
        if m:
            try:
                result[col] = cast(m.group(1))
            except (ValueError, OverflowError):
                pass
    return result


def build_row(source_key, fn, params, metrics):
    row = {col: "" for col in COLUMNS}
    row["source_dir"] = source_key
    row["filename"] = fn
    for k, v in params.items():
        if k in row and v is not None:
            row[k] = v
    for col in metrics:
        if col in row:
            row[col] = metrics[col]
    # Convert the per-crossing m^-2 values to cm^-2 s^-1 with the run's own
    # f_rep / n_b. Missing factors -> lumi_ee left blank (reported), never raw.
    f_rep, n_b = metrics.get("f_rep"), metrics.get("n_b")
    if f_rep is not None and n_b is not None:
        for raw, conv in (("lumi_ee_m2", "lumi_ee"), ("lumi_fine_m2", "lumi_fine")):
            if raw in metrics:
                row[conv] = lumi_rate(metrics[raw], f_rep, n_b)
    return row


def scan_dir(source_key, src_dir):
    """One pass over a source dir -> list of (key, fn, path, mtime_ns, size).
    os.scandir carries stat data in the dirent, so this is one syscall per entry
    instead of glob + a separate stat() per file."""
    out = []
    if not os.path.isdir(src_dir):
        print(f"  WARNING: missing source dir {src_dir}", file=sys.stderr)
        return out
    with os.scandir(src_dir) as it:
        for e in it:
            if not e.name.endswith(".ref"):
                continue
            try:
                st = e.stat()
            except OSError:
                continue
            out.append((f"{source_key}/{e.name}", e.name, e.path,
                        st.st_mtime_ns, st.st_size))
    return out


def load_cache(path, disabled=False):
    if disabled or not os.path.exists(path):
        return {}
    try:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            blob = json.load(fh)
        if blob.get("version") != CACHE_VERSION:
            print(f"  cache version {blob.get('version')} != {CACHE_VERSION}; ignoring it")
            return {}
        return blob.get("entries", {})
    except (OSError, ValueError, EOFError) as exc:
        print(f"  cache unreadable ({exc}); ignoring it")
        return {}


def save_cache(path, entries):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump({"version": CACHE_VERSION, "entries": entries}, fh)
    os.replace(tmp, path)


def sanity_check(records):
    """Print mean lumi_ee by emitt_y."""
    from collections import defaultdict
    buckets = defaultdict(list)
    for r in records:
        ey, lee = r.get("emitt_y", ""), r.get("lumi_ee", "")
        if ey != "" and lee != "":
            try:
                buckets[float(ey)].append(float(lee))
            except (ValueError, TypeError):
                pass
    print("\n--- Sanity check: mean lumi_ee by emitt_y  [1e34 cm^-2 s^-1] ---")
    for ey in sorted(buckets):
        vals = buckets[ey]
        print(f"  emitt_y={ey:.4f} mm.mrad  n={len(vals):5d}  "
              f"mean_lumi_ee={sum(vals)/len(vals)/1e34:.4f}e34")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--full", action="store_true",
                    help="ignore the read cache and re-read every .ref file")
    ap.add_argument("--workers", type=int,
                    default=min(32, (os.cpu_count() or 4) * 4),
                    help="thread-pool size for file reads (default: %(default)s)")
    ap.add_argument("--out", action="append", metavar="PATH",
                    help="CSV destination; repeatable. Default: data/lumi_extracted.csv")
    ap.add_argument("--restrict-to", metavar="CSV",
                    help="extract only the runs (source_dir, filename) listed in CSV, e.g. the "
                         "committed snapshot data/lumi_extracted.csv; runs added to the raw "
                         "tree later are ignored")
    ap.add_argument("--test", action="store_true",
                    help="parse 10 files per dir, print a sample, write nothing")
    args = ap.parse_args()

    out_csvs = args.out if args.out else DEFAULT_OUT_CSVS
    t0 = time.time()

    # ── inventory ──
    print("=== INVENTORY ===")
    listing = []
    for key, src in SOURCE_DIRS.items():
        found = scan_dir(key, src)
        print(f"  {key}: {len(found)} .ref files")
        listing.extend(found)
    print(f"  total: {len(listing)} .ref files  ({time.time()-t0:.1f}s)")

    if args.restrict_to:
        with open(args.restrict_to, newline="") as fh:
            wanted = {f"{r['source_dir']}/{r['filename']}" for r in csv.DictReader(fh)}
        on_disk = {x[0] for x in listing}
        listing = [x for x in listing if x[0] in wanted]
        absent = sorted(wanted - on_disk)
        print(f"  --restrict-to {args.restrict_to}: {len(wanted)} runs listed, "
              f"{len(listing)} found on disk, {len(absent)} missing, "
              f"{len(on_disk - wanted)} later runs ignored")
        for key in absent[:20]:
            print(f"    missing: {key}")

    if args.test:
        print("\n=== TEST RUN (10 files per dir) ===")
        for key in SOURCE_DIRS:
            sample = [x for x in listing if x[0].startswith(key + "/")][:10]
            ok = fail = 0
            for _, fn, path, _, _ in sample:
                p = parse_filename(fn)
                if p is None:
                    fail += 1
                    continue
                ok += 1
                if ok == 1:
                    m = parse_content(path) or {}
                    print(f"  {key} sample: emitt_y={p['emitt_y']} "
                          f"lumi_ee_m2={m.get('lumi_ee_m2')} f_rep={m.get('f_rep')} "
                          f"n_b={m.get('n_b')}")
            print(f"  {key}: {ok} parsed, {fail} failed")
        print("\n--test: nothing written.")
        return

    # ── decide what needs re-reading ──
    cache = load_cache(CACHE_PATH, disabled=args.full)
    reuse, todo = {}, []
    for key, fn, path, mtime_ns, size in listing:
        hit = cache.get(key)
        if hit and hit.get("mtime_ns") == mtime_ns and hit.get("size") == size:
            reuse[key] = hit
        else:
            todo.append((key, fn, path, mtime_ns, size))
    print(f"\n=== READ PLAN ===")
    print(f"  reused from cache (unchanged): {len(reuse)}")
    print(f"  to read (new / changed / uncached): {len(todo)}")

    # ── read what changed, in parallel (I/O bound -> threads) ──
    unparsed, no_lumi, no_rate, read_err = [], [], [], []
    fresh = {}
    if todo:
        def work(item):
            key, fn, path, mtime_ns, size = item
            params = parse_filename(fn)
            if params is None:
                return key, fn, None, None, mtime_ns, size
            return key, fn, params, parse_content(path), mtime_ns, size

        done = 0
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for key, fn, params, metrics, mtime_ns, size in pool.map(work, todo):
                done += 1
                if done % 2000 == 0 or done == len(todo):
                    print(f"  ... {done}/{len(todo)} read ({time.time()-t0:.1f}s)")
                if params is None:
                    unparsed.append(fn)
                    continue
                if metrics is None:
                    read_err.append(fn)
                    continue
                if "lumi_ee_m2" not in metrics:
                    no_lumi.append(fn)
                elif "f_rep" not in metrics or "n_b" not in metrics:
                    no_rate.append(fn)
                fresh[key] = {"mtime_ns": mtime_ns, "size": size,
                              "row": build_row(key.split("/", 1)[0], fn,
                                               params, metrics)}

    # ── rebuild the full record set from cache + fresh ──
    entries = dict(reuse)
    entries.update(fresh)
    src_order = {k: i for i, k in enumerate(SOURCE_DIRS)}
    records = [entries[k]["row"] for k in
               sorted(entries,
                      key=lambda k: (src_order.get(k.split("/", 1)[0], 99), k))]

    print(f"\n=== RESULT ===")
    print(f"  rows written:          {len(records)}")
    print(f"  unparsed filenames:    {len(unparsed)}")
    print(f"  parsed but no lumi_ee: {len(no_lumi)}")
    print(f"  lumi but no f_rep/n_b (lumi_ee left blank): {len(no_rate)}")
    print(f"  unreadable:            {len(read_err)}")

    # ── write CSVs atomically ──
    for out_csv in out_csvs:
        os.makedirs(os.path.dirname(out_csv), exist_ok=True)
        tmp = out_csv + ".tmp"
        with open(tmp, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLUMNS)
            w.writeheader()
            w.writerows(records)
        os.replace(tmp, out_csv)
        print(f"  saved: {out_csv}")

    # A restricted run must not shrink the read cache for later full runs.
    save_cache(CACHE_PATH, {**load_cache(CACHE_PATH), **entries} if args.restrict_to else entries)

    # ── failure report ──
    os.makedirs(os.path.dirname(REPORT_PATH), exist_ok=True)
    with open(REPORT_PATH, "w") as fh:
        fh.write("extract_lumi_data.py report\n")
        fh.write(f"files seen: {len(listing)}   rows written: {len(records)}\n")
        fh.write(f"re-read this run: {len(todo)}   reused from cache: {len(reuse)}\n\n")
        for title, items in (("UNPARSED FILENAMES", unparsed),
                             ("PARSED BUT NO lumi_ee", no_lumi),
                             ("lumi_ee PRESENT BUT NO f_rep/n_b ECHO (lumi_ee left blank, lumi_ee_m2 kept)", no_rate),
                             ("UNREADABLE", read_err)):
            fh.write(f"== {title} ({len(items)}) ==\n")
            for fn in sorted(items):
                fh.write(f"  {fn}\n")
            fh.write("\n")
    print(f"  report: {REPORT_PATH}")

    sanity_check(records)
    print(f"\nelapsed: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
