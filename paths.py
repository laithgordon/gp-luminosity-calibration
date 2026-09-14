"""Where this repository finds its inputs.

Everything the analysis reads from the repository itself is located relative
to REPO_ROOT, so a fresh clone runs from any directory without edits.

The raw GUINEA-PIG++ output (.ref files and pairs.dat dumps, several GB) is not
in git. Set GP_RAW_ROOT to the directory that contains output/, output_nm/ and
output_no_pairs/, each with a C3_250/ subdirectory:

    export GP_RAW_ROOT=/path/to/GuineaPig_Feb_2025

If GP_RAW_ROOT is unset, the parent of the repository is used, which is the
layout the study was run in. This module is the only place that path is set.
"""
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
DATA_DIR = REPO_ROOT / "data"
PLOTS_DIR = REPO_ROOT / "plots"
RAW_ROOT = Path(os.environ.get("GP_RAW_ROOT", str(REPO_ROOT.parent))).expanduser().resolve()

COLLIDER = "C3_250"


def raw_dir(source_dir: str) -> Path:
    """Directory holding one campaign's raw output, e.g. raw_dir('output_nm')."""
    return RAW_ROOT / source_dir / COLLIDER
