#!/usr/bin/env python3
"""Convert TCGA WGS filtered SV tables into SVsig input CSVs.

Reads   <in_dir>/<COHORT>.new.merged.filtered.SV.tsv        (columns chr1,pos1,str1,chr2,pos2,str2,individual)
Writes  <out_dir>/<COHORT>.new.merged.filtered.SV_svsig.csv (the 15-column SVsig format)

Chromosomes: X->23, Y->24; SVs touching the mitochondrion (M) are dropped.
Strands: 0 -> '+', anything else -> '-'.

Usage:
  python3 data/datasets_to_proper_format_SVsig.py [-i IN_DIR] [-o OUT_DIR] [COHORT ...]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_IN_DIR = "/Users/akowalew/Documents/Getz Lab/1. TCGA-WGS/Data/NEW_filtered/"
DEFAULT_OUT_DIR = Path(__file__).resolve().parent / "input_data"
SUFFIX = ".new.merged.filtered.SV.tsv"


def chr_to_num(x):
    if x == "X":
        return 23
    if x == "Y":
        return 24
    return int(float(x))


def convert_to_svsig(in_path: Path, cohort_name: str, out_dir: Path) -> Path:
    df = pd.read_table(in_path)

    # Drop SVs with a mitochondrial breakpoint (compare as strings: the column may be mixed-type)
    mito = (df["chr1"].astype(str) == "M") | (df["chr2"].astype(str) == "M")
    df = df.loc[~mito].reset_index(drop=True)
    n = df.shape[0]

    out = pd.DataFrame({
        "seqnames": df["chr1"].astype(str).map(chr_to_num),
        "start": df["pos1"],
        "strand": np.where(df["str1"] == 0, "+", "-"),
        "altchr": df["chr2"].astype(str).map(chr_to_num),
        "altpos": df["pos2"],
        "altstrand": np.where(df["str2"] == 0, "+", "-"),
        "dcc_project_code": cohort_name,
        "sv_id": np.arange(n),
        "sid": df["individual"],
        "donor_unique_id": df["individual"],
        "weights": np.ones(n),
        "topo": np.zeros(n),
        "topo_n": np.zeros(n),
        "mech": np.zeros(n),
        "homseq": np.zeros(n),
    })

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / (in_path.name[:-len(".tsv")] + "_svsig.csv")
    out.to_csv(out_path, sep=",", index=False)
    return out_path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--in_dir", default=DEFAULT_IN_DIR, help="folder with the *.new.merged.filtered.SV.tsv files")
    ap.add_argument("-o", "--out_dir", default=str(DEFAULT_OUT_DIR), help="where to write the SVsig CSVs")
    ap.add_argument("cohorts", nargs="*", help="only convert these cohorts (default: all found)")
    args = ap.parse_args(argv)

    in_dir, out_dir = Path(args.in_dir).expanduser(), Path(args.out_dir).expanduser()
    if not in_dir.is_dir():
        print(f"Input directory not found: {in_dir}", file=sys.stderr)
        return 2  # distinct code so callers can fall back to already-converted files

    files = sorted(p for p in in_dir.iterdir() if p.name.endswith(SUFFIX))
    if args.cohorts:
        wanted = set(args.cohorts)
        files = [p for p in files if p.name.split(".")[0] in wanted]
        missing = wanted - {p.name.split(".")[0] for p in files}
        if missing:
            print(f"No input file for: {', '.join(sorted(missing))}", file=sys.stderr)
            return 1
    if not files:
        print(f"No *{SUFFIX} files in {in_dir}", file=sys.stderr)
        return 1

    for p in files:
        cohort = p.name.split(".")[0]
        out_path = convert_to_svsig(p, cohort, out_dir)
        print(f"[{cohort}] wrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
