#!/usr/bin/env python3
import argparse, json, gzip, os, sys
from pathlib import Path

def total_bases_fastq(fq_path: Path) -> int:
    """
    Sum the lengths of ALL sequence lines in a FASTQ (gz or plain).
    Returns 0 if missing or unreadable.
    """
    if not fq_path.exists():
        return 0
    opener = gzip.open if fq_path.suffix == ".gz" or fq_path.name.endswith(".gz") else open
    total = 0
    try:
        with opener(fq_path, "rt", encoding="utf-8", errors="ignore") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    total += len(line.rstrip("\n\r"))
    except Exception:
        return 0
    return total

def genome_length(fa_dir: Path) -> int:
    candidates = []
    for pat in ("*.fa", "*.fasta", "*.fna", "*.FA", "*.FASTA", "*.FNA"):
        candidates += list(fa_dir.glob(pat))
    if not candidates:
        return 0
    fa = candidates[0]
    total = 0
    with open(fa, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                continue
            total += len(line.strip())
    return total

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--acc2classified", required=True, help="JSON: {acc: out_dir}")
    ap.add_argument("--out-root", required=True, help="Root where split FASTQs live")
    ap.add_argument("--pathogen-root", required=True, help="Root with data/pathogens_for_wepp/<out_dir>/*.fa")
    ap.add_argument("--min-depth", type=float, required=True)
    ap.add_argument("--out-json", required=True, help="Output JSON of eligible {acc: out_dir}")
    args = ap.parse_args()

    with open(args.acc2classified) as fh:
        acc2dir = json.load(fh)

    out_root = Path(args.out_root)
    pathogen_root = Path(args.pathogen_root)
    eligible = {}

    for acc, out_dir in acc2dir.items():
        acc_stem = acc.split(".", 1)[0]
        fq_r1 = out_root / out_dir / f"{acc_stem}_R1.fq.gz"

        total_bases = total_bases_fastq(fq_r1)
        fq_r2 = fq_r1.parent / fq_r1.name.replace("_R1.fq.gz", "_R2.fq.gz")
        if fq_r2.exists():
            total_bases += total_bases_fastq(fq_r2)

        genome_len = genome_length(pathogen_root / out_dir)

        depth = (total_bases / genome_len) if genome_len > 0 else 0.0

        if depth >= args.min_depth:
            eligible[acc] = out_dir

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_json, "w") as out:
        json.dump(eligible, out, indent=2)

if __name__ == "__main__":
    main()
