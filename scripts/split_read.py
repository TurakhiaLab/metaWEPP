#!/usr/bin/env python3
"""
Split reads classified by Kraken2 into per-accession FASTQs.

• Reads whose accession *is* in the reference panel are written to
    <OUT_ROOT>/<accession>/<accession>_{R1,R2}.fq.gz
• Reads whose accession is *not* in the panel are written to
    <OUT_ROOT>/other_pathogens/<accession>/<accession>_{R1,R2}.fq.gz
"""

import argparse, gzip, os, sys
from pathlib import Path

# ────────────────────────── helpers ──────────────────────────
def open_in(fn):
    fn = str(fn)
    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn, "r")

def open_out(fn):
    fn = str(fn)
    Path(fn).parent.mkdir(parents=True, exist_ok=True)
    return gzip.open(fn, "wt") if fn.endswith(".gz") else open(fn, "w")
def read_id(header):
    """
    Return read name without leading @ and without /1 or /2 suffix,
    so it matches what Kraken writes.
    """
    return header.split()[0].lstrip("@").rstrip("/12")

def load_taxid_map(path):
    d = {}
    with open(path) as f:
        for line in f:
            if line.strip():
                acc, taxid = line.split(None, 1)
                d[taxid.strip()] = acc.strip()
    return d                                            # taxid → accession

def load_classifications(path, tax2acc):
    m = {}
    with open_in(path) as f:
        for line in f:
            if line.startswith("C\t"):
                _, rid, taxid, *_ = line.rstrip("\n").split("\t")
                acc = tax2acc.get(taxid)
                if acc:
                    m[rid] = acc
    return m                                           # readID → accession

def load_ref_accessions(s):
    """`s` may be a path to a file or a comma-separated string."""
    if s is None:
        return set()
    p = Path(s)
    if p.exists():
        return {l.strip() for l in p.read_text().splitlines() if l.strip()}
    # otherwise treat as comma-list
    return {x.strip() for x in s.split(",") if x.strip()}

# decide where to write
def out_root_for(acc, ref_accessions, base_out):
    if acc in ref_accessions:
        return Path(base_out) / acc
    return Path(base_out) / "other_pathogens" / acc

# ───────────────────────── splitting ─────────────────────────
def split_fastq(fq, mate, read2acc, refs, out_root):
    writers, seen_accs = {}, set()
    ext = ".fq.gz" if fq.endswith(".gz") else ".fq"

    with open_in(fq) as ih:
        while True:
            h = ih.readline()
            if not h:
                break
            s, p, q = ih.readline(), ih.readline(), ih.readline()
            acc = read2acc.get(read_id(h))
            if not acc:
                continue

            seen_accs.add(acc)
            root = out_root_for(acc, refs, out_root)
            if acc not in writers:
                fn = root / f"{acc.split('.')[0]}_{mate}{ext}"
                writers[acc] = open_out(fn)
            writers[acc].writelines([h, s, p, q])

    for w in writers.values(): w.close()

    # ensure stub for mate that got zero reads
    for acc in seen_accs:
        stub = out_root_for(acc, refs, out_root) / f"{acc.split('.')[0]}_{mate}{ext}"
        if not stub.exists():
            with open_out(stub): pass

# ─────────────────────────── CLI ────────────────────────────
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("-k", "--kraken-out", required=True)
    ap.add_argument("-m", "--mapping",    required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2")
    ap.add_argument("-o", "--out-dir", required=True)
    ap.add_argument("--ref-accessions",
                    help="File or comma-list of accessions in the reference panel")
    return ap.parse_args()

def main():
    a = parse_args()
    Path(a.out_dir).mkdir(parents=True, exist_ok=True)

    refs      = load_ref_accessions(a.ref_accessions)
    tax2acc   = load_taxid_map(a.mapping)
    read2acc  = load_classifications(a.kraken_out, tax2acc)

    if not read2acc:
        sys.exit("No classified reads matched the mapping file – nothing to split.")

    split_fastq(a.r1, "R1", read2acc, refs, a.out_dir)
    if a.r2:
        split_fastq(a.r2, "R2", read2acc, refs, a.out_dir)

if __name__ == "__main__":
    main()
