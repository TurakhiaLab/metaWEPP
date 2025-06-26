#!/usr/bin/env python3
import argparse, gzip, os, sys
from pathlib import Path

# ───────────── IO helpers ──────────────
def open_in(fn):   return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn, "r")
def open_out(fn):  Path(fn).parent.mkdir(parents=True, exist_ok=True); return gzip.open(fn, "wt") if fn.endswith(".gz") else open(fn, "w")

def read_id(header):
    """Return the read name without leading @ (keeps /1 or /2)."""
    return header.split()[0].lstrip("@")

# ─────── load mapping & classifications ───────
def load_taxid_map(path):
    tax2acc = {}
    with open(path) as f:
        for line in f:
            if line.strip():
                acc, taxid = line.split(None, 1)
                tax2acc[taxid.strip()] = acc.strip()
    return tax2acc                                # taxid → accession

def load_classifications(path, tax2acc):
    read2acc = {}
    with open_in(path) as f:
        for line in f:
            if not line.startswith("C\t"):
                continue
            _c, rid, taxid, *_ = line.rstrip("\n").split("\t")
            acc = tax2acc.get(taxid)
            if acc:
                read2acc[rid] = acc
    return read2acc                               # readID → accession


def split_fastq(fq, mate_label, read2acc, out_root):
    writers = {}
    ext = ".fq.gz" if fq.endswith(".gz") else ".fq"
    processed = 0
    seen_accs = set()                         # ← NEW

    with open_in(fq) as ih:
        while True:
            h = ih.readline()
            if not h:
                break
            s, p, q = ih.readline(), ih.readline(), ih.readline()

            acc = read2acc.get(read_id(h))
            if not acc:
                continue

            seen_accs.add(acc)                # ← track real hits
            if acc not in writers:
                prefix = acc.split(".")[0]
                fn = os.path.join(out_root, acc, f"{prefix}_{mate_label}{ext}")
                writers[acc] = open_out(fn)
            writers[acc].writelines([h, s, p, q])

            processed += 1
            if processed % 1_000_000 == 0:
                print(f"[{mate_label}] wrote {processed:,} reads", file=sys.stderr)

    # close all real writers
    for w in writers.values():
        w.close()

    # ───── create placeholders for mate-label that got zero reads ─────
    for acc in seen_accs:                     # only those that appeared
        prefix = acc.split(".")[0]
        fn = os.path.join(out_root, acc, f"{prefix}_{mate_label}{ext}")
        if not os.path.exists(fn):            # still missing? create stub
            with (gzip.open(fn, "wt") if fn.endswith(".gz") else open(fn, "w")):
                pass


# ───────── CLI / main ─────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-k", "--kraken-out", required=True)
    p.add_argument("-m", "--mapping",    required=True)
    p.add_argument("--r1", required=True)
    p.add_argument("--r2")
    p.add_argument("-o", "--out-dir", required=True)
    return p.parse_args()

def main():
    a = parse_args()
    Path(a.out_dir).mkdir(parents=True, exist_ok=True)

    tax2acc  = load_taxid_map(a.mapping)
    read2acc = load_classifications(a.kraken_out, tax2acc)

    if not read2acc:
        sys.exit("No classified reads matched the mapping file – nothing to split.")

    split_fastq(a.r1, "R1", read2acc, a.out_dir)
    if a.r2:
        split_fastq(a.r2, "R2", read2acc, a.out_dir)

if __name__ == "__main__":
    main()

