#!/usr/bin/env python3
import argparse, gzip, os, sys

def parse_args():
    p = argparse.ArgumentParser(
        description="Split paired FASTQs into per-taxid (now per-accession) folders."
    )
    p.add_argument("-k", "--kraken-out", required=True,
                   help="Kraken2 --output file (only ‘C’ lines are read)")
    p.add_argument("-m", "--mapping", required=True,
                   help="2-column file: <accession>\t<taxid>")
    p.add_argument("--r1", required=True, help="R1 FASTQ (.fq or .fq.gz)")
    p.add_argument("--r2", required=True, help="R2 FASTQ (.fq or .fq.gz)")
    p.add_argument("-o", "--out-dir", required=True,
                   help="Root directory that will receive the per-accession folders")
    return p.parse_args()

def open_in(fn):  return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)
def open_out(fn): return gzip.open(fn, "wt") if fn.endswith(".gz") else open(fn)

def load_taxid_map(path):
    d = {}
    with open(path) as f:
        for acc, tax in (l.rstrip().split(None, 1) for l in f if l.strip()):
            d.setdefault(tax, acc)          # first accession wins if duplicates
    return d                                # taxid → accession

def load_classifications(path):
    m = {}
    with open_in(path) as f:
        for line in f:
            if not line.startswith("C\t"): continue
            _, read_id, taxid, *_ = line.rstrip("\n").split("\t")
            m[read_id] = taxid
    return m                                # read_id → taxid

def split_reads(fq, read2tax, tax2acc, out_root, mate):
    writers, ext = {}, ".fq.gz" if fq.endswith(".gz") else ".fq"
    with open_in(fq) as ih:
        while True:
            h = ih.readline()
            if not h: break
            s, p, q = ih.readline(), ih.readline(), ih.readline()
            rid = h.strip().split()[0].lstrip("@").rstrip("/12")
            taxid = read2tax.get(rid)
            acc   = tax2acc.get(taxid)
            if acc is None: continue

            folder = os.path.join(out_root, acc)
            os.makedirs(folder, exist_ok=True)
            if acc not in writers:
                prefix = acc.split(".")[0]
                fn = f"{folder}/{prefix}_{mate}{ext}"
                writers[acc] = open_out(fn)
            w = writers[acc]
            w.writelines([h, s, p, q])

    for w in writers.values(): w.close()

def main():
    a          = parse_args()
    tax2acc    = load_taxid_map(a.mapping)
    read2taxid = load_classifications(a.kraken_out)
    os.makedirs(a.out_dir, exist_ok=True)
    split_reads(a.r1, read2taxid, tax2acc, a.out_dir, "R1")
    split_reads(a.r2, read2taxid, tax2acc, a.out_dir, "R2")

if __name__ == "__main__":
    main()
