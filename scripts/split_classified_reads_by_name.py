#!/usr/bin/env python3
import argparse
import gzip
import os
import sys

def parse_args():
    p = argparse.ArgumentParser(
        description="Split paired FASTQs into per‑taxid sub‑folders based on Kraken2 output"
    )
    p.add_argument(
        "-k", "--kraken-out", required=True,
        help="Kraken2 --output file (tab-delimited, 'C' lines only)"
    )
    p.add_argument(
        "--r1", required=True,
        help="Raw R1 FASTQ file (can be .fastq or .fastq.gz)"
    )
    p.add_argument(
        "--r2", required=True,
        help="Raw R2 FASTQ file (can be .fastq or .fastq.gz)"
    )
    p.add_argument(
        "-o", "--out-dir", required=True,
        help="Directory in which to write per‑taxid FASTQs"
    )

    p.add_argument(
    "--report", required=True,
    help="Kraken2 report file used to name taxid fq files"
  )

    return p.parse_args()

def open_infile(fn):
    if fn.endswith(".gz"):
        return gzip.open(fn, "rt")
    else:
        return open(fn, "r")

def open_outfile(fn):
    if fn.endswith(".gz"):
        return gzip.open(fn, "wt")
    else:
        return open(fn, "w")

def parse_kraken(fn):
    mapping = {}
    with open_infile(fn) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] != "C":
                continue
            read_id = parts[1]
            taxid   = parts[2]
            mapping[read_id] = taxid
    return mapping

def split_reads(fastq_path, mapping, out_dir, mate_label, taxid_to_name):
    writers = {}
    ext = ".fq.gz" if fastq_path.endswith(".gz") else ".fq"

    with open_infile(fastq_path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq   = fh.readline()
            plus  = fh.readline()
            qual  = fh.readline()
            if not qual:
                break

            hdr = header.strip().split()[0]
            if hdr.startswith("@"):
                hdr = hdr[1:]
            if hdr.endswith("/1") or hdr.endswith("/2"):
                base_id = hdr[:-2]
            else:
                base_id = hdr

            taxid = mapping.get(base_id)
            if taxid is None:
                continue

            # ensure taxid sub‑directory exists
            taxon_dir = os.path.join(out_dir, taxid)
            os.makedirs(taxon_dir, exist_ok=True)

            if taxid not in writers:
                out_path = os.path.join(taxon_dir, f"{taxid}_{mate_label}{ext}")
                writers[taxid] = open_outfile(out_path)

            w = writers[taxid]
            w.write(header)
            w.write(seq)
            w.write(plus)
            w.write(qual)

    for w in writers.values():
        w.close()

  def parse_report(report_fn):
    taxid_to_name = {}
    with open(report_fn) as f:
        for line in f:
            parts = line.strip().split(None, 5)  # whitespace-split, at most 6 parts
            if len(parts) < 6:
                continue
            taxid = parts[4]
            name = parts[5].strip()
            name_clean = name.replace(" ", "_")
            taxid_to_name[taxid] = name_clean
    return taxid_to_name

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    mapping = parse_kraken(args.kraken_out)
    taxid_to_name = parse_report(args.report)

    split_reads(args.r1, mapping, args.out_dir, "R1", taxid_to_name)
    split_reads(args.r2, mapping, args.out_dir, "R2", taxid_to_name)

if __name__ == "__main__":
    main()
