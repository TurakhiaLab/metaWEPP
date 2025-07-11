#!/usr/bin/env python3
"""
Generate acc2classified_dir.json

Keeps only accessions that
  • have ≥1 read classified by Kraken (status ‘C’)
  • are present in config/acc2dirname.json     (your reference panel)

The result maps accession  →  directory-name (e.g.  'NC_038235.1' → 'rsv_a').
"""

import argparse, gzip, json, sys
from pathlib import Path
from collections import defaultdict

# ── helpers ────────────────────────────────────────────────────────────────
def open_in(fn):
    return gzip.open(fn, "rt") if str(fn).endswith(".gz") else open(fn, "r")

def load_tax_map(path):
    """
    seqid2taxid.map  ➜  dict[taxid] -> list[accession]
    (one taxid can have several accessions)
    """
    tax2acc = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.strip():
                acc, taxid = line.split(None, 1)
                tax2acc[taxid.strip()].append(acc.strip())
    return tax2acc

def classified_taxids(kraken_out):
    """Return the set of taxids that appear in ‘C’ rows of kraken output."""
    taxids = set()
    with open_in(kraken_out) as f:
        for line in f:
            if line.startswith("C\t"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    taxids.add(parts[2])
    return taxids

# ── main ───────────────────────────────────────────────────────────────────
def cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kraken-out", required=True)
    ap.add_argument("--mapping",    required=True, help="seqid2taxid.map")
    ap.add_argument("--acc2dir",    required=True, help="config/acc2dirname.json")
    ap.add_argument("-o", "--outfile",
                    default="config/acc2classified_dir.json")
    return ap.parse_args()

def main():
    a          = cli()
    tax2acc    = load_tax_map(a.mapping)
    good_taxid = classified_taxids(a.kraken_out)
    acc2dir    = json.load(open(a.acc2dir))

    result = {}
    for taxid in good_taxid:
        for acc in tax2acc.get(taxid, []):
            if acc in acc2dir:                    # only keep panel accessions
                result[acc] = acc2dir[acc]

    # write
    Path(a.outfile).parent.mkdir(parents=True, exist_ok=True)
    with open(a.outfile, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"Wrote {a.outfile} with {len(result)} accessions.")

if __name__ == "__main__":
    main()
