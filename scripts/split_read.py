#!/usr/bin/env python3
"""
Split reads classified by Kraken2 into per-accession FASTQs.

• Reads whose accession is in the reference panel are written to
    <OUT_ROOT>/<accession>/<accession>_{R1,R2}.fq.gz
• Reads whose accession is not in the panel are written to
    <OUT_ROOT>/other_pathogens/<accession>/<accession>_{R1,R2}.fq.gz
"""

import argparse, gzip, json, sys, subprocess, tempfile, shutil, re, multiprocessing, os, random
from multiprocessing import Process
from pathlib import Path

def split_fastq(fq, mate, read2taxid, refs, acc2dir, acc2taxid, taxid2name,
                out_root, pigz_threads, batch_dir, recipients_by_taxid, tax2acc, canonical_policy):
    writers = {}
    processes = {}
    ext = ".fq.gz"
    with open_in(fq) as ih:
        it = iter(ih)
        while True:
            h = next(it, None)
            if h is None:
                break
            s = next(it, "")
            p = next(it, "")
            q = next(it, "")
            if not q:
                break

            rid = read_id(h)
            t = read2taxid.get(rid)
            if not t:
                continue

            # Panel routing via tree (own + descendants + immediate parent)
            acc_targets = list(recipients_by_taxid.get(t, []))

            # Fallback: normal split by direct taxid→accession mapping (single canonical)
            if not acc_targets:
                accs = tax2acc.get(t, [])
                if not accs:
                    continue
                acc_targets = [accs[0] if canonical_policy == "first" else accs[-1]]

            for acc in acc_targets:
                w, _, _ = writer_for(
                    acc, mate, ext, out_root, refs, acc2dir, acc2taxid, taxid2name,
                    writers, processes, pigz_threads, batch_dir
                )
                w.write(h); w.write(s); w.write(p); w.write(q)

    close_writers(writers, processes)
    if batch_dir is not None:
        batch_compress(writers, pigz_threads)

def get_read_classified_species(kraken_out, kraken_report):
    taxons_wepp_pathogens = get_taxid_of_pathogens_for_wepp()
    wepp_species_taxon_name_mapping = {}
    
    # 1. Get taxid→names from kraken_report
    taxid_to_names = {}
    with open(kraken_report, "r") as f:
        lines = [line.rstrip("\n") for line in f]

    i = 0
    while i < len(lines):
        line = lines[i]
        parts = line.split("\t")
        if len(parts) < 6:
            i += 1
            continue

        taxid = parts[4].strip()
        name = parts[5].strip()
        rank = parts[3].strip()

        # In case of genus, if the member species are in wepp_pathogens, then we need to map its taxon with those of wepp species
        if rank.startswith('G'):
            species_under_genus = []
            j = i + 1
            while j < len(lines):
                next_parts = lines[j].split("\t")
                if len(next_parts) >= 6 and next_parts[3].strip().startswith('S'):
                    species_under_genus.append((next_parts[4].strip(), next_parts[5].strip()))
                    j += 1
                else:
                    break
            
            wepp_species = []
            for sp_taxid, sp_name in species_under_genus:
                if sp_taxid in taxons_wepp_pathogens.keys():
                    wepp_species.append(sp_name)
                    wepp_species_taxon_name_mapping[sp_name] = sp_taxid

            if wepp_species:
                wepp_species_names = set(wepp_species)
                taxid_to_names.setdefault(taxid, set()).update(wepp_species_names)
                for sp_taxid, _ in species_under_genus:
                    taxid_to_names.setdefault(sp_taxid, set()).update(wepp_species_names)
            else:
                taxid_to_names.setdefault(taxid, set()).add(name)
                for sp_taxid, sp_name in species_under_genus:
                    taxid_to_names.setdefault(sp_taxid, set()).add(sp_name)

            i = j # Move to the line after the last species
        
        elif taxid:
            taxid_to_names.setdefault(taxid, set()).add(name)
            i += 1
        else:
            i += 1


    # 2. Get read→tax_id from kraken_report
    read_to_name = {}
    with open(kraken_out, "r") as f:
            for line in f:
                if line.startswith("C\t"):
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 3:
                        _, rid, taxid = parts[:3]

                        # Map read → species name(s) using the first table
                        if taxid in taxid_to_names:
                            names = taxid_to_names[taxid]
                            if len(names) > 1:
                                # In case of wepp species, prefer all those names
                                wepp_names = set(wepp_species_taxon_name_mapping.keys())
                                wepp_intersection = names.intersection(wepp_names)
                                if wepp_intersection:
                                    read_to_name[rid] = list(wepp_intersection)
                                else:
                                    # If no intersection, pick a random name
                                    read_to_name[rid] = [random.choice(list(names))]
                            else:
                                read_to_name[rid] = list(names)

    return read_to_name, wepp_species_taxon_name_mapping

def get_taxid_of_pathogens_for_wepp():
    added_taxons_path = "data/pathogens_for_wepp/added_taxons.tsv"
    added_taxons = {}
    if os.path.exists(added_taxons_path):
        with open(added_taxons_path, "r") as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        taxid, folder = parts
                        added_taxons[taxid] = folder
    return added_taxons
    
# ─────────────────────────── CLI ────────────────────────────
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("-k", "--kraken-out", required=True)
    ap.add_argument("-r", "--kraken-report", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2")
    ap.add_argument("-o", "--out-dir", required=True)
    ap.add_argument("-t","--threads", type=int, default=4)

    return ap.parse_args()

def main():
    a = parse_args()
    
    Path(a.out_dir).mkdir(parents=True, exist_ok=True)
    # Per-read classified taxid from Kraken output
    read2species = get_read_classified_species(a.kraken_out, a.kraken_report)
    if not read2species:
        sys.exit("No classified reads in Kraken output.")

    if a.r2:
        p = Process(
            target=run_split,
            args=(
                a.r2, "R2", read2species, refs, acc2dir, acc2taxid, taxid2name,
                a.out_dir, a.pigz_threads, batch_dir, recipients_by_taxid, tax2acc, a.canonical_policy
            )
        )
        p.start()

    # Always process R1 in this process
    run_split(
        a.r1, "R1", read2species, refs, acc2dir, acc2taxid, taxid2name,
        a.out_dir, a.pigz_threads, batch_dir, recipients_by_taxid, tax2acc, a.canonical_policy
    )

    # Wait for R2 if running
    if a.r2:
        p.join()

if __name__ == "__main__":
    main()