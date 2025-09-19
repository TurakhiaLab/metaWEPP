#!/usr/bin/env python3
"""
Helper script that reuses viral_usher's taxonomy→RefSeq interactive lookup,
then fetches the RefSeq FASTA, optionally adds it to a Kraken2 DB, and handles MAT.
"""

import sys
import os
import re
import logging
import subprocess
import shutil
import time
import argparse
import requests

# Require viral_usher so we reuse its HTTP/helper behavior exactly.
try:
    from viral_usher import ncbi_helper as vu_ncbi_helper
except Exception as e:
    print("Error: unable to import 'viral_usher.ncbi_helper'. Install viral-usher and ensure it's importable.")
    print(f"Import error: {e}")
    sys.exit(1)

NC_RE  = re.compile(r'NC_\d+(?:\.\d+)?', re.IGNORECASE)
GCF_RE = re.compile(r'GCF_\d+(?:\.\d+)?', re.IGNORECASE)

REQUESTS_TIMEOUT = 60
SLEEP_BETWEEN = 0.34

ROOT_PATHOGENS_DIR = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "../data", "pathogens_for_wepp")
)

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def efetch_fasta_by_accession(accession: str) -> str:
    r = requests.get(
        f"{NCBI_BASE}/efetch.fcgi",
        params={"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text"},
        timeout=REQUESTS_TIMEOUT,
    )
    r.raise_for_status()
    time.sleep(SLEEP_BETWEEN)
    txt = r.text
    if not txt.strip().startswith(">"):
        raise RuntimeError(f"FASTA not returned for {accession}")
    return txt

def sanitize_folder_name(name):
    return re.sub(r"[^A-Za-z0-9._-]", "", name)

def default_species_folder(species_name):
    return sanitize_folder_name(species_name.replace(" ", ""))

def ensure_folder(path):
    os.makedirs(path, exist_ok=True)
    return path

def save_fasta(path, accession, fasta_text):
    file_path = os.path.join(path, f"{accession}.fna")
    with open(file_path, "w", encoding="utf-8") as fh:
        fh.write(fasta_text)
    return file_path

def add_file_to_kraken(db_dir, fasta_path):
    cmd = ["k2", "add-to-library", "--db", db_dir, "--file", fasta_path]
    print("[CMD]", " ".join(cmd))
    subprocess.run(cmd, check=True)

def copy_mat(src, dst_dir):
    dst = os.path.join(dst_dir, os.path.basename(src))
    shutil.copy2(src, dst)
    return dst

def refseq_seems_in_db(db_dir, refseq_id):
    lib = os.path.join(db_dir, "library", "added")
    if not os.path.isdir(lib):
        return False
    prefix = f"{refseq_id}_"
    for fname in os.listdir(lib):
        if fname.startswith(prefix):
            return True
    return False

# --- viral_usher lookup functions (verbatim control flow) ---

def get_input(prompt):
    try:
        text = input(prompt)
    except (KeyboardInterrupt, EOFError):
        print("")
        sys.exit(1)
    if re.match(r"^[Qq](uit)?", text):
        print("Exiting...")
        sys.exit(0)
    return text

def prompt_int_choice(prompt, min_value=1, max_value=None):
    while True:
        try:
            choice = get_input(f"{prompt} [{min_value}]: ")
            choice = int(choice) if choice else min_value
            if min_value <= choice and (max_value is None or choice <= max_value):
                return choice
            else:
                print(f"Please enter a number between {min_value} and {max_value}.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def prompt_taxonomy_id(species_term, tax_entries):
    print(f"\nFound the following matches for '{species_term}':")
    for idx, entry in enumerate(tax_entries):
        print(f"{idx+1}. {entry['sci_name']}")
    print(f"{len(tax_entries)+1}. None of the above, go back")
    return prompt_int_choice("Enter a number to choose an option", 1, len(tax_entries)+1)

def prompt_refseq_id(taxid, refseq_entries):
    print(f"\nFound the following RefSeq IDs for Taxonomy ID {taxid}:")
    for idx, entry in enumerate(refseq_entries):
        label = f"{entry['accession']}: {entry['title']}"
        if (entry["strain"] and entry["strain"] != "No strain" and entry["strain"] not in label):
            label += f" ({entry['strain']})"
        print(f"{idx+1}. {label}")
    print(f"{len(refseq_entries)+1}. None of the above, go back")
    return prompt_int_choice("Enter a number to choose an option", 1, len(refseq_entries)+1)

def get_species_taxonomy_refseq(ncbi):
    species_term = ""
    while not species_term:
        species_term = get_input("\nWhat is your virus of interest? ")
    return lookup_taxonomy_refseq(ncbi, species_term)

def lookup_taxonomy_refseq(ncbi, species_term):
    print(f"Looking up NCBI Taxonomy entries for '{species_term}'...")
    tax_entries = ncbi.get_taxonomy_entries('"' + species_term + '"')
    if tax_entries:
        return get_taxonomy_refseq(ncbi, species_term, tax_entries)
    else:
        print(f"\nNo matches found for '{species_term}'.")
        return get_species_taxonomy_refseq(ncbi)

def get_taxonomy_refseq(ncbi, species_term, tax_entries):
    choice = prompt_taxonomy_id(species_term, tax_entries)
    if choice > len(tax_entries):
        return get_species_taxonomy_refseq(ncbi)
    else:
        taxid = tax_entries[choice-1]["tax_id"]
        return lookup_refseq(ncbi, species_term, tax_entries, taxid)

def lookup_refseq(ncbi, species_term, tax_entries, taxid):
    print(f"Looking up RefSeqs associated with Taxonomy ID {taxid}...")
    refseq_entries = ncbi.get_refseqs_for_taxid(taxid)
    if refseq_entries:
        return get_refseq(ncbi, species_term, tax_entries, taxid, refseq_entries)
    else:
        print(f"No RefSeqs found for Taxonomy ID '{taxid}'.  Try a different match from species search.")
        return get_taxonomy_refseq(ncbi, species_term, tax_entries)

def get_refseq(ncbi, species_term, tax_entries, taxid, refseq_entries):
    choice = prompt_refseq_id(taxid, refseq_entries)
    if choice > len(refseq_entries):
        return get_taxonomy_refseq(ncbi, species_term, tax_entries)
    else:
        refseq_id = refseq_entries[choice-1]["accession"]
        print(f"Looking up NCBI Assembly accession for selected RefSeq '{refseq_id}'...")
        assembly_id = ncbi.get_assembly_acc_for_refseq_acc(refseq_id)
        if not assembly_id:
            print(f"Could not find assembly ID for RefSeq ID {refseq_id} -- can't download RefSeq.")
            print("Try choosing a different RefSeq.")
            return get_refseq(ncbi, species_term, tax_entries, taxid, refseq_entries)
        return species_term, taxid, refseq_id, assembly_id

# --- local workflow ---

def yes_no(prompt, default_yes=True):
    d = "Y/n" if default_yes else "y/N"
    while True:
        ans = (get_input(f"{prompt} [{d}]: ") or "").strip().lower()
        if not ans:
            return default_yes
        if ans in ("y", "yes"): return True
        if ans in ("n", "no"):  return False
        print("Please answer yes or no.")

def run_viral_usher(species, refseq_id, workdir):
    config_path = os.path.join(workdir, f"viral_usher_config_{refseq_id}.toml")
    cmd_init = [
        "viral_usher", "init",
        "--refseq", refseq_id,
        "--workdir", workdir,
        "--config", config_path
    ]
    cmd_build = ["viral_usher", "build", "--config", config_path]

    print("[CMD]", " ".join(cmd_init))
    subprocess.run(cmd_init, check=True)
    print("[CMD]", " ".join(cmd_build))
    subprocess.run(cmd_build, check=True)

    candidates = [
        os.path.join(workdir, "mat.pb.gz"),
        os.path.join(workdir, "optimized.pb.gz"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    print(f"[INFO] MAT will be located at:\n  {os.path.join(species_dir, 'mat.pb.gz')}")
    return None

def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    ap = argparse.ArgumentParser(
        description="Use viral_usher's exact lookup; then download FASTA, add to Kraken2, and handle MAT."
    )
    ap.add_argument("--db", required=True, help="Path to Kraken2 DB directory.")
    ap.add_argument("--skip-existing-files", action="store_true",
                    help="Skip FASTA download if the output file already exists.")
    args = ap.parse_args()

    ensure_folder(ROOT_PATHOGENS_DIR)

    while True:
        ncbi = vu_ncbi_helper.NcbiHelper()
        species, taxid, refseq_id, assembly_id = get_species_taxonomy_refseq(ncbi)

        default_name = default_species_folder(species)
        user_folder = get_input(f"Folder name inside pathogens_for_wepp [{default_name}]: ").strip() or default_name
        user_folder = sanitize_folder_name(user_folder)
        species_dir = ensure_folder(os.path.join(ROOT_PATHOGENS_DIR, user_folder))
        print(f"[INFO] Working directory: {species_dir}")

        fasta_path = os.path.join(species_dir, f"{refseq_id}.fna")
        if args.skip_existing_files and os.path.isfile(fasta_path):
            print(f"[SKIP] FASTA exists: {fasta_path}")
        else:
            print(f"[FETCH] {refseq_id} via NCBI efetch...")
            fasta = efetch_fasta_by_accession(refseq_id)
            fasta_path = save_fasta(species_dir, refseq_id, fasta)
            print(f"[OK] FASTA saved: {fasta_path}")

        if os.path.isfile(fasta_path):
            if refseq_seems_in_db(args.db, refseq_id):
                print(f"[INFO] {refseq_id} already present in DB library. Skipping add-to-library.")
            else:
                if yes_no("Add this FASTA to the Kraken2 DB now?", default_yes=True):
                    add_file_to_kraken(args.db, fasta_path)
                    if refseq_seems_in_db(args.db, refseq_id):
                        print(f"[INFO] Added to Kraken2 DB: {refseq_id} (confirmed)")
                    else:
                        print(f"[INFO] Added to Kraken2 DB: {refseq_id}")
                else:
                    print("[INFO] Skipped adding to Kraken2 DB.")
        else:
            print("[WARN] FASTA missing; cannot add to Kraken2 DB.")

        print("\n(Optional) Provide path to an existing MAT (Mutation-Annotated Tree).")
        print("If omitted, a MAT will be built with viral_usher. Ensure the MAT reference matches the RefSeq.")
        mat_path_in = (get_input("Path to MAT (or press Enter to build): ") or "").strip()

        if mat_path_in:
            if not os.path.isfile(mat_path_in):
                print(f"[WARN] MAT not found: {mat_path_in}")
            else:
                copied = copy_mat(mat_path_in, species_dir)
                print(f"[OK] MAT copied: {copied}")
        else:
            print("[INFO] Building MAT via viral_usher (init → build)…")
            print(f"[INFO] MAT will be at:\n  {os.path.join(species_dir, 'mat.pb.gz')}")
            try:
                built_mat = run_viral_usher(species, refseq_id, species_dir)
                if built_mat:
                    final_mat = copy_mat(built_mat, species_dir) if os.path.dirname(built_mat) != species_dir else built_mat
                    print(f"[OK] MAT ready: {final_mat}")
                else:
                    print("[WARN] viral_usher finished but MAT file not found. Check logs/outputs.")
            except FileNotFoundError:
                print("[ERROR] 'viral_usher' executable not found on PATH.")
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] viral_usher failed with exit code {e.returncode}.")

        if not yes_no("\nProcess another virus?", default_yes=False):
            break

if __name__ == "__main__":
    main()