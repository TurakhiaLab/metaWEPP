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
import glob

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


def get_taxid_for_accession(accession: str) -> str:
    """
    Use NCBI E-utils to get the taxonomy ID for a given RefSeq accession.
    """
    try:
        r = requests.get(
            f"{NCBI_BASE}/esummary.fcgi",
            params={"db": "nuccore", "id": accession, "retmode": "json"},
            timeout=REQUESTS_TIMEOUT,
        )
        r.raise_for_status()
        time.sleep(SLEEP_BETWEEN)
        data = r.json()
        result = data.get("result", {})
        uids = result.get("uids", [])
        if not uids:
            raise RuntimeError(f"No result found for {accession}")
        taxid = result[uids[0]].get("taxid")
        if not taxid:
            raise RuntimeError(f"No taxid found for {accession}")
        return str(taxid)
    except (requests.RequestException, RuntimeError, KeyError, IndexError) as e:
        logging.error(f"Failed to get taxid for {accession}: {e}")
        return None


def sanitize_folder_name(name):
    return re.sub(r"[^A-Za-z0-9._-]", "", name)

def default_species_folder(species_name):
    return sanitize_folder_name(species_name.replace(" ", "_"))

def ensure_folder(path):
    os.makedirs(path, exist_ok=True)
    return path

def save_fasta(path, accession, fasta_text):
    file_path = os.path.join(path, f"{accession}.fna")
    with open(file_path, "w", encoding="utf-8") as fh:
        fh.write(fasta_text)
    return file_path

def add_file_to_kraken(db_dir, fasta_path):
    cmd = ["kraken2", "add-to-library", "--db", db_dir, "--file", fasta_path]
    print("[CMD]", " ".join(cmd))
    subprocess.run(cmd, check=True)

def build_kraken_db(db_dir, threads=8):
    cmd = ["kraken2", "build", "--db", db_dir, f"--threads={threads}"]
    print("[CMD]", " ".join(cmd))
    subprocess.run(cmd, check=True)

def copy_viz_and_jsonl(workdir, pathogen_dir):
    os.makedirs(pathogen_dir, exist_ok=True)

    # Copy viz.pb.gz
    viz_src = os.path.join(workdir, "viz.pb.gz")
    if os.path.isfile(viz_src):
        dst = os.path.join(pathogen_dir, "viz.pb.gz")
        shutil.copy2(viz_src, dst)
        print(f"[INFO] viz.pb.gz copied → {dst}")
    else:
        print("\n[WARN] viz.pb.gz not found in workdir\n")

    # Copy all *.jsonl or *.jsonl.gz files
    jsonls = glob.glob(os.path.join(workdir, "*.jsonl")) + glob.glob(os.path.join(workdir, "*.jsonl.gz"))
    if jsonls:
        for src in jsonls:
            dst = os.path.join(pathogen_dir, os.path.basename(src))
            shutil.copy2(src, dst)
        print(f"[INFO] Copied {len(jsonls)} JSONL files → {pathogen_dir}")
    else:
        print("\n[WARN] No JSONL files found in workdir\n")


def refseq_seems_in_db(db_dir, refseq_id):
    """
    Checks if a RefSeq ID appears to be in the Kraken2 database.
    It does this by fetching the taxonomy ID for the RefSeq accession,
    and then checking if that taxid is mapped to any sequence in the Kraken2 DB's
    seqid2taxid.map file. This indicates that the taxon is present in the database.
    As a fallback, it also checks for the refseq in the 'library/added' folder.
    """
    taxid = get_taxid_for_accession(refseq_id)
    if not taxid:
        print(f"\n[WARN] Could not determine taxid for {refseq_id}. Assuming it's not in the DB to be safe.")
        return False

    # 1. Check if the taxon is in the built database by checking seqid2taxid.map
    seqid2taxid_path = os.path.join(db_dir, "seqid2taxid.map")
    if os.path.isfile(seqid2taxid_path):
        # Use grep for a much faster check on potentially huge map files
        cmd = ["grep", "-q", "-m", "1", f"\t{taxid}$\", seqid2taxid_path"]
        result = subprocess.run(cmd)
        if result.returncode == 0:
            print(f"[INFO] Found taxid {taxid} (for {refseq_id}) in the Kraken2 DB.")
            return True

    # 2. Fallback: check library/added for files that have been added but not yet built
    lib_added_dir = os.path.join(db_dir, "library", "added")
    if os.path.isdir(lib_added_dir):
        for fname in os.listdir(lib_added_dir):
            if refseq_id in fname:
                print(f"[INFO] Found {refseq_id} in the 'added' library (may not be built into DB yet).")
                return True
    
    print(f"[INFO] Neither taxid {taxid} nor RefSeq {refseq_id} found in Kraken DB '{db_dir}'.")
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
        assembly_id = None
        try:
            assembly_id = ncbi.get_assembly_acc_for_refseq_acc(refseq_id)
        except Exception as assembly_err:
            logging.warning(
                "Assembly lookup failed for RefSeq %s (%s); continuing with RefSeq only.",
                refseq_id,
                assembly_err,
            )
        if not assembly_id:
            logging.warning(
                "Assembly ID unavailable for RefSeq %s; continuing with RefSeq only.",
                refseq_id,
            )
        return species_term, taxid, refseq_id, assembly_id

# --- local workflow ---

def yes_no(prompt):
    d = "y/N"
    while True:
        ans = (get_input(f"{prompt} [{d}]: ") or "").strip().lower()
        if not ans:
            return False
        if ans in ("y", "yes"): return True
        if ans in ("n", "no"):  return False
        print("Please answer yes or no.")


def get_mat_interactively(species_dir: str) -> str:
    """
    Ask user for MAT path. If provided, copy it.
    If user presses Enter, return empty string to signal `viral_usher` build.
    """
    while True:
        mat_path = get_input("Path to MAT file (or press Enter to build with viral_usher): ").strip()
        if not mat_path:
            return ""  # User wants to build

        if not os.path.isfile(mat_path):
            print("File not found; try again.")
            continue
        if not (mat_path.endswith(".pb") or mat_path.endswith(".pb.gz")):
            print("MAT file must end with .pb or .pb.gz; try again.")
            continue

        dest = os.path.join(species_dir, os.path.basename(mat_path))
        if os.path.abspath(mat_path) == os.path.abspath(dest):
            print("[INFO] MAT already present in species directory; using existing file.")
            return dest
        try:
            shutil.copy2(mat_path, dest)
            print(f"[INFO] Copied MAT → {dest}")
            return dest
        except Exception as copy_err:
            print(f"[ERROR] Failed to copy MAT to {dest}: {copy_err}")
            continue

def mat_present(directory: str) -> bool:
    """Check whether the directory contains a MAT (.pb or .pb.gz)."""
    for pattern in ("*.pb", "*.pb.gz"):
        if glob.glob(os.path.join(directory, pattern)):
            return True
    return False


def any_mat_available(root: str) -> bool:
    """Return True if any pathogen directory under root has a MAT."""
    if not os.path.isdir(root):
        return False
    for entry in os.scandir(root):
        if entry.is_dir() and mat_present(entry.path):
            return True
    return False

def run_viral_usher_interactive(species, refseq_id, workdir):
    """
    Run viral_usher in the given workdir, letting it prompt for arguments interactively.
    """
    os.makedirs(workdir, exist_ok=True)
    config_path = os.path.join(workdir, f"viral_usher_config_{refseq_id}.toml")

    print("\n--- viral_usher setup ---")
    print("Running 'viral_usher init' interactively.")
    print(f"\n\n!!![IMPORTANT] When prompted for a 'workdir' to download and build the MAT, please enter the following path to ensure outputs are saved correctly:")
    print(f"  {workdir}\n")


    args = ["--refseq", refseq_id, "--config", config_path]
    # By not passing --workdir, viral_usher init remains interactive for most options.

    cmd_init = ["viral_usher", "init"] + args
    cmd_build = ["viral_usher", "build", "--config", config_path]

    print("\nPlease follow the prompts from 'viral_usher init':\n")
    #print("[CMD]", " ".join(cmd_init), "\n")
    subprocess.run(cmd_init, check=True)

    print("[CMD]", " ".join(cmd_build), "\n")
    subprocess.run(cmd_build, check=True)

    viz_path = os.path.join(workdir, "viz.pb.gz")
    if os.path.isfile(viz_path):
        print(f"[INFO] Found viz.pb.gz → {viz_path}")
        return viz_path
    else:
        print("\n[WARN] viz.pb.gz not found after build")
        return None

def add_pathogens_workflow(args):
    try:
        from viral_usher import ncbi_helper as vu_ncbi_helper
    except ImportError as e:
        print("Error: unable to import 'viral_usher.ncbi_helper'. Install viral-usher and ensure it's importable.", file=sys.stderr)
        print(f"Import error: {e}", file=sys.stderr)
        sys.exit(1)

    while True:
        ncbi = vu_ncbi_helper.NcbiHelper()
        species, taxid, refseq_id, assembly_id = get_species_taxonomy_refseq(ncbi)

        # Create folder name from species name by replacing spaces with underscores
        folder_name = species.replace(" ", "_")
        print(f"[INFO] The folder for '{species}' will be named '{folder_name}'.")
        species_dir = os.path.join(ROOT_PATHOGENS_DIR, folder_name)

        if os.path.isdir(species_dir):
            print(f"\n[WARN] Folder already exists: {species_dir}\n")
        else:
            os.makedirs(species_dir)

        # --- Fetch/save FASTA ---
        fasta_path = os.path.join(species_dir, f"{refseq_id}.fna")
        print(f"[FETCH] {refseq_id} via NCBI efetch...")
        fasta = efetch_fasta_by_accession(refseq_id)
        fasta_path = save_fasta(species_dir, refseq_id, fasta)
        print(f"[INFO] FASTA saved: {fasta_path}")

        # --- Auto-add FASTA to Kraken2 DB if not already there ---
        if args.db:
            if os.path.isfile(fasta_path):
                if refseq_seems_in_db(args.db, refseq_id):
                    print(f"[INFO] {refseq_id} or its taxon already present in DB. Skipping add-and-build.")
                else:
                    print(f"[INFO] Pathogen {refseq_id} is absent. Adding FASTA to the Kraken2 DB...")
                    add_file_to_kraken(args.db, fasta_path)
                    print(f"[INFO] Building Kraken2 DB to incorporate {refseq_id}...")
                    build_kraken_db(args.db, threads=8)
                    # Confirm after build
                    if refseq_seems_in_db(args.db, refseq_id):
                        print(f"[INFO] Successfully added and built DB with {refseq_id} (confirmed).")
                    else:
                        print(f"\n[WARN] Added {refseq_id} and rebuilt DB, but confirmation check failed. Please inspect the DB manually.\n")
            else:
                print("\n[WARN] FASTA missing; cannot add to Kraken2 DB.\n")

        existing_mat = get_mat_interactively(species_dir)

        # --- Build with viral_usher and copy viz/jsonl outputs ---
        if existing_mat:
            print(f"[INFO] Using provided MAT: {existing_mat}")
        else:
            print("[INFO] Running viral_usher (init → build)…")
            try:
                build_dir = os.path.join(species_dir, "viral_usher_build")
                viz_path = run_viral_usher_interactive(species, refseq_id, build_dir)
                # Copy only viz.pb.gz and JSONL outputs into the pathogen folder
                copy_viz_and_jsonl(build_dir, species_dir)

                if viz_path:
                    print(f"[INFO] Final viz.pb.gz available at: {viz_path}")
                else:
                    print("\n[WARN] viz.pb.gz not found after viral_usher build.\n")

                try:
                    shutil.rmtree(build_dir)
                    print(f"[INFO] Cleaned build dir: {build_dir}")
                except Exception as cleanup_err:
                    print(f"\n[WARN] Could not remove build dir {build_dir}: {cleanup_err}\n")

            except FileNotFoundError:
                print("[ERROR] 'viral_usher' executable not found on PATH.")
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] viral_usher failed with exit code {e.returncode}.")

        if not mat_present(species_dir):
            print(
                "[ERROR] No MAT (.pb/.pb.gz) found in this pathogen folder. "
                "Please provide a MAT before continuing."
            )
            continue

        if not yes_no("Would you like to add another pathogen?"):
            print("[INFO] Finished adding pathogens.")
            break

def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    ap = argparse.ArgumentParser(
        description="Use viral_usher's lookup; then download FASTA, auto-add to Kraken2, and copy viz/jsonl outputs."
    )
    ap.add_argument("--db", help="Path to Kraken2 DB directory.")
    args = ap.parse_args()

    ensure_folder(ROOT_PATHOGENS_DIR)

    if yes_no("Do you want to add a new species for haplotype-level analysis?"):
        add_pathogens_workflow(args)

    if not any_mat_available(ROOT_PATHOGENS_DIR):
        print(
            "\n\n[WARN] No MAT (.pb/.pb.gz) found in data/pathogens_for_wepp. "
            "At least one pathogen with a MAT is required for the pipeline to run."
        )
    else:
        print("\n[INFO] Setup complete. Pathogen data is ready.")


if __name__ == "__main__":
    main()
