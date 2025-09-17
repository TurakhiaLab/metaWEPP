#!/usr/bin/env python3
import os, re, sys, time, argparse, subprocess, requests, logging, shutil

# -------- Config --------
NC_RE  = re.compile(r'NC_\d+(?:\.\d+)?', re.IGNORECASE)
GCF_RE = re.compile(r'GCF_\d+(?:\.\d+)?', re.IGNORECASE)

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_TOOL = "wepp-nc-fetch"
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "you@example.org")  # set env var for etiquette
REQUESTS_TIMEOUT = 60
SLEEP_BETWEEN = 0.34  # ~3 req/s (NCBI-friendly)

ROOT_PATHOGENS_DIR = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "data", "pathogens_for_wepp")
)

# -------- I/O helpers (viral_usher-like) --------
def get_input(prompt):
    try:
        text = input(prompt)
    except (KeyboardInterrupt, EOFError):
        print("")
        sys.exit(1)
    if re.match(r"^[Qq](uit)?", (text or "")):
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
            print(f"Please enter a number between {min_value} and {max_value}.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def yes_no(prompt, default_yes=True):
    d = "Y/n" if default_yes else "y/N"
    while True:
        ans = (get_input(f"{prompt} [{d}]: ") or "").strip().lower()
        if not ans:
            return default_yes
        if ans in ("y", "yes"): return True
        if ans in ("n", "no"):  return False
        print("Please answer yes or no.")

# -------- NCBI HTTP helpers --------
def _get(url, params):
    params = dict(params)
    params.setdefault("tool", NCBI_TOOL)
    params.setdefault("email", NCBI_EMAIL)
    r = requests.get(url, params=params, timeout=REQUESTS_TIMEOUT)
    r.raise_for_status()
    time.sleep(SLEEP_BETWEEN)
    return r

def ncbi_search_taxonomy(term, retmax=50):
    """
    Broadened taxonomy search:
      - search [All Names] (scientific + synonyms)
      - try simple variants (spaces<->hyphens)
      - for single tokens, also try a prefix wildcard (e.g., 'sars*')
    Returns a list of {'tax_id', 'sci_name'} up to retmax.
    """
    term = (term or "").strip()
    if not term:
        return []

    # Build query variants
    variants = {term}
    variants.add(term.replace("-", " "))
    variants.add(term.replace(" ", "-"))

    queries = set()
    for v in variants:
        v = v.strip()
        if not v:
            continue
        # Search across all names (scientific, common, synonyms)
        queries.add(f'{v}[All Names]')
        # If it's a single "word" (no spaces/hyphens), try a prefix wildcard
        if re.match(r"^[A-Za-z0-9]+$", v):
            queries.add(f'{v}*[All Names]')

    # Collect IDs from all queries
    idset = []
    seen_ids = set()
    for q in queries:
        r = _get(f"{NCBI_BASE}/esearch.fcgi", {
            "db": "taxonomy",
            "term": q,
            "retmax": retmax,
            "retmode": "json",
        })
        ids = r.json().get("esearchresult", {}).get("idlist", []) or []
        for tid in ids:
            if tid not in seen_ids:
                seen_ids.add(tid)
                idset.append(tid)
        # Stop early if we've collected enough
        if len(idset) >= retmax:
            break

    if not idset:
        return []

    # Summarize (limit to retmax)
    idlist = idset[:retmax]
    r2 = _get(f"{NCBI_BASE}/esummary.fcgi", {
        "db": "taxonomy",
        "id": ",".join(idlist),
        "retmode": "json",
    })
    res = r2.json().get("result", {})
    out = []
    for tid in res.get("uids", []):
        rec = res.get(tid, {}) or {}
        sci = rec.get("scientificname") or rec.get("title") or f"taxid {tid}"
        out.append({"tax_id": tid, "sci_name": sci})
    return out

def ncbi_refseqs_for_taxid(taxid, retmax=20):
    term = f"txid{taxid}[Organism:exp] AND refseq[filter]"
    r = _get(f"{NCBI_BASE}/esearch.fcgi", {
        "db": "nuccore",
        "term": term,
        "retmax": retmax,
        "retmode": "json",
    })
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    r2 = _get(f"{NCBI_BASE}/esummary.fcgi", {
        "db": "nuccore",
        "id": ",".join(ids),
        "retmode": "json",
    })
    res = r2.json().get("result", {})
    out = []
    for uid in res.get("uids", []):
        rec = res.get(uid, {})
        accver = rec.get("accessionversion") or rec.get("caption") or ""
        title  = rec.get("title") or ""
        strain = None
        m = re.search(r"strain ([^,;]+)", title, re.IGNORECASE)
        if m: strain = m.group(1).strip()
        if NC_RE.match(accver):
            out.append({"accession": accver.upper(), "title": title, "strain": strain})
        else:
            for key in ("caption", "title", "extra"):
                val = str(rec.get(key, ""))
                mm = NC_RE.search(val)
                if mm:
                    out.append({"accession": mm.group(0).upper(), "title": title, "strain": strain})
                    break
    # Dedup preserve order
    seen, dedup = set(), []
    for r in out:
        if r["accession"] not in seen:
            seen.add(r["accession"])
            dedup.append(r)
    return dedup

def ncbi_assembly_for_refseq(refseq_acc):
    r = _get(f"{NCBI_BASE}/esearch.fcgi", {
        "db": "assembly",
        "term": refseq_acc,
        "retmode": "json",
        "retmax": 5,
    })
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None
    r2 = _get(f"{NCBI_BASE}/esummary.fcgi", {
        "db": "assembly",
        "id": ",".join(ids),
        "retmode": "json",
    })
    res = r2.json().get("result", {})
    for uid in res.get("uids", []):
        acc = res.get(uid, {}).get("assemblyaccession") or ""
        if GCF_RE.match(acc):
            return acc
    return None

def efetch_fasta_by_accession(accession):
    r = _get(f"{NCBI_BASE}/efetch.fcgi", {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    })
    txt = r.text
    if not txt.strip().startswith(">"):
        raise RuntimeError(f"FASTA not returned for {accession}")
    return txt

# -------- Interactive selection (viral_usher-like) --------
def prompt_taxonomy_id(species_term, tax_entries):
    print(f"\nFound the following matches for '{species_term}':")
    for idx, entry in enumerate(tax_entries):
        print(f"{idx+1}. {entry['sci_name']} (taxid {entry['tax_id']})")
    print(f"{len(tax_entries)+1}. None of the above, go back")
    return prompt_int_choice("Enter a number to choose an option", 1, len(tax_entries)+1)

def prompt_refseq_id(taxid, refseq_entries):
    print(f"\nFound the following RefSeq IDs for Taxonomy ID {taxid}:")
    for idx, entry in enumerate(refseq_entries):
        label = f"{entry['accession']}: {entry.get('title','').strip()}"
        strain = entry.get("strain")
        if strain and strain != "No strain" and strain not in label:
            label += f" ({strain})"
        print(f"{idx+1}. {label}")
    print(f"{len(refseq_entries)+1}. None of the above, go back")
    return prompt_int_choice("Enter a number to choose an option", 1, len(refseq_entries)+1)

def get_species_taxonomy_refseq():
    species_term = ""
    while not species_term:
        species_term = get_input("\nWhat is your virus of interest? ")
    return lookup_taxonomy_refseq(species_term)

def lookup_taxonomy_refseq(species_term):
    print(f"Looking up NCBI Taxonomy entries for '{species_term}'...")
    tax_entries = ncbi_search_taxonomy(species_term) #('"' + species_term + '"') <-- more precise, but less results
    if tax_entries:
        return get_taxonomy_refseq(species_term, tax_entries)
    print(f"\nNo matches found for '{species_term}'.")
    return get_species_taxonomy_refseq()

def get_taxonomy_refseq(species_term, tax_entries):
    choice = prompt_taxonomy_id(species_term, tax_entries)
    if choice > len(tax_entries):
        return get_species_taxonomy_refseq()
    taxid = tax_entries[choice-1]["tax_id"]
    return lookup_refseq(species_term, tax_entries, taxid)

def lookup_refseq(species_term, tax_entries, taxid):
    print(f"Looking up RefSeqs associated with Taxonomy ID {taxid}...")
    refseq_entries = ncbi_refseqs_for_taxid(taxid)
    if refseq_entries:
        return get_refseq(species_term, tax_entries, taxid, refseq_entries)
    print(f"No RefSeqs found for Taxonomy ID '{taxid}'. Try a different match from species search.")
    return get_taxonomy_refseq(species_term, tax_entries)

def get_refseq(species_term, tax_entries, taxid, refseq_entries):
    choice = prompt_refseq_id(taxid, refseq_entries)
    if choice > len(refseq_entries):
        return get_taxonomy_refseq(species_term, tax_entries)
    refseq_id = refseq_entries[choice-1]["accession"]
    print(f"Looking up NCBI Assembly accession for selected RefSeq '{refseq_id}'...")
    assembly_id = ncbi_assembly_for_refseq(refseq_id)
    if not assembly_id:
        print(f"Could not find assembly ID for RefSeq ID {refseq_id} -- can't download RefSeq.")
        print("Try choosing a different RefSeq.")
        return get_refseq(species_term, tax_entries, taxid, refseq_entries)
    return species_term, str(taxid), refseq_id, assembly_id

# -------- Filesystem & Kraken helpers --------
def sanitize_folder_name(name):
    # keep alnum, dot, underscore, dash
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

import glob, gzip

def refseq_seems_in_db(db_dir, refseq_id):
    """
    Detect whether a RefSeq is already present in the Kraken2 DB library.
    Looks for any file in <db>/library/ that starts with '<RefSeq>_'.
    """
    lib = os.path.join(db_dir, "library/added")
    if not os.path.isdir(lib):
        return False

    prefix = f"{refseq_id}_"
    for fname in os.listdir(lib):
        if fname.startswith(prefix):
            return True
    return False

# -------- viral_usher integration (init -> build) --------
def run_viral_usher(species, refseq_id, workdir):
    """
    Use viral_usher init (non-interactive) to create a config, then build.
    """
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

    # Accept either mat.pb.gz (what you want) or optimized.pb.gz (older/alt naming)
    candidates = [
        os.path.join(workdir, "mat.pb.gz"),
        os.path.join(workdir, "optimized.pb.gz"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return None

# -------- CLI --------
def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    ap = argparse.ArgumentParser(
        description="Interactive species→taxid→RefSeq selection; download FASTA; handle MAT via viral_usher; and add to Kraken2 DB."
    )
    ap.add_argument("--db", required=True, help="Path to Kraken2 DB directory.")
    ap.add_argument("--skip-existing-files", action="store_true",
                    help="If <ACCESSION>.fna already exists locally, skip re-download.")
    args = ap.parse_args()

    ensure_folder(ROOT_PATHOGENS_DIR)
    print("=== metaWEPP helper (viral_usher-like init) ===")

    while True:
        # Step 1–3: interactive selection
        species, taxid, refseq_id, assembly_id = get_species_taxonomy_refseq()

        # Ask user what to name the folder inside pathogens_for_wepp (default to species-based)
        default_name = default_species_folder(species)
        user_folder = get_input(f"Folder name inside pathogens_for_wepp [{default_name}]: ").strip() or default_name
        user_folder = sanitize_folder_name(user_folder)
        species_dir = ensure_folder(os.path.join(ROOT_PATHOGENS_DIR, user_folder))
        print(f"[INFO] Working directory for this pathogen: {species_dir}")

        # Download FASTA into the chosen folder
        fasta_path = os.path.join(species_dir, f"{refseq_id}.fna")
        if args.skip_existing_files and os.path.isfile(fasta_path):
            print(f"[SKIP] Reference FASTA already exists: {fasta_path}")
        else:
            print(f"[FETCH] {refseq_id} via NCBI efetch...")
            fasta = efetch_fasta_by_accession(refseq_id)
            fasta_path = save_fasta(species_dir, refseq_id, fasta)
            print(f"[OK] Saved reference FASTA: {fasta_path}")

        # Step 4: presence/DB check -> only prompt/add if not already present in DB library
        if os.path.isfile(fasta_path):
            already_in_db = refseq_seems_in_db(args.db, refseq_id)

            if already_in_db:
                print(f"[INFO] {refseq_id} already present in DB library (prefix match). Skipping add-to-library.")
            else:
                if yes_no("Add this reference FASTA to the Kraken2 DB now?", default_yes=True):
                    add_file_to_kraken(args.db, fasta_path)
                    if refseq_seems_in_db(args.db, refseq_id):
                        print(f"[INFO] Added to Kraken2 DB: {refseq_id} (confirmed in library)")
                    else:
                        print(f"[INFO] Added to Kraken2 DB: {refseq_id}")
                else:
                    print("[INFO] Skipped adding to Kraken2 DB.")
        else:
            print("[WARN] Reference FASTA missing; cannot add to Kraken2 DB.")

        # Step 5: MAT handling
        print("\n(Optional) Provide path to an existing MAT (Mutation-Annotated Tree).")
        print("If omitted, we'll build it with viral_usher. "
              "NOTE: Verify the MAT’s reference matches the downloaded RefSeq.")
        mat_path_in = (get_input("Path to MAT (or press Enter to build): ") or "").strip()

        if mat_path_in:
            if not os.path.isfile(mat_path_in):
                print(f"[WARN] MAT not found: {mat_path_in}")
            else:
                copied = copy_mat(mat_path_in, species_dir)
                print(f"[OK] Copied MAT into {species_dir}: {copied}")
        else:
            print("[INFO] No MAT provided; building via viral_usher (init → build)…")
            try:
                built_mat = run_viral_usher(species, refseq_id, species_dir)
                if built_mat:
                    final_mat = copy_mat(built_mat, species_dir) if os.path.dirname(built_mat) != species_dir else built_mat
                    print(f"[OK] MAT ready: {final_mat}")
                else:
                    print("[WARN] viral_usher completed but MAT file not found. Please check logs and outputs.")
            except FileNotFoundError:
                print("[ERROR] viral_usher executable not found on PATH. Please install/configure viral_usher.")
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] viral_usher failed with exit code {e.returncode}. Check the console output.")

        # Step 6: loop
        if not yes_no("\nDo you want to process another virus?", default_yes=False):
            print("\nAll done. Bye!")
            break

if __name__ == "__main__":
    main()