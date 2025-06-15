import os
import sys
import csv
import subprocess
import hashlib

# Start dummy taxid at a high reserved range to avoid real taxids
next_dummy_taxid = 9000000
assigned_dummy_taxids = {}

def get_fasta_length(fasta_path):
    total_length = 0
    try:
        with open(fasta_path, 'r') as fasta_file:
            for line in fasta_file:
                if not line.startswith(">"):
                    total_length += len(line.strip())
    except FileNotFoundError:
        print(f"-------- FASTA file does not exist: {fasta_path} --------")
        return 0
    return total_length

def get_taxid_from_name(name):
    """Try taxonkit lookup; fallback to stable pseudo-taxid."""
    try:
        result = subprocess.run(
            ["taxonkit", "name2taxid"],
            input=name + "\n",
            text=True,
            capture_output=True,
            check=True
        )
        taxid_line = result.stdout.strip()
        if not taxid_line or "\t" not in taxid_line:
            raise ValueError("No valid taxid returned")
        taxid = taxid_line.split("\t")[1]
        if not taxid.isdigit():
            raise ValueError("Non-numeric taxid")
        return taxid
    except Exception:
        # Fallback: generate deterministic pseudo-taxid using hash
        hash_digest = hashlib.md5(name.encode()).hexdigest()
        dummy_taxid = int(hash_digest[:6], 16) + 9000000  # stays within dummy range
        print(f"-------- TaxID not found for {name}. Assigned dummy taxid: {dummy_taxid} --------")
        return str(dummy_taxid)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python generate_mess_tsv.py <output_tsv> <fasta_file> <coverage>")
        sys.exit(1)

    output_tsv = sys.argv[1]
    fasta_path = sys.argv[2]
    coverage = sys.argv[3]

    os.makedirs(os.path.dirname(output_tsv), exist_ok=True)

    fasta_name = os.path.basename(fasta_path).replace(".fa", "").replace(".fasta", "")
    total_length = get_fasta_length(fasta_path)
    tax_id = get_taxid_from_name(fasta_name)

    with open(output_tsv, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerow(["fasta", "path", "cov_sim", "tax_id", "total_sequence_length"])
        writer.writerow([fasta_name, fasta_path, coverage, tax_id, total_length])
        print(f"-------- Wrote TSV for {fasta_name} → {output_tsv} --------")
