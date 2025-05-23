import os
import sys
import csv

def get_fasta_length(fasta_path):
    total_length = 0
    try:
        with open(fasta_path, 'r') as fasta_file:
            for line in fasta_file:
                if not line.startswith(">"):
                    total_length += len(line.strip())
    except FileNotFoundError:
        print(f"-------- FASTA file does not exist: {fasta_path} --------")
        sys.exit(1)
    return total_length

def update_tsv(tsv_path, fasta_path, coverage, tax_id):
    """Creates or updates the TSV file for MeSS. Skips if entry already exists."""
    fasta_name = os.path.basename(fasta_path).replace(".fa", "").replace(".fasta", "")
    total_length = get_fasta_length(fasta_path)

    new_entry = [fasta_name, fasta_path, coverage, tax_id, total_length]

    file_exists = os.path.exists(tsv_path)
    entry_exists = False

    # Check if entry exist
    if file_exists:
        with open(tsv_path, 'r') as tsv_file:
            reader = csv.reader(tsv_file, delimiter='\t')
            for row in reader:
                if len(row) >= 2 and row[1] == fasta_path:
                    entry_exists = True
                    break

    # If entry exists, skip
    if entry_exists:
        print(f"-------- Entry for {fasta_path} already exists. Skip adding. --------")
        return

    # Write the new entry
    with open(tsv_path, 'a', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        # If file does not exist, write the header first
        if not file_exists:
            writer.writerow(["fasta", "path", "cov_sim", "tax_id", "total_sequence_length"])
        writer.writerow(new_entry)

    print(f"-------- Added new entry for: {fasta_path} --------")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python generate_mess_tsv.py <tsv_path> <fasta_path> <coverage> <tax_id>")
        sys.exit(1)

    tsv_path = sys.argv[1]
    fasta_path = sys.argv[2]
    coverage = sys.argv[3]
    tax_id = sys.argv[4]

    update_tsv(tsv_path, fasta_path, coverage, tax_id)

