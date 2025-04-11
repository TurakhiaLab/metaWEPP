import sys
import gzip
from collections import defaultdict

def load_kraken_classification(kraken_file):
    """
    Reads the Kraken classification file and groups reads by Taxon ID.
    Returns a dictionary {taxon_id: set(read_names)}
    """
    taxon_reads = defaultdict(set)

    with open(kraken_file, "r") as file:
        for line in file:
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue  # Skip malformed lines

            read_name = fields[1]  # Read name
            taxon_id = fields[2]   # Taxon ID

            # Store read name under its assigned Taxon ID
            taxon_reads[taxon_id].add(read_name)

    return taxon_reads

def split_fasta_reads(fq_file, taxon_reads):
    """
    Reads the FASTQ file and writes separate FASTQ files for each Taxon ID.
    """
    taxon_files = {taxon_id: open(f"{taxon_id}.fq", "w") for taxon_id in taxon_reads}

    with gzip.open(fq_file, "rt") as file:
        while True:
            # Read four lines per FASTQ entry
            header = file.readline().strip()
            sequence = file.readline().strip()
            plus_line = file.readline().strip()
            quality = file.readline().strip()

            if not header:
                break  # End of file

            read_name = header.split("/")[0][1:]  # Extract read name from FASTQ

            # Check which taxon this read belongs to
            for taxon_id, reads in taxon_reads.items():
                if read_name in reads:
                    taxon_files[taxon_id].write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")

    # Close all output FASTQ files
    for file in taxon_files.values():
        file.close()

def main(kraken_file, fq_file):
    # Load Kraken classified reads
    taxon_reads = load_kraken_classification(kraken_file)


    # Split FASTQ reads based on Taxon ID
    split_fasta_reads(fq_file, taxon_reads)

    print(f"FASTQ reads split into separate files for each Taxon ID.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python split_fq_by_taxon.py <kraken_output_file> <fastq_file>")
        sys.exit(1)

    kraken_file = sys.argv[1]
    fq_file = sys.argv[2]

    main(kraken_file, fq_file)