#!/usr/bin/env python3
"""
Split reads classified by Kraken2 into per-accession FASTQs.

• Reads whose accession is in the reference panel are written to
    <OUT_ROOT>/<accession>/<accession>_{R1,R2}.fq.gz
• Reads whose accession is not in the panel are written to
    <OUT_ROOT>/other_pathogens/<accession>/<accession>_{R1,R2}.fq.gz
"""

import argparse, gzip, json, sys, subprocess, tempfile, shutil, re, multiprocessing, os, random, glob
from multiprocessing import Process
from pathlib import Path
from multiprocessing import Queue

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

def get_genome_length(pathogen_name):
    """
    Calculates the genome length of a pathogen by reading its FASTA file.
    """
    base_path = "data/pathogens_for_wepp"
    pathogen_path = os.path.join(base_path, pathogen_name)
    
    # Find FASTA files (.fna or .fasta or .fa)
    fasta_files = glob.glob(os.path.join(pathogen_path, "*.fna")) + glob.glob(os.path.join(pathogen_path, "*.fasta")) + glob.glob(os.path.join(pathogen_path, "*.fa"))

    total_length = 0
    for f_path in fasta_files:
        with open(f_path, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    total_length += len(line.strip())
    return total_length

def split_fastq(fq_path, mate, read2species, out_dir, threads, queue=None):
    """
    Splits a FASTQ file based on read classification.

    For each read, it looks up its assigned species name from `read2species`.
    - If the species is a WEPP pathogen, reads are written to
      `out_dir/<species_name>/<species_name>_{mate}.fastq.gz`.
    - Otherwise, reads are written to
      `out_dir/Other_Pathogens/<species_name>/<species_name>_{mate}.fastq.gz`.

    This function uses pigz for parallel gzip compression to improve speed.
    It returns a dictionary with the total read length for each WEPP pathogen.
    """
    # Get the set of WEPP pathogen names
    taxons_wepp_pathogens = get_taxid_of_pathogens_for_wepp()
    wepp_pathogen_names = set(taxons_wepp_pathogens.values())

    # Dictionary to hold pigz processes and their stdin pipes
    writers = {}
    # Dictionary to cache the output directory for each species
    species_to_dir = {}
    # Dictionary to hold the process objects
    procs = {}
    # Dictionary to store read lengths for coverage calculation
    pathogen_read_lengths = {}

    try:
        # Use pigz for parallel decompression, which is faster than gzip module
        pigz_proc = subprocess.Popen(["pigz", "-dc", fq_path], stdout=subprocess.PIPE, text=True)
        f_in = pigz_proc.stdout

        while True:
            header = f_in.readline()
            if not header:
                break
            seq = f_in.readline()
            plus = f_in.readline()
            qual = f_in.readline()

            # Extract read ID from the header, e.g., '@A00... 1:N:0:1' -> 'A00...'
            read_id = header.split()[0][1:]
            if read_id.endswith('/1') or read_id.endswith('/2'):
                read_id = read_id[:-2]

            if read_id in read2species:
                species_names = read2species[read_id]

                # Reads that map to a single species
                if len(species_names) == 1:
                    species_name = species_names[0]
                    # Removing blank spaces and other such characters in species_name for the filename
                    sanitized_name = re.sub(r'[^a-zA-Z0-9_-]', '_', species_name)

                    # Update read lengths for coverage calculation
                    if species_name in wepp_pathogen_names:
                        pathogen_read_lengths.setdefault(species_name, 0)
                        pathogen_read_lengths[species_name] += len(seq.strip())

                    # Determine output directory for the species if not already seen
                    if species_name not in species_to_dir:
                        if species_name in wepp_pathogen_names:
                            target_dir = Path(out_dir) / sanitized_name
                        else:
                            target_dir = Path(out_dir) / "Other_Pathogens" / sanitized_name
                        target_dir.mkdir(parents=True, exist_ok=True)
                        species_to_dir[species_name] = target_dir

                    # Get or create a writer for the species
                    if species_name not in writers:
                        target_dir = species_to_dir[species_name]
                        out_path = target_dir / f"{sanitized_name}_{mate}.fastq.gz"
                        
                        # Use pigz for faster compression
                        proc = subprocess.Popen(
                            ["pigz", "-p", str(threads)],
                            stdin=subprocess.PIPE,
                            stdout=open(out_path, "wb")
                        )
                        writers[species_name] = proc.stdin
                        procs[species_name] = proc

                    # Write the FASTQ record
                    writer = writers[species_name]
                    record = f"{header}{seq}{plus}{qual}"
                    writer.write(record.encode('utf-8'))
                
                # Reads that map to multiple species -> Use minimap2 to resolve ambiguity
                else:
                    pass
                    
    finally:
        # Ensure the decompression process is terminated
        if 'pigz_proc' in locals() and pigz_proc.poll() is None:
            pigz_proc.terminate()
            pigz_proc.wait()

        # Close all file writers and wait for pigz compression processes to finish
        for writer in writers.values():
            if not writer.closed:
                writer.close()
        for proc in procs.values():
            proc.wait()

    if queue:
        queue.put(pathogen_read_lengths)
    else:
        return pathogen_read_lengths

def get_read_classified_species(kraken_out, kraken_report):
    taxons_wepp_pathogens = get_taxid_of_pathogens_for_wepp()
    
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
                if sp_taxid in taxons_wepp_pathogens:
                    wepp_species.append(taxons_wepp_pathogens[sp_taxid])

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
                                wepp_names = set(taxons_wepp_pathogens.values())
                                wepp_intersection = names & wepp_names
                                if wepp_intersection:
                                    read_to_name[rid] = list(wepp_intersection)
                                else:
                                    # If no intersection, pick a random name
                                    read_to_name[rid] = [random.choice(tuple(names))]
                            else:
                                read_to_name[rid] = list(names)

    return read_to_name
    
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

    if not shutil.which("pigz"):
        sys.exit("Error: 'pigz' is not installed or not in PATH. Please install it to continue.")
    
    Path(a.out_dir).mkdir(parents=True, exist_ok=True)
    # Per-read classified taxid from Kraken output
    read_to_name = get_read_classified_species(a.kraken_out, a.kraken_report)
    if not read_to_name:
        sys.exit("No classified reads in Kraken output.")

    r2_process = None
    queue = None
    if a.r2:
        queue = Queue()
        r2_process = Process(
            target=split_fastq,
            args=(
                a.r2, "R2", read_to_name,
                a.out_dir, a.threads, queue
            )
        )
        r2_process.start()

    # Always process R1 in this process
    r1_lengths = split_fastq(
        a.r1, "R1", read_to_name,
        a.out_dir, a.threads
    )

    # Wait for R2 if running and combine results
    total_pathogen_read_lengths = r1_lengths
    if r2_process:
        r2_lengths = queue.get()
        r2_process.join()
        for pathogen, length in r2_lengths.items():
            total_pathogen_read_lengths[pathogen] = total_pathogen_read_lengths.get(pathogen, 0) + length

    # Calculate and write coverage for each pathogen_for_wepp
    taxons_wepp_pathogens = get_taxid_of_pathogens_for_wepp()
    wepp_pathogen_names = set(taxons_wepp_pathogens.values())
    coverage_file_path = Path(a.out_dir) / "pathogen_coverage.tsv"
    with open(coverage_file_path, "w") as f_cov:
        for pathogen_name in wepp_pathogen_names:
            genome_length = get_genome_length(pathogen_name)
            if genome_length > 0:
                total_read_length = total_pathogen_read_lengths.get(pathogen_name, 0)
                coverage = total_read_length / genome_length
                f_cov.write(f"{pathogen_name}\t{coverage}\n")

if __name__ == "__main__":
    main()