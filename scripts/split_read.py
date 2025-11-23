#!/usr/bin/env python3
"""
Split reads classified by Kraken2 into per-accession FASTQs.

• Reads whose accession is in the reference panel are written to
    <OUT_ROOT>/<accession>/<accession>_{R1,R2}.fq.gz
• Reads whose accession is not in the panel are written to
    <OUT_ROOT>/other_pathogens/<accession>/<accession>_{R1,R2}.fq.gz
"""

import argparse, gzip, json, sys, subprocess, tempfile, shutil, re, concurrent, os, random, glob
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

def get_taxid_of_pathogens_for_wepp():
    added_taxons_path = "data/pathogens_for_wepp/added_taxons.csv"
    added_taxons = {}
    if os.path.exists(added_taxons_path):
        with open(added_taxons_path, "r") as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split(',')
                    if len(parts) == 2:
                        taxid, folder = parts
                        added_taxons.setdefault(taxid, set()).add(folder)
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

def get_wepp_pathogen_fasta_files(species_name):
    """
    Finds FASTA files for a given WEPP pathogen species.
    """
    base_path = "data/pathogens_for_wepp"
    pathogen_path = os.path.join(base_path, species_name)
    
    # Find FASTA files (.fna or .fasta or .fa)
    fasta_files = glob.glob(os.path.join(pathogen_path, "*.fna")) + glob.glob(os.path.join(pathogen_path, "*.fasta")) + glob.glob(os.path.join(pathogen_path, "*.fa"))
    return fasta_files

def split_fastq(fq_path, mate, read2species, out_dir, threads):
    """
    Splits a FASTQ file based on read classification.
    This is a single-threaded implementation.
    """
    taxons_wepp_pathogens = get_taxid_of_pathogens_for_wepp()
    wepp_pathogen_names = set()
    for s in taxons_wepp_pathogens.values():
        wepp_pathogen_names.update(s)

    pathogen_read_lengths = {}
    ambiguous_reads_info = {}
    
    writers = {}
    ambiguous_writers = {}

    try:
        pigz_proc = subprocess.Popen(["pigz", "-dc", fq_path], stdout=subprocess.PIPE, text=True)
        f_in = pigz_proc.stdout

        while True:
            header = f_in.readline()
            if not header:
                break
            seq = f_in.readline()
            plus = f_in.readline()
            qual = f_in.readline()

            read_id = header.split()[0][1:]
            # Handle cases where read IDs don't match the read names from kraken output
            if read_id.endswith('/1') or read_id.endswith('/2'):
                read_id = read_id[:-2]

            # Skip unclassified reads
            if read_id not in read2species:
                continue
            
            species_names = read2species[read_id]
            record = f"{header}{seq}{plus}{qual}"

            # Unique mapping -> Assign to that species
            if len(species_names) == 1:
                species_name = species_names[0]
                is_wepp_pathogen = species_name in wepp_pathogen_names

                # Update pathogen_read_lengths
                if is_wepp_pathogen:
                    pathogen_read_lengths.setdefault(species_name, 0)
                    pathogen_read_lengths[species_name] += len(seq.strip())

                if species_name not in writers:
                    if is_wepp_pathogen:
                        target_dir = Path(out_dir) / species_name
                    else:
                        target_dir = Path(out_dir) / "Other_Pathogens" / species_name
                    target_dir.mkdir(parents=True, exist_ok=True)
                    out_path = target_dir / f"{species_name}_{mate}.fastq.gz"

                    out_f = open(out_path, "ab")
                    proc = subprocess.Popen(
                        ["pigz", "-p", str(threads)],
                        stdin=subprocess.PIPE,
                        stdout=out_f
                    )
                    writers[species_name] = (proc, out_f)

                proc, out_f = writers[species_name]
                proc.stdin.write(record.encode("utf-8"))

            # Ambiguous mapping -> Store in temp files for later resolution
            else:  
                group_name = "_".join(sorted(species_names))
                
                if group_name not in ambiguous_reads_info:
                    fasta_files = []
                    for name in species_names:
                        fasta_files.extend(get_wepp_pathogen_fasta_files(name))
                    
                    temp_dir = Path(out_dir) / "temp_ambiguous" / group_name
                    temp_dir.mkdir(parents=True, exist_ok=True)
                    
                    ambiguous_reads_info[group_name] = {
                        "fastas": fasta_files,
                        "dir": temp_dir
                    }

                if group_name not in ambiguous_writers:
                    out_path = ambiguous_reads_info[group_name]["dir"] / f"ambiguous_{mate}.fastq.gz"
                    out_f = open(out_path, "ab")
                    proc = subprocess.Popen(
                        ["pigz", "-p", str(threads)],
                        stdin=subprocess.PIPE,
                        stdout=out_f
                    )
                    ambiguous_writers[group_name] = (proc, out_f)
                
                proc, out_f = ambiguous_writers[group_name]
                proc.stdin.write(record.encode("utf-8"))

    finally:
        for proc, out_f in writers.values():
            proc.stdin.close()      
            proc.wait()
            out_f.close()

        for proc, out_f in ambiguous_writers.values():
            proc.stdin.close()
            proc.wait()
            out_f.close()

        if 'pigz_proc' in locals():
            pigz_proc.terminate()
            pigz_proc.wait()

    return pathogen_read_lengths, ambiguous_reads_info

def resolve_ambiguous_reads(ambiguous_reads_info, out_dir, seq_type, threads):
    # Get the set of WEPP pathogen names
    taxons_wepp_pathogens = get_taxid_of_pathogens_for_wepp()
    wepp_pathogen_names = set()
    for s in taxons_wepp_pathogens.values():
        wepp_pathogen_names.update(s)
    
    writers = {}
    procs = {}
    resolved_pathogen_read_lengths = {}
    
    try:
        for group_name, info in ambiguous_reads_info.items():
            temp_dir = info['dir']
            
            # Check for R1 and R2 reads
            r1_path = temp_dir / "ambiguous_R1.fastq.gz"
            r2_path = temp_dir / "ambiguous_R2.fastq.gz"
            
            # Concatenate all reference fastas for this group into a single file
            ref_fasta = temp_dir / "ref.fasta"
            with open(ref_fasta, "wb") as f_out:
                for fasta_path in sorted(set(info["fastas"])):
                    with open(fasta_path, "rb") as f_in:
                        shutil.copyfileobj(f_in, f_out)
            
            sam_path = temp_dir / "alignments.sam"
            
            # Decide the preset
            if seq_type in ("d", "s"):
                preset = "sr"          # short reads
            elif seq_type == "n":
                preset = "map-ont"     # nanopore
            else:
                raise ValueError(f"Invalid seq_type '{seq_type}'. Expected 'd', 's', or 'n'.")
            
            # Build minimap2 command
            minimap_cmd = [
                "minimap2",
                "-a", # Output SAM format
                "--secondary=no", # Do not output secondary alignments
                "-x", preset,
                "-t", str(threads),
                str(ref_fasta),
            ]
            
            reads_to_process = []
            if r1_path.exists() and r2_path.exists():
                minimap_cmd.extend([str(r1_path), str(r2_path)])
                reads_to_process.append(("R1", r1_path))
                reads_to_process.append(("R2", r2_path))
            elif r1_path.exists():
                minimap_cmd.extend([str(r1_path)])
                reads_to_process.append(("R1", r1_path))
            else:
                raise ValueError(f"Fastq file not found in '{temp_dir}'.")

            with open(sam_path, "w") as f:
                subprocess.run(minimap_cmd, stdout=f, check=True)

            # Build a map from reference sequence ID to species name
            ref_to_species = {}
            for fasta_path in info["fastas"]:
                path_parts = Path(fasta_path).parts
                species_name = path_parts[-2]
                
                _open = gzip.open if fasta_path.endswith(".gz") else open
                with _open(fasta_path, "rt") as f:
                    for line in f:
                        if line.startswith(">"):
                            ref_id = line.split()[0][1:]
                            ref_to_species[ref_id] = species_name

            # Parse SAM and assign reads to species
            read_to_best_species = {}
            with open(sam_path, "r") as f:
                for line in f:
                    if line.startswith('@'):
                        continue
                    parts = line.strip().split("\t")
                    read_id = parts[0]
                    ref_name = parts[2]
                    # unmapped read
                    if ref_name == '*' or ref_name not in ref_to_species:
                        continue
                    
                    species_name = ref_to_species[ref_name]

                    align_score = 0
                    # Look for AS:i: tag for alignment score
                    for i in range(11, len(parts)):
                        if parts[i].startswith("AS:i:"):
                            align_score = int(parts[i][5:])
                            break
                    
                    # Paired reads may have the same read_id (read name), so picking the reference with best alignment score
                    if read_id not in read_to_best_species or align_score > read_to_best_species[read_id][1]:
                        read_to_best_species[read_id] = (species_name, align_score)

            for mate, fq_path in reads_to_process:
                with gzip.open(fq_path, "rt") as f_in:
                    while True:
                        header = f_in.readline()
                        if not header: break
                        seq = f_in.readline()
                        plus = f_in.readline()
                        qual = f_in.readline()
                        
                        read_id = header.split()[0][1:]
                        # Handle cases where read IDs don't match the read names from SAM output
                        if read_id.endswith('/1') or read_id.endswith('/2'):
                            read_id = read_id[:-2]
                        if read_id in read_to_best_species:
                            species_name, _ = read_to_best_species[read_id]
                            resolved_pathogen_read_lengths.setdefault(species_name, 0)
                            resolved_pathogen_read_lengths[species_name] += len(seq.strip())
                            
                            writer_key = (species_name, mate)
                            if writer_key not in writers:
                                target_dir = Path(out_dir) / species_name
                                target_dir.mkdir(parents=True, exist_ok=True)
                                out_path = target_dir / f"{species_name}_{mate}.fastq.gz"
                                
                                # Use pigz for faster compression
                                out_f = open(out_path, "ab")  # Append mode
                                proc = subprocess.Popen(
                                    ["pigz", "-p", str(threads)],
                                    stdin=subprocess.PIPE,
                                    stdout=out_f
                                )
                                writers[writer_key] = proc.stdin
                                procs[writer_key] = (proc, out_f)
                            
                            writer = writers[writer_key]
                            record = f"{header}{seq}{plus}{qual}"
                            writer.write(record.encode('utf-8'))
                            writer.flush() 
            
        # Clean up the parent temp_ambiguous directory
        temp_ambiguous_dir = Path(out_dir) / "temp_ambiguous"
        if temp_ambiguous_dir.exists():
            shutil.rmtree(temp_ambiguous_dir)

    finally:
        for (proc, out_f) in procs.values():
            proc.stdin.close()
            proc.wait()
            out_f.close()

    return resolved_pathogen_read_lengths

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
        # Sanitize species name from kraken report
        name = re.sub(r'[^a-zA-Z0-9_-]', '_', name)
        rank = parts[3].strip()

        # In case of genus, if the member species are in wepp_pathogens, then we need to map its taxon with those of wepp species
        if rank.startswith('G'):
            species_under_genus = []
            j = i + 1
            while j < len(lines):
                next_parts = lines[j].split("\t")
                if len(next_parts) >= 6 and next_parts[3].strip().startswith('S'):
                    # Sanitize species name from kraken report
                    sp_name = re.sub(r'[^a-zA-Z0-9_-]', '_', next_parts[5].strip())
                    species_under_genus.append((next_parts[4].strip(), sp_name))
                    j += 1
                else:
                    break
            
            wepp_species = []
            for sp_taxid, _ in species_under_genus:
                if sp_taxid in taxons_wepp_pathogens:
                    wepp_species.append(taxons_wepp_pathogens[sp_taxid])

            # if there are wepp species under this genus then all the reads of genus and the species_under_genus should map to all the wepp species
            if wepp_species: 
                wepp_species_names = set()
                for item in wepp_species:
                    if isinstance(item, set):
                        wepp_species_names.update(item)
                    else:
                        wepp_species_names.add(item)
                taxid_to_names.setdefault(taxid, set()).update(wepp_species_names)
                for sp_taxid, _ in species_under_genus:
                    taxid_to_names.setdefault(sp_taxid, set()).update(wepp_species_names)
            # if no wepp species under this genus, just map the genus to the single species when (species_under_genus) == 1
            else:   
                if len(species_under_genus) == 1:
                    taxid_to_names.setdefault(taxid, set()).add(species_under_genus[0][1])
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

    # 2. Get read→tax_id from kraken_output
    read_to_name = {}
    with open(kraken_out, "r") as f:
            for line in f:
                if line.startswith("C\t"):
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 3:
                        _, rid, taxid = parts[:3]

                        # If the taxid is for a WEPP pathogen, use the folder name directly.
                        if taxid in taxons_wepp_pathogens:
                            read_to_name[rid] = list(taxons_wepp_pathogens[taxid])
                            continue

                        # Otherwise, use the scientific name from the kraken report.
                        if taxid in taxid_to_names:
                            names = taxid_to_names[taxid]
                            if len(names) > 1:
                                # In case of wepp species, prefer all those names
                                wepp_names = set()
                                for s in taxons_wepp_pathogens.values():
                                    wepp_names.update(s)
                                wepp_intersection = names & wepp_names
                                if wepp_intersection:
                                    read_to_name[rid] = list(wepp_intersection)
                                else:
                                    # If no wepp_pathogen_names, pick a random name from names
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
    ap.add_argument("-s","--sequencing-type", required=True)

    return ap.parse_args()

def main():
    a = parse_args()

    if not shutil.which("pigz"):
        sys.exit("Error: 'pigz' is not installed or not in PATH. Please install it to continue.")
    
    if not shutil.which("minimap2"):
        sys.exit("Error: 'minimap2' is not installed or not in PATH. Please install it to continue.")

    Path(a.out_dir).mkdir(parents=True, exist_ok=True)
    # Per-read classified taxid from Kraken output
    read_to_name = get_read_classified_species(a.kraken_out, a.kraken_report)
    if not read_to_name:
        sys.exit("No classified reads in Kraken output.")

    results = {}
    with ProcessPoolExecutor(max_workers=2) as exe:
        futures = {}

        # Launch R1
        futures['R1'] = exe.submit(
            split_fastq, a.r1, "R1", read_to_name, a.out_dir, a.threads
        )

        # Launch R2 (if present)
        if a.r2:
            futures['R2'] = exe.submit(
                split_fastq, a.r2, "R2", read_to_name, a.out_dir, a.threads
            )

        # Collect results
        for mate, fut in futures.items():
            results[mate] = fut.result()

    # Combine R1 + R2 results
    total_pathogen_read_lengths, total_ambiguous_reads = results['R1']

    if 'R2' in results:
        r2_lengths, r2_ambiguous = results['R2']

        for pathogen, length in r2_lengths.items():
            total_pathogen_read_lengths[pathogen] = \
                total_pathogen_read_lengths.get(pathogen, 0) + length

        for group, info in r2_ambiguous.items():
            if group not in total_ambiguous_reads:
                total_ambiguous_reads[group] = info
    
    if total_ambiguous_reads:
        resolved_lengths = resolve_ambiguous_reads(total_ambiguous_reads, a.out_dir, a.sequencing_type, a.threads)
        for pathogen, length in resolved_lengths.items():
            total_pathogen_read_lengths[pathogen] = total_pathogen_read_lengths.get(pathogen, 0) + length

    # Calculate and write coverage for each pathogen_for_wepp
    taxons_wepp_pathogens = get_taxid_of_pathogens_for_wepp()
    wepp_pathogen_names = set()
    for s in taxons_wepp_pathogens.values():
        wepp_pathogen_names.update(s)
        
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