#!/usr/bin/env python3
"""
Improved Evaluation Script for Kraken Classifications

This script computes overall and per-genome evaluation metrics using
the ground truth FASTQ files and Kraken output. It prints a summary
of the overall performance and also a breakdown by genome.
"""

from collections import defaultdict

# File paths (hardcoded)
GROUND_TRUTH_FILE_1 = "/home/jseangmany@AD.UCSD.EDU/art/snakemake_art/output/combined1_R1.fq"
GROUND_TRUTH_FILE_2 = "/home/jseangmany@AD.UCSD.EDU/art/snakemake_art/output/combined1_R2.fq"
KRAKEN_OUTPUT_FILE = "/home/jseangmany@AD.UCSD.EDU/art/snakemake_art/kraken_reports/kraken_output2.txt"

# Ground truth mapping from genome accession to taxonomic ID
GROUND_TRUTH_TAXIDS = {
    "NC_000913.3|kraken:taxid|511145": "511145",   # E. coli K-12 substr. MG1655
    "AF013254.1": "11250",     # Some genome
    "ON811098.1": "2697049",   # Some genome
    "OQ557947.1": "10244"      # Some genome
}


def parse_fastq_genomes(fq_file):
    """
    Parse a FASTQ file and return a dictionary mapping read ID to genome accession.
    
    Expected header format: '@<genome_accession>-<other_info>/[1|2]'
    Example: '@AF013254.1-8000/1' yields a mapping: key = 'AF013254.1-8000', value = 'AF013254.1'
    """
    read_to_genome = {}
    with open(fq_file, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            if not qual:
                break
            if header.startswith('@'):
                # Remove initial '@'
                full_id = header[1:]
                # Remove read pair indicator by splitting at '/'
                core_id = full_id.split('/')[0]
                # Extract genome accession (assumes accession comes before any '-' delimiter)
                genome_accession = core_id.split('-')[0]
                read_to_genome[core_id] = genome_accession
    return read_to_genome


def evaluate_kraken(kraken_file, fq1, fq2, genome_to_taxid):
    """
    Evaluate Kraken predictions using ground truth from FASTQ files.
    
    Computes:
      - Overall metrics (precision, recall, accuracy)
      - Per-genome metrics
    """
    # Build a combined mapping of read ID to genome accession from both FASTQ files.
    read_to_genome = parse_fastq_genomes(fq1)
    read_to_genome.update(parse_fastq_genomes(fq2))
    
    # Overall counters
    total = 0
    correct = 0
    false_positives = 0
    false_negatives = 0
    
    # Per-genome statistics: each key maps to a dict of counters.
    per_genome_stats = defaultdict(lambda: {"total": 0, "correct": 0, "false_positive": 0, "false_negative": 0})
    
    with open(kraken_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue  # Skip malformed lines
            
            status, read_id, predicted_taxid_str = parts[0], parts[1], parts[2]
            try:
                predicted_taxid = int(predicted_taxid_str)
            except ValueError:
                continue  # Skip if predicted taxid isn't an integer
            
            # Get the genome accession based on the read ID.
            genome_name = read_to_genome.get(read_id)
            if not genome_name:
                continue  # Skip if read not found in ground truth
            
            true_taxid_str = genome_to_taxid.get(genome_name)
            if not true_taxid_str:
                continue  # Skip if genome not in our mapping
            
            try:
                true_taxid = int(true_taxid_str)
            except ValueError:
                continue  # Skip if true taxid isn't valid
            
            # Update counters.
            total += 1
            per_genome_stats[genome_name]["total"] += 1
            
            if predicted_taxid == 0:
                false_negatives += 1
                per_genome_stats[genome_name]["false_negative"] += 1
            elif predicted_taxid == true_taxid:
                correct += 1
                per_genome_stats[genome_name]["correct"] += 1
            else:
                false_positives += 1
                per_genome_stats[genome_name]["false_positive"] += 1
    
    precision = correct / (correct + false_positives) if (correct + false_positives) > 0 else 0
    recall = correct / (correct + false_negatives) if (correct + false_negatives) > 0 else 0
    accuracy = correct / total if total > 0 else 0
    
    # Print overall evaluation.
    print("\nOverall Evaluation Summary")
    print("--------------------------")
    print(f"Total evaluated reads: {total}")
    print(f"Correct classifications: {correct}")
    print(f"False positives: {false_positives}")
    print(f"False negatives (unclassified): {false_negatives}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"Accuracy:  {accuracy:.4f}")
    
    # Print per-genome evaluation.
    print("\nPer-Genome Evaluation Summary")
    print("-----------------------------")
    for genome, stats in per_genome_stats.items():
        total_genome = stats["total"]
        correct_genome = stats["correct"]
        fp_genome = stats["false_positive"]
        fn_genome = stats["false_negative"]
        precision_genome = correct_genome / (correct_genome + fp_genome) if (correct_genome + fp_genome) > 0 else 0
        recall_genome = correct_genome / (correct_genome + fn_genome) if (correct_genome + fn_genome) > 0 else 0
        accuracy_genome = correct_genome / total_genome if total_genome > 0 else 0
        
        print(f"Genome: {genome}")
        print(f"  Total evaluated reads: {total_genome}")
        print(f"  Correct classifications: {correct_genome}")
        print(f"  False positives: {fp_genome}")
        print(f"  False negatives (unclassified): {fn_genome}")
        print(f"  Precision: {precision_genome:.4f}")
        print(f"  Recall:    {recall_genome:.4f}")
        print(f"  Accuracy:  {accuracy_genome:.4f}\n")


if __name__ == "__main__":
    evaluate_kraken(KRAKEN_OUTPUT_FILE, GROUND_TRUTH_FILE_1, GROUND_TRUTH_FILE_2, GROUND_TRUTH_TAXIDS)
