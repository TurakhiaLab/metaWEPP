#!/usr/bin/env python3
"""
Evaluate Kraken classifications using FASTQ files and a genome-to-taxid mapping.
Supports either a user-supplied JSON map or Kraken's seqid2taxid.map.
"""

import argparse
import json
import os
import sys
from collections import defaultdict


def parse_fastq_genomes(fq_file):
    """
    Parse a FASTQ file and return a mapping from read ID to genome accession.
    Assumes read ID format: <accession>-<suffix>/1 or /2
    Extracts accession by removing suffix after the last '-'.
    """
    read_to_genome = {}
    try:
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
                    full_id = header[1:].split()[0]
                    core_id = full_id.rsplit('/', 1)[0]
                    genome_accession = core_id.rsplit('-', 1)[0]  # split on last dash
                    read_to_genome[core_id] = genome_accession
    except Exception as e:
        print(f"[ERROR] Failed to parse FASTQ file: {fq_file} — {str(e)}", file=sys.stderr)
        sys.exit(1)

    if not read_to_genome:
        print(f"[ERROR] No reads found in FASTQ file: {fq_file}", file=sys.stderr)
        sys.exit(1)

    return read_to_genome


def load_seqid2taxid_map(kraken_db_path):
    """
    Load seqid2taxid.map from a Kraken database and return a dict.
    """
    seqid2taxid_path = os.path.join(kraken_db_path, "seqid2taxid.map")
    if not os.path.isfile(seqid2taxid_path):
        print(f"[ERROR] seqid2taxid.map not found in {kraken_db_path}", file=sys.stderr)
        sys.exit(1)

    mapping = {}
    try:
        with open(seqid2taxid_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 2:
                    seqid, taxid = parts
                    mapping[seqid] = taxid
    except Exception as e:
        print(f"[ERROR] Failed to read seqid2taxid.map — {str(e)}", file=sys.stderr)
        sys.exit(1)

    return mapping


def evaluate_kraken(kraken_file, fq1, fq2, genome_to_taxid):
    read_to_genome = parse_fastq_genomes(fq1)
    read_to_genome.update(parse_fastq_genomes(fq2))

    total = correct = false_positives = false_negatives = 0
    skipped_reads = unmatched_reads = 0
    per_genome_stats = defaultdict(lambda: {"total": 0, "correct": 0, "false_positive": 0, "false_negative": 0})

    try:
        with open(kraken_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    skipped_reads += 1
                    continue

                status, read_id, predicted_taxid_str = parts[0], parts[1], parts[2]
                try:
                    predicted_taxid = int(predicted_taxid_str)
                except ValueError:
                    skipped_reads += 1
                    continue

                genome_name = read_to_genome.get(read_id)
                if not genome_name:
                    unmatched_reads += 1
                    print(f"[WARNING] Read ID '{read_id}' not found in FASTQ mapping.", file=sys.stderr)
                    continue

                true_taxid_str = genome_to_taxid.get(genome_name)
                if not true_taxid_str:
                    print(f"[ERROR] Genome ID '{genome_name}' not found in taxid mapping.", file=sys.stderr)
                    sys.exit(1)

                try:
                    true_taxid = int(true_taxid_str)
                except ValueError:
                    print(f"[ERROR] Invalid taxid '{true_taxid_str}' for genome '{genome_name}'", file=sys.stderr)
                    sys.exit(1)

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
    except FileNotFoundError:
        print(f"[ERROR] Kraken file '{kraken_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Failed to process Kraken file: {kraken_file} — {str(e)}", file=sys.stderr)
        sys.exit(1)

    if total == 0:
        print("[ERROR] No valid Kraken classifications were evaluated. Exiting.", file=sys.stderr)
        sys.exit(1)

    precision = correct / (correct + false_positives) if (correct + false_positives) > 0 else 0
    recall = correct / (correct + false_negatives) if (correct + false_negatives) > 0 else 0
    accuracy = correct / total if total > 0 else 0

    print("\nOverall Evaluation Summary")
    print("--------------------------")
    print(f"Total evaluated reads: {total}")
    print(f"Correct classifications: {correct}")
    print(f"False positives: {false_positives}")
    print(f"False negatives (unclassified): {false_negatives}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"Accuracy:  {accuracy:.4f}")
    print(f"Skipped malformed Kraken lines: {skipped_reads}")
    print(f"Unmatched read IDs: {unmatched_reads}")

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


def main():
    parser = argparse.ArgumentParser(description="Evaluate Kraken classifications using FASTQ and taxid mapping.")
    parser.add_argument('--fq1', required=True, help="FASTQ file R1")
    parser.add_argument('--fq2', required=True, help="FASTQ file R2")
    parser.add_argument('--kraken', required=True, help="Kraken output file")
    parser.add_argument('--taxid_map', help="JSON file mapping genome IDs to taxids")
    parser.add_argument('--kraken_db', help="Path to Kraken database (uses seqid2taxid.map)")

    args = parser.parse_args()

    for file_path in [args.fq1, args.fq2, args.kraken]:
        if not os.path.isfile(file_path):
            print(f"[ERROR] File not found: {file_path}", file=sys.stderr)
            sys.exit(1)

    if args.taxid_map:
        try:
            with open(args.taxid_map, 'r') as f:
                genome_to_taxid = json.load(f)
        except Exception as e:
            print(f"[ERROR] Failed to load JSON taxid map — {str(e)}", file=sys.stderr)
            sys.exit(1)
    elif args.kraken_db:
        genome_to_taxid = load_seqid2taxid_map(args.kraken_db)
    else:
        print("[ERROR] Either --taxid_map or --kraken_db must be provided.", file=sys.stderr)
        sys.exit(1)

    evaluate_kraken(args.kraken, args.fq1, args.fq2, genome_to_taxid)


if __name__ == "__main__":
    main()