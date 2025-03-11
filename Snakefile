# Snakefile

rule all:
    input:
        "split_fq.done"

###############################################################################
# Rule 1: Split the mixed genomes FASTA into individual FASTA files
# Input: Mixed genome fasta file
# Output: Split genome fasta files in the `individual_genomes` directory.
###############################################################################
rule split_genomes:
    input:
        "/data/qix007/DBkraken2/filtered_genomes.fa" # Replace with directory that contains the mixed genomes fasta file.
    output:
        "split_genomes.done"
    shell:
        """
        mkdir -p individual_genomes
        awk '/^>/{{ 
            if (f) close(f);
            id = substr($1, 2);
            gsub(/[^A-Za-z0-9._-]/, "_", id);
            f = "individual_genomes/" id ".fa";
        }} 
        {{ print > f }}' {input}
        touch {output}
        """
###############################################################################
# Rule 2: Simulate paired-end reads for each individual genome using ART
# Input: Fasta files that were split previously in the `individual_genomes` directory.
# Output: Simulated paired reads in the `art` directory.
###############################################################################
rule simulate_reads:
    input:
        "split_genomes.done"  # Wait for splitting to finish.
    output:
        done = touch("simulate_reads.done")
    shell:
        r"""
        mkdir -p art
        # Loop over each individual FASTA file and simulate reads.
        for fasta in individual_genomes/*.fa; do
            genome=$(basename "$fasta" .fa)
            art_illumina --paired --rndSeed 0 --noALN --maskN 0 --seqSys MSv3 \
                --len 150 --rcount 4000 -m 200 -s 10 --in "$fasta" --out art/"$genome"
        done
        touch {output.done}
        """

###############################################################################
# Rule 3: Merge all simulated read files into two combined FASTQ files.
# Input: Reads in the `art` directory.
# Output: Merged reads as `combined1_R1` or `combined1_R2` in the `output` directory.
###############################################################################
rule merge_reads:
    input:
        "simulate_reads.done"  # Ensure simulations are complete.
    output:
        R1 = "output/combined1_R1.fq",
        R2 = "output/combined1_R2.fq",
        done = touch("merge_reads.done")
    shell:
        r"""
        mkdir -p output
        # Merge the R1 files (assumes ART produces files ending with "1.fq")
        cat art/*1.fq > {output.R1}
        # Merge the R2 files (assumes ART produces files ending with "2.fq")
        cat art/*2.fq > {output.R2}
        touch {output.done}
        """

###############################################################################
# Rule 4: Run Kraken2 on the merged paired-end reads.
# Input: Combined reads in the `output` directory.
# Output: In the `kraken_reports` directory, Kraken output file `kraken_output2.txt`, and Kraken report `kraken_report2.txt`. Also outputs which reads were classified, `classified_reads2` and which were unclassified, `unclassified_reads2`.
###############################################################################
rule kraken:
    input:
        "merge_reads.done"  # Wait until merge is complete.
    output:
        report = "kraken_reports/kraken_report2.txt",
        kraken_out = "kraken_reports/kraken_output2.txt",
        # Kraken outputs paired files using a '#' placeholder in the name.
        done = touch("kraken.done")
    shell:
        r"""
        mkdir -p kraken_reports
        kraken2 --db /data/qix007/newDBkraken2 --threads 4 \
            --paired output/combined1_R1.fq output/combined1_R2.fq \
            --report {output.report} --output {output.kraken_out} \
            --classified-out classified_reads2#.fq --unclassified-out unclassified_reads2#.fq
        touch {output.done}
        """

###############################################################################
# Rule 5: Evaluate Kraken using the provided Python script.
# Input: Kraken report file, `kraken_report2.txt` in the `kraken_reports` directory.
# Output: `accuracy.txt` in the `kraken_reports` directory.
###############################################################################
rule evaluate_kraken:
    input:
        "kraken.done"  # Ensure Kraken is complete.
    output:
        accuracy = "kraken_reports/accuracy.txt",
        done = touch("evaluate_kraken.done")
    shell:
        r"""
        python kraken_reports/evaluate_kraken.py > {output.accuracy}
        touch {output.done}
        """

###############################################################################
# Rule 6: Split the FASTQ files by running the external script.
# Input: `combined1_R1`, and `combined1_R2` reads in the `output` directory. Definition of mapping is also needed. Example is listed in the README documentation.
# Output: split fastq files in the directory it is ran from.
###############################################################################
rule split_fastq:
    input:
        "evaluate_kraken.done"  # Wait until Kraken evaluation is done.
    output:
        done = touch("split_fq.done")
    shell:
        r"""
        bash split_fq.sh
        touch {output.done}
        """
