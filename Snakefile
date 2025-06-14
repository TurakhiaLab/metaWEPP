# Snakefile

# 1) Load METAWEPP config
configfile: "config.yaml"

# imports
import os

# 2) Constants from config
WEPP_DIR         = config["wepp_results_dir"]
SIM_TOOL         = config.get("simulation_tool", "none").upper() 
TAXIDS = config["target_taxids"]
DATA_DIR = config["wepp_data_dir"]

# Error checking
if SIM_TOOL not in ("ART", "MESS", "NONE"):
    raise ValueError(f"Unknown SIM_TOOL: {SIM_TOOL}. Must be 'ART', 'MESS', or 'NONE'.")

if SIM_TOOL == "NONE":
    if "fq1" not in config or "fq2" not in config:
        raise ValueError(
            "SIM_TOOL is NONE, but config is missing 'fq1' or 'fq2'.\n"
            "Please provide real FASTQ files in your config.yaml like:\n"
            "  fq1: data/sample_R1.fastq\n"
            "  fq2: data/sample_R2.fastq"
        )

# Automatically define where merged reads should go
if SIM_TOOL in ("ART", "MESS"):
    FQ1 = "simulated/merged_R1.fq"
    FQ2 = "simulated/merged_R2.fq"
else:
    FQ1 = config["fq1"]
    FQ2 = config["fq2"]

# 3) Top‐level: expect one WEPP/run.txt per taxid under WEPP_DIR

rule all:
    input:
        expand(
            "{wepp_dir}/{dataset}/{dataset}_run.txt",
            wepp_dir = WEPP_DIR,
            dataset  = TAXIDS
        )

# 4) Optional ART simulation

# Rule to split genomes for ART:
rule split_genomes:
    input:
        fasta = config["mixed_genomes_fasta"]
    output:
        expand("individual_genomes/{genome}.fa", genome=GENOMES)
    shell:
        """
        mkdir -p individual_genomes
        awk '/^>/{gsub(/^>/,""); fname="individual_genomes/"$1".fa"; print > fname; next} {print >> fname}' {input.fasta}
        """

if SIM_TOOL == "ART":
    rule simulate_reads_art:
        input:
            genomes = rules.split_genomes.output
        output:
            expand("art/{genome}1.fq", genome=GENOMES),
            expand("art/{genome}2.fq", genome=GENOMES)
        shell:
            """
            mkdir -p art
            for fasta in {input.genomes}; do
                genome=$(basename "$fasta" .fa)
                art_illumina --paired --rndSeed 0 --noALN --maskN 0 --seqSys MSv3 \
                            --len 150 --rcount 800000 -m 200 -s 10 \
                            --in "$fasta" --out art/"$genome"
            done
            """

    merge_input = {
        "R1s": expand("art/{genome}1.fq", genome=GENOMES),
        "R2s": expand("art/{genome}2.fq", genome=GENOMES)
    }

# 5) Optional MeSS simulation

# Generate TSVs for MeSS:
rule generate_mess_tsv:
    input:
        fasta = "individual_genomes/{genome}.fa"
    output:
        tsv = "individual_genomes/{genome}.tsv"
    params:
        script = "scripts/generate_mess_tsv.py",
        coverage = config["coverage"]
    shell:
        """
        python {params.script} {output.tsv} {input.fasta} {params.coverage} {wildcards.genome}
        """

if SIM_TOOL == "MESS":
    rule simulate_reads_mess:
        input:
            expand("individual_genomes/{genome}.tsv", genome=GENOMES)
        output:
            expand("simulated_reads/{genome}_R1.fq", genome=GENOMES),
            expand("simulated_reads/{genome}_R2.fq", genome=GENOMES)
        shell:
            """
            mkdir -p simulated_reads
            for tsv in individual_genomes/*.tsv; do
                mess simulate --input "$tsv" \
                            --tech illumina \
                            --output simulated_reads/ \
                            --fasta individual_genomes/
            done
            """

    merge_input = {
        "R1s": expand("simulated_reads/{genome}_R1.fq", genome=GENOMES),
        "R2s": expand("simulated_reads/{genome}_R2.fq", genome=GENOMES)
    }

# Merge reads
if SIM_TOOL in ("ART", "MESS"):
    rule merge_reads:
        input:
            **merge_input
        output:
            R1 = FQ1,
            R2 = FQ2
        shell:
            """
            mkdir -p $(dirname {output.R1})
            cat {input.R1s} > {output.R1}
            cat {input.R2s} > {output.R2}
            """

# 6) Kraken2 classification 
rule kraken:
    input:
        r1 = FQ1,
        r2 = FQ2
    output:
        report     = config["kraken_report"],
        kraken_out = config["kraken_output"]
    threads: config.get("kraken_threads", 4)
    params:
        db = config["kraken_db"]
    shell:
        """
        mkdir -p $(dirname {output.report})
        kraken2 --db {params.db} --threads {threads} \
                --paired {input.r1} {input.r2} \
                --report {output.report} \
                --output {output.kraken_out}
        """

# 7) Split Kraken output into per‐taxid FASTQs
rule split_per_taxid:
    input:
        kraken_out = config["kraken_output"],
        r1         = FQ1,
        r2         = FQ2
    output:
        # one pair per taxid
        expand("{data_dir}/{taxid}/{taxid}_R1.fq", data_dir=config["wepp_data_dir"], taxid=config["target_taxids"]),
        expand("{data_dir}/{taxid}/{taxid}_R2.fq", data_dir=config["wepp_data_dir"], taxid=config["target_taxids"])
    params:
        script = "scripts/split_classified_reads.py",
        outdir = config["wepp_data_dir"]
    shell:
        """
        python {params.script} \
          --kraken-out {input.kraken_out} \
          --r1 {input.r1} \
          --r2 {input.r2} \
          --out-dir {params.outdir} 
        """

# 7.5) Prepare wepp input directory
# Copies necessary fasta file and MAT into WEPP's DIR (the data directory) as it is required for WEPP to run properly.
rule prepare_wepp_inputs:
    input:
        r1 = "{data_dir}/{taxid}/{taxid}_R1.fq",
        r2 = "{data_dir}/{taxid}/{taxid}_R2.fq",
        fasta = config["REF"],
        pb = config["TREE"]
    output:
        fasta_out = "{data_dir}/{taxid}/" + os.path.basename(config["REF"]),
        pb_out    = "{data_dir}/{taxid}/" + os.path.basename(config["TREE"])
    params:
        data_dir = config["wepp_data_dir"]
    shell:
        """
        cp {input.fasta} {output.fasta_out}
        cp {input.pb} {output.pb_out}
        """

# 8) Invoke WEPP’s Snakefile for each taxid

# Make sure everything is compressed (Compresses fq files, turns them into fastq.gz files, removes fq files because WEPP expects a specific length of fq files)
rule compress_fastqs_for_wepp:
    input:
        r1 = f"{DATA_DIR}" + "/{taxid}/{taxid}_R1.fq",
        r2 = f"{DATA_DIR}" + "/{taxid}/{taxid}_R2.fq"
    output:
        r1_gz = f"{DATA_DIR}" + "/{taxid}/{taxid}_R1.fastq.gz",
        r2_gz = f"{DATA_DIR}" + "/{taxid}/{taxid}_R2.fastq.gz"
    run:
        import os
        import shutil
        for in_f, out_f in zip([input.r1, input.r2], [output.r1_gz, output.r2_gz]):
            if not os.path.exists(out_f):
                if in_f.endswith(".gz"):
                    shutil.copy(in_f, out_f)
                else:
                    shell(f"gzip -c {in_f} > {out_f}")
            # After confirming output exists, remove original .fq
            if os.path.exists(out_f) and os.path.exists(in_f) and not in_f.endswith(".gz"):
                os.remove(in_f)

# Run WEPP
rule run_wepp:
    input:
        r1    = f"{DATA_DIR}" + "/{taxid}/{taxid}_R1.fastq.gz",
        r2    = f"{DATA_DIR}" + "/{taxid}/{taxid}_R2.fastq.gz",
        ready = f"{DATA_DIR}" + "/{taxid}/wepp_input_ready.txt"
    output:
        run_txt = config["wepp_results_dir"] + "/{taxid}/{taxid}_run.txt"
    threads: config.get("wepp_threads", 32)
    params:
        dir    = lambda wildcards: f"{wildcards.taxid}",
        prefix = lambda wildcards: f"{wildcards.taxid}"
    shell:
        """
        python ./scripts/run_inner.py \
            --snakefile {config[wepp_workflow]} \
            --workdir {config[wepp_root]} \
            --dir {params.dir} \
            --prefix {params.prefix} \
            --primer_bed {config[PRIMER_BED]} \
            --tree {config[TREE]} \
            --ref {config[REF]} \
            --clade_idx {config[CLADE_IDX]} \
            --configfile {config[wepp_config]} \
            --cores {threads}
        """
