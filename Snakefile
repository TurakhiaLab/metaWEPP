# Snakefile

# 1) Load METAWEPP config
configfile: "config/config.yaml"

# imports
import os

# 2) Constants from config
WEPP_DIR         = config["wepp_results_dir"]
SIM_TOOL         = config.get("simulation_tool", "none").upper() 
TAXIDS = config["target_taxids"]
DATA_DIR = config["wepp_data_dir"]
REF_BASENAME = os.path.basename(config["REF"])
PB_BASENAME = os.path.basename(config["TREE"])

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
if SIM_TOOL in ("MESS"):
    FQ1 = "simulated_reads/fastq/merged_R1.fq.gz"
    FQ2 = "simulated_reads/fastq/merged_R2.fq.gz"
else:
    FQ1 = config["fq1"]
    FQ2 = config["fq2"]

# Helper function for the merge reads rule:

# def existing_merge_genomes():
#     all_genomes, = glob_wildcards("simulated_reads/fastq/{genome}_R1.fq.gz")
#     return [g for g in all_genomes if os.path.exists(f"simulated_reads/fastq/{g}_R2.fq.gz")]

#MERGE_GENOMES = existing_merge_genomes()

# 3) Top‐level: these are all wildcards required for the pipeline to run smoothly
GENOME_NAME, = glob_wildcards("individual_genomes/{genome}.fa")

rule all:
    input:
        # expand(
        #     "{wepp_dir}/{dataset}/{dataset}_run.txt",
        #     wepp_dir = WEPP_DIR,
        #     dataset  = TAXIDS,
        # ),
        FQ1,
        FQ2

if SIM_TOOL == "MESS":

    # 4) Split genomes to individual fasta files to input to MeSS (their multi genome doesn't work for some reason)

    checkpoint split_genomes:
        input:
            fasta = config["mixed_genomes_fasta"]
        output:
            directory("individual_genomes")
        shell:
            """
            mkdir -p {output}
            awk '/^>/{{split($1, id, ">"); fname="{output}/"id[2]".fa"; print > fname; next}} {{print >> fname}}' {input.fasta}
            """
            
    def get_genomes_from_checkpoint(wc):
        ckpt = checkpoints.split_genomes.get(**wc)
        genomes, = glob_wildcards(os.path.join(ckpt.output[0], "{genome}.fa"))
        return genomes

    # 5) Optional MeSS simulation

    # Generate TSVs for MeSS:
    rule generate_mess_tsv:
        input:
            fasta = "individual_genomes/{genome}.fa"
        output:
            tsv = "individual_genomes_tsvs/{genome}.tsv"
        params:
            script = "scripts/generate_mess_tsv.py",
            coverage = config["coverage"]
        shell:
            """
            python {params.script} {output.tsv} {input.fasta} {params.coverage}
            """


            
    rule simulate_reads_mess:
        input:
            tsv = "individual_genomes_tsvs/{genome}.tsv"
        output:
            r1 = "simulated_reads/fastq/{genome}_R1.fq.gz", 
            r2 = "simulated_reads/fastq/{genome}_R2.fq.gz", 
        params:
            genome_base = lambda wildcards: wildcards.genome.split(".")[0]
        resources:
            mess_slots=1
        shell:
            """
            mess simulate \
                --input {input.tsv} \
                --tech illumina \
                --output simulated_reads/ \
                --fasta ./individual_genomes

            mv simulated_reads/fastq/{params.genome_base}_R1.fq.gz {output.r1}
            mv simulated_reads/fastq/{params.genome_base}_R2.fq.gz {output.r2}
            """

    rule merge_reads:
        input:
            sim_reads_1 = lambda wc: expand("simulated_reads/fastq/{genome}_R1.fq.gz", genome=get_genomes_from_checkpoint(wc)),
            sim_reads_2 = lambda wc: expand("simulated_reads/fastq/{genome}_R2.fq.gz", genome=get_genomes_from_checkpoint(wc))
        output:
            R1 = FQ1,
            R2 = FQ2

        shell:
            """
            zcat simulated_reads/fastq/*_R1.fq.gz | gzip > {output.R1}
            zcat simulated_reads/fastq/*_R2.fq.gz | gzip > {output.R2}
            """


# Make sure everything is compressed (Compresses fq files, turns them into fastq.gz files, removes fq files because WEPP expects a specific length of fq files)
if SIM_TOOL != "MESS":
    def resolve_fastq(wildcards, read):
        fq = f"{wildcards.taxid}_{read}.fq"
        fq_gz = f"{wildcards.taxid}_{read}.fq.gz"
        return fq_gz if os.path.exists(fq_gz) else fq

    rule compress_fastqs_for_wepp:
        input:
            r1=lambda wc: resolve_fastq(wc, "R1"),
            r2=lambda wc: resolve_fastq(wc, "R2")
        output:
            r1_gz=FQ1,
            r2_gz=FQ2
        run:
            import os
            import shutil

            for in_f, out_f in zip([input.r1, input.r2], [output.r1_gz, output.r2_gz]):
                if os.path.abspath(in_f) == os.path.abspath(out_f):
                    continue  # Already correct format
                elif in_f.endswith(".gz"):
                    shutil.copy(in_f, out_f)
                else:
                    shell(f"gzip -c {in_f} > {out_f}")
                    os.remove(in_f)

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
        expand("{data_dir}/{taxid}/{taxid}_R1.fq.gz", data_dir=config["wepp_data_dir"], taxid=config["target_taxids"]),
        expand("{data_dir}/{taxid}/{taxid}_R2.fq.gz", data_dir=config["wepp_data_dir"], taxid=config["target_taxids"])
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
# Copies necessary fasta file and MAT into WEPP's DIR (the data directory) as it is required for WEPP to run properly

rule prepare_wepp_inputs:
    input:
        r1 = "{data_dir}/{taxid}/{taxid}_R1.fq.gz",
        r2 = "{data_dir}/{taxid}/{taxid}_R2.fq.gz",
        fasta = config["REF"],
        pb = config["TREE"]
    output:
        fasta_out = "{data_dir}/{taxid}/" + REF_BASENAME,
        pb_out    = "{data_dir}/{taxid}/" + PB_BASENAME
    params:
        data_dir = config["wepp_data_dir"]
    shell:
        """
        cp {input.fasta} {output.fasta_out}
        cp {input.pb} {output.pb_out}
        """

# 8) Invoke WEPP’s Snakefile for each taxid

# Run WEPP
rule run_wepp:
    input:
        r1    = f"{DATA_DIR}" + "/{taxid}/{taxid}_R1.fq.gz",
        r2    = f"{DATA_DIR}" + "/{taxid}/{taxid}_R2.fq.gz"
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
