# Snakefile

# 1) Load METAWEPP config
configfile: "config/config.yaml"

# imports
import os
import csv, re, itertools

# 2) Constants from config
WEPP_DIR         = config["wepp_results_dir"]
SIM_TOOL         = config.get("simulation_tool", "none").upper() 
TAXIDS = config["target_taxids"]
DATA_DIR = config["wepp_data_dir"]
REF_BASENAME = os.path.basename(config["REF"])
PB_BASENAME = os.path.basename(config["TREE"])
TAXID_MAP = os.path.join(config["kraken_db"], "seqid2taxid.map")

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

# ────────────────────────────────────────────────────────────────
# Auto-discover reference genomes placed in:
#   data/pathogens_for_detailed_analysis/pathogen*/<anything>.fa*
# ────────────────────────────────────────────────────────────────
from pathlib import Path
import re

PATHOGEN_ROOT = Path("data/pathogens_for_detailed_analysis")

# Collect all *.fa / *.fasta files one level below each pathogenX dir
FASTAS = sorted(PATHOGEN_ROOT.glob("*/**/*.fa*"))

if not FASTAS:
    raise ValueError(
        f"No reference FASTA files were found under data/pathogens_for_detailed_analysis"
    )

# Full paths
REF_FASTA_FILES = [str(p) for p in FASTAS]

# Accessions (file stem); e.g. "AF013254.1" from ".../AF013254.1.fa"
REF_ACCESSIONS = sorted({p.parent.name for p in FASTAS})

# Optional: a mapping  accession → absolute path
ACC2FASTA = {p.stem: str(p.resolve()) for p in FASTAS}

# Make them visible in the log when Snakefile loads (optional)
print("Found reference genomes:", ", ".join(REF_ACCESSIONS))
# ────────────────────────────────────────────────────────────────


# ────────────────────────────────────────────────────────────────
# Helper function for splitting reads
# ────────────────────────────────────────────────────────────────
def detect_out_root(fq1):
    if "simulated_reads/fastq" in fq1:
        return "results/fastq"
    sample = os.path.basename(os.path.dirname(fq1)) # e.g. data/sample_1/…
    return f"results/{sample}"

OUT_ROOT = detect_out_root(FQ1)

def accessions_from_map(map_file, target_taxids):
    keep = set(str(t) for t in target_taxids)
    accs = []
    with open(map_file) as f:
        for acc, tax in (l.split() for l in f if l.strip()):
            if tax in keep:
                accs.append(acc)
    return accs

# ACCESSIONS = accessions_from_map(config["taxid_map"], config["target_taxids"]) 

def split_fastq_outputs():
    out_files = []
    for acc in REF_ACCESSIONS:
        pfx = acc.split(".")[0]
        out_files.extend([
            f"{OUT_ROOT}/{acc}/{pfx}_R1.fq.gz",
            f"{OUT_ROOT}/{acc}/{pfx}_R2.fq.gz",
        ])
    return out_files

SPLIT_FASTQS = split_fastq_outputs()
# ────────────────────────────────────────────────────────────────


# ────────────────────────────────────────────────────────────────
# Helper function to find .pb.gz to each reference pathogen
ACC2PB = {}
for acc, fasta_path in ACC2FASTA.items():
    pdir = Path(fasta_path).parent
    # allow either *.pb.gz or *.mat – first match wins
    pb_candidates = list(itertools.chain(
        pdir.glob("*.pb.gz"),
        pdir.glob("*.mat"),
        pdir.glob("*.pb")          # just in case
    ))
    if not pb_candidates:
        raise ValueError(f"No PB/MAT file found next to {fasta_path}")
    ACC2PB[acc] = str(pb_candidates[0])

def wepp_done_files():
    return [
        f"{config['wepp_results_dir']}/{acc}/{acc}_run.txt"
        for acc in REF_ACCESSIONS
    ]

# 3) Top‐level: these are all wildcards required for the pipeline to run smoothly
GENOME_NAME, = glob_wildcards("individual_genomes/{genome}.fa")

rule all:
    input:
        wepp_done_files(),
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
rule split_per_accession:
    input:
        kraken_out = config["kraken_output"],
        r1         = FQ1,
        r2         = FQ2,
        mapping    = TAXID_MAP
    output:
        fastqs = SPLIT_FASTQS
    params:
        script = "scripts/split_classified_reads_2.py",
        outdir = OUT_ROOT
    shell:
        """
        python {params.script} \
            --kraken-out {input.kraken_out} \
            --mapping {input.mapping} \
            --r1 {input.r1} \
            --r2 {input.r2} \
            --out-dir {params.outdir}
            
        """


# 7.5) Prepare wepp input directory
# For every accession (wildcard: {acc}) copy the reference *.fa and the matching 
# *.pb.gz / *.mat into the same OUT_ROOT/<acc>/ directory that already contains 
# the split reads.
rule prepare_wepp_inputs:
    """
    Copy R1/R2, reference FASTA, and PB/MAT tree for each accession
    from OUT_ROOT/{acc}/  →  DATA_DIR/{acc}/
    """
    input:
        r1    = lambda wc: f"{OUT_ROOT}/{wc.acc}/{wc.acc.split('.')[0]}_R1.fq.gz",
        r2    = lambda wc: f"{OUT_ROOT}/{wc.acc}/{wc.acc.split('.')[0]}_R2.fq.gz",
        fasta = lambda wc: ACC2FASTA[wc.acc],
        pb    = lambda wc: ACC2PB[wc.acc]
    output:
        new_r1    = DATA_DIR + "/{acc}/{acc}_R1.fastq.gz",
        new_r2    = DATA_DIR + "/{acc}/{acc}_R2.fastq.gz",
        fasta_out = DATA_DIR + "/{acc}/{acc}.fa",
        pb_out    = DATA_DIR + "/{acc}/{acc}.pb.gz"
    params:
        data_dir = lambda wc: f"{DATA_DIR}/{wc.acc}"
    shell:
        """
        mkdir -p {params.data_dir}
        # copy reference files
        cp {input.fasta} {output.fasta_out}
        cp {input.pb}    {output.pb_out}
        # copy split reads
        cp {input.r1}    {output.new_r1}
        cp {input.r2}    {output.new_r2}
        """



# 8) Invoke WEPP’s Snakefile for each taxid
# Run WEPP
rule run_wepp:
    """
    Launch the inner WEPP workflow once for every accession (`{acc}`).
    """
    input:
        r1    = lambda wc: f"{DATA_DIR}/{wc.acc}/{wc.acc.split('.')[0]}_R1.fastq.gz",
        r2    = lambda wc: f"{DATA_DIR}/{wc.acc}/{wc.acc.split('.')[0]}_R2.fastq.gz",
        fasta = lambda wc: f"{DATA_DIR}/{wc.acc}/{wc.acc}.fa",
        tree  = lambda wc: f"{DATA_DIR}/{wc.acc}/{wc.acc}.pb.gz",

    output:
        run_txt = config["wepp_results_dir"] + "/{acc}/{acc}_run.txt",

    threads: config.get("wepp_threads", 32)

    params:
        snakefile  = config["wepp_workflow"],
        workdir    = config["wepp_root"],
        primer_bed = config["PRIMER_BED"],
        clade_idx  = config["CLADE_IDX"],
        cfgfile    = config["wepp_config"],
        resultsdir = config["wepp_results_dir"],

    shell:
        r"""
        # ensure results folder exists
        mkdir -p {params.resultsdir}/{wildcards.acc}

        # strip version suffix once here
        PREFIX=$(echo "{wildcards.acc}" | cut -d '.' -f 1)

        # run the inner Snakemake workflow
        python ./scripts/run_inner.py \
            --snakefile  {params.snakefile} \
            --workdir    {params.workdir}  \
            --dir        {wildcards.acc}   \
            --prefix     ${{PREFIX}}       \
            --primer_bed {params.primer_bed} \
            --tree       {input.tree}      \
            --ref        {input.fasta}     \
            --clade_idx  {params.clade_idx} \
            --configfile {params.cfgfile}  \
            --cores      {threads}

        # mark completion for the outer DAG
        touch {output.run_txt}
        """