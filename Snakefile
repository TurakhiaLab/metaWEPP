# Snakefile

# imports
import os
import csv, re, itertools

from pathlib import Path
import sys, glob

from collections import defaultdict
import itertools

# 1) Load METAWEPP config
configfile: "config/config.yaml"

# 2) Constants from config
WEPP_DIR         = config["wepp_results_dir"]
SIM_TOOL         = config.get("SIMULATION_TOOL", "none").upper() 
DATA_DIR = config["wepp_data_dir"]
#REF_BASENAME = os.path.basename(config["REF"])
#PB_BASENAME = os.path.basename(config["TREE"])
TAXID_MAP = os.path.join(config["KRAKEN_DB"], "seqid2taxid.map")
IS_SINGLE_END = config.get("IS_SINGLE_END", "p").lower() == "s"


# Get FQ1 and FQ2
if SIM_TOOL == "MESS":
    FQ1 = "simulated_reads/fastq/merged_R1.fq.gz"
    FQ2 = "simulated_reads/fastq/merged_R2.fq.gz"

else:
    if "fq_dir" not in config:
        raise ValueError(
            "SIM_TOOL is NONE, but config is missing 'fq_dir'.\n"
            "Please add, e.g.,\n"
            "  fq_dir: RSVA_test   # relative to the data/ folder\n"
            "to your config.yaml."
        )

    fq_dir = Path("data") / config["FQ_DIR"]
    if not fq_dir.exists():
        raise FileNotFoundError(f"Input folder {fq_dir} does not exist")

    # always look for R1
    r1_files = sorted(fq_dir.glob("*.fastq*"))
    if len(r1_files) != 1:
        raise RuntimeError(
            f"{fq_dir} must contain exactly one *.fastq* file "
            f"(found {len(r1_files)})."
        )

    # look for R2 only when paired-end
    if IS_SINGLE_END:
        r2_files = []
        FQ2 = ""
    else:
        r2_files = sorted(fq_dir.glob("*_R2.fastq*"))
        if len(r2_files) != 1:
            raise RuntimeError(
                f"{fq_dir} must contain exactly one *_R2.fastq* file "
                f"(found {len(r2_files)})."
            )
        FQ2 = str(r2_files[0])

    FQ1 = str(r1_files[0])

print(f"Input FASTQs chosen:\n  FQ1 = {FQ1}\n  FQ2 = {FQ2}", file=sys.stderr)


# Helper function for the merge reads rule:

# def existing_merge_genomes():
#     all_genomes, = glob_wildcards("simulated_reads/fastq/{genome}_R1.fq.gz")
#     return [g for g in all_genomes if os.path.exists(f"simulated_reads/fastq/{g}_R2.fq.gz")]

#MERGE_GENOMES = existing_merge_genomes()

# ────────────────────────────────────────────────────────────────
# Auto-discover reference genomes placed in:
#   data/pathogens_for_detailed_analysis/pathogen*/<anything>.fa*
# ────────────────────────────────────────────────────────────────
PATHOGEN_ROOT = Path("data/pathogens_for_detailed_analysis")

# find any *.fa / *.fasta below the first directory level
FASTAS = sorted(PATHOGEN_ROOT.glob("*/*.fa*"))

if not FASTAS:
    raise ValueError("No reference FASTA files found under data/pathogens_for_detailed_analysis")

# build:  accession list,  accession to fasta,  accession to pb ----
REF_ACCESSIONS = sorted({fa.parent.name for fa in FASTAS})      # ['AF013254.1', …]

ACC2FASTA = {}          # accession to single FASTA path
for fa in FASTAS:
    acc = fa.parent.name
    ACC2FASTA.setdefault(acc, str(fa.resolve()))   

ACC2PB = {}
for acc in REF_ACCESSIONS:
    pdir = PATHOGEN_ROOT / acc
    cand = next(itertools.chain(
        pdir.glob("*.pb.gz"),
        pdir.glob("*.mat"),
        pdir.glob("*.pb"),
    ), None)
    if cand is None:
        raise ValueError(f"No *.pb.gz or *.mat found in {pdir}")
    ACC2PB[acc] = str(cand.resolve())

print("Found reference folders:", ", ".join(REF_ACCESSIONS))
print(f"IS_SINGLE_END = {IS_SINGLE_END}", file=sys.stderr)

# ────────────────────────────────────────────────────────────────
# Helper function for splitting reads
def detect_out_root(fq1):
    if "simulated_reads/fastq" in fq1:
        return "results/fastq"
    sample = os.path.basename(os.path.dirname(fq1)) # e.g. data/sample_1/…
    return f"results/{sample}"

OUT_ROOT = detect_out_root(FQ1)


def split_fastq_outputs():
    out_files = []
    for acc in REF_ACCESSIONS:
        pfx = acc.split(".")[0]
        out_files.extend([
            f"{OUT_ROOT}/{acc}/{pfx}_R1.fq.gz",
            #f"{OUT_ROOT}/{acc}/{pfx}_R2.fq.gz",
        ])
    return out_files

SPLIT_FASTQS = split_fastq_outputs()


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

# ────────────────────────────────────────────────────────────────
# Helper function for prepare wepp
def optional_file(path):
    """
    Return the path as a string if it exists,
    otherwise return an empty list (→ Snakemake ignores it).
    """
    return path if os.path.exists(path) else []

# ────────────────────────────────────────────────────────────────
# Helper function for wepp done
def wepp_done_files():
    return [
        f"{config['wepp_results_dir']}/{acc}/{acc}_run.txt"
        for acc in REF_ACCESSIONS
    ]

# ────────────────────────────────────────────────────────────────
# Helper function for run wepp
def has_reads(fq):
    """
    Return shell code that exits 0 (true) when the FASTQ has ≥ 1 read.
    Works for .gz and plain fastq. 4 lines = 1 read.
    """
    return f'( (gzip -cd {fq} 2>/dev/null || cat {fq}) | head -4 | wc -l ) -ge 4'


# 3) Top‐level: these are all wildcards required for the pipeline to run smoothly
GENOME_NAME, = glob_wildcards("individual_genomes/{genome}.fa")

ALL_INPUT = wepp_done_files() + [FQ1] + ([] if IS_SINGLE_END else [FQ2])
rule all:
    input: ALL_INPUT

if SIM_TOOL == "MESS":

    # 4) Split genomes to individual fasta files to input to MeSS (their multi genome doesn't work for some reason)

    checkpoint split_genomes:
        input:
            fasta = config["METAGENOMIC_REF"]
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
            coverage = config["COVERAGE"]
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
            source $(conda info --base)/etc/profile.d/conda.sh
            conda activate mess

            mess simulate \
                --input {input.tsv} \
                --tech illumina \
                --output simulated_reads/ \
                --fasta individual_genomes

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

# ─── helper function for names ──────────────────────────────────────────
FQ1_GZ = FQ1 if FQ1.endswith(".gz") else FQ1 + ".gz"
FQ2_GZ = (
    "" if IS_SINGLE_END
    else (FQ2 if FQ2.endswith(".gz") else FQ2 + ".gz")
)

# Make sure everything is compressed (Compresses fq files, turns them into fastq.gz files, 
# removes fq files because WEPP expects a specific length of fq files)
if SIM_TOOL != "MESS":

    if IS_SINGLE_END:
        rule compress_fastqs_for_wepp_se:
            input:
                r1 = FQ1
            output:
                r1_gz = FQ1_GZ
            shell:
                """
                mkdir -p $(dirname {output.r1_gz})
                if [ ! -f {output.r1_gz} ]; then gzip -c {input.r1} > {output.r1_gz}; fi
                """
    else:
        rule compress_fastqs_for_wepp_pe:
            input:
                r1 = FQ1,
                r2 = FQ2
            output:
                r1_gz = FQ1_GZ,
                r2_gz = FQ2_GZ
            shell:
                """
                mkdir -p $(dirname {output.r1_gz})
                if [ ! -f {output.r1_gz} ]; then gzip -c {input.r1} > {output.r1_gz}; fi
                if [ ! -f {output.r2_gz} ]; then gzip -c {input.r2} > {output.r2_gz}; fi
                """

# 6) Kraken2 classification 
rule kraken:
    input:
        r1 = FQ1,
        r2 = (lambda wc: [] if IS_SINGLE_END else FQ2),
    output:
        report     = config["KRAKEN_REPORT"],
        kraken_out = config["KRAKEN_OUTPUT"],
    threads: config.get("kraken_threads", 4)
    params:
        db        = config["KRAKEN_DB"],
        mode_flag = "--single" if IS_SINGLE_END else "--paired",
        mate2     = (lambda wc: "" if IS_SINGLE_END else FQ2),   # ← FIX
    shell:
        r"""
        mkdir -p $(dirname {output.report})
        kraken2 --db {params.db} --threads {threads} {params.mode_flag} \
                {input.r1} {params.mate2} \
                --report {output.report} \
                --output {output.kraken_out}
        """

# 7) Split Kraken output into per-taxid FASTQs 
rule split_per_accession:
    input:
        kraken_out = config["KRAKEN_OUTPUT"],
        mapping    = TAXID_MAP,
        r1         = FQ1,
        r2         = (lambda wc: [] if IS_SINGLE_END else FQ2),
    output:
        dir = directory(OUT_ROOT)
    params:
        script = "scripts/split_read.py",
        r2_arg = (lambda wc: "" if IS_SINGLE_END else f"--r2 {FQ2}"),
        dir    =  OUT_ROOT,
    shell:
        r"""
        python {params.script} \
            --kraken-out {input.kraken_out} \
            --mapping    {input.mapping}    \
            --r1         {input.r1}         \
            {params.r2_arg} \
            --out-dir    {params.dir}
        """

# 7.5) Prepare wepp input directory
# For every accession (wildcard: {acc}) copy the reference *.fa and the matching 
# *.pb.gz / *.mat into the same OUT_ROOT/<acc>/ directory that already contains 
# the split reads.
if IS_SINGLE_END:
    # ───────── single-end ────────────────────────────────────────────────────
    rule prepare_wepp_inputs:
        input:
            fastq_dir = directory(OUT_ROOT),     # ensures split_per_accession ran
            r1    = lambda wc: optional_file(
                    f"{OUT_ROOT}/{wc.acc}/{wc.acc.split('.')[0]}_R1.fq.gz"),
            fasta = lambda wc: ACC2FASTA[wc.acc],
            pb    = lambda wc: ACC2PB[wc.acc]
        output:
            new_r1    = DATA_DIR + "/{acc}/{acc}.fastq.gz",
            fasta_out = DATA_DIR + "/{acc}/{acc}.fa",
            pb_out    = DATA_DIR + "/{acc}/{acc}.pb.gz"
        params:
            data_dir = lambda wc: f"{DATA_DIR}/{wc.acc}"
        shell:
            r"""
            mkdir -p {params.data_dir}
            cp {input.fasta} {output.fasta_out}
            cp {input.pb}    {output.pb_out}

            if [ -f "{input.r1}" ]; then
                cp {input.r1} {output.new_r1}
            else
                echo | gzip -c > {output.new_r1}
            fi
            """


else:
    # ───────── paired-end ───────────────────────────────────────────────────
    rule prepare_wepp_inputs:
        input:
            fastq_dir = directory(OUT_ROOT),     # ensures split_per_accession ran
            r1  = lambda wc: optional_file(
                    f"{OUT_ROOT}/{wc.acc}/{wc.acc.split('.')[0]}_R1.fq.gz"),
            r2  = lambda wc: optional_file(
                    f"{OUT_ROOT}/{wc.acc}/{wc.acc.split('.')[0]}_R2.fq.gz"),
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
            r"""
            mkdir -p {params.data_dir}
            cp {input.fasta} {output.fasta_out}
            cp {input.pb}    {output.pb_out}

            if [ -f "{input.r1}" ]; then
                cp {input.r1} {output.new_r1}
            else
                echo | gzip -c > {output.new_r1}
            fi

            if [ -f "{input.r2}" ]; then
                cp {input.r2} {output.new_r2}
            else
                echo | gzip -c > {output.new_r2}
            fi
            """


# 8) Invoke WEPP’s Snakefile for each taxid
# Run WEPP
rule run_wepp:
    input:
        r1  = lambda wc: (
                f"{DATA_DIR}/{wc.acc}/{wc.acc}.fastq.gz"       
                if IS_SINGLE_END else
                f"{DATA_DIR}/{wc.acc}/{wc.acc}_R1.fastq.gz"), 
        r2  = (lambda wc: [] if IS_SINGLE_END else 
                f"{DATA_DIR}/{wc.acc}/{wc.acc}_R2.fastq.gz"),
        fasta_full = lambda wc: f"{DATA_DIR}/{wc.acc}/{wc.acc}.fa",
        tree_full  = lambda wc: f"{DATA_DIR}/{wc.acc}/{wc.acc}.pb.gz",
    output:
        run_txt = config["wepp_results_dir"] + "/{acc}/{acc}_run.txt",
    threads: config.get("wepp_threads", 32)
    params:
        snakefile   = config["wepp_workflow"],
        workdir     = config['wepp_root'],
        primer_bed  = config["PRIMER_BED"],
        clade_idx   = config["CLADE_IDX"],
        cfgfile     = config["wepp_config"],
        resultsdir  = config["wepp_results_dir"],
        prefix      = lambda wc: wc.acc.split('.')[0],
        ref_name    = lambda wc: f"{wc.acc}.fa",
        tree_name   = lambda wc: f"{wc.acc}.pb.gz",
        se_flag     = "--is_single_end" if IS_SINGLE_END else "", 
    conda:
        "/home/qix007@AD.UCSD.EDU/metagenomic-WBE/env/wepp.yaml"
    shell:
        r"""
        # ─── skip WEPP when both FASTQs are empty ──────────────────────
        has_reads() {{
            local f=$1
            [ "$( (gzip -cd "$f" 2>/dev/null || cat "$f") | head -4 | wc -l )" -ge 4 ]
        }}

        if ! has_reads "{input.r1}" && ( [ -z "{input.r2}" ] || ! has_reads "{input.r2}" ); then
            echo "No reads for {wildcards.acc}; removing data folder and skipping WEPP." >&2
            rm -rf "{DATA_DIR}/{wildcards.acc}"
            mkdir -p {params.resultsdir}/{wildcards.acc}
            touch {output.run_txt}
            exit 0
        fi
        # ─── continue with normal WEPP run ─────────────────────────────
        mkdir -p {params.resultsdir}/{wildcards.acc}

        python ./scripts/run_inner.py \
            --snakefile  {params.snakefile} \
            --workdir    {params.workdir} \
            --dir        {wildcards.acc} \
            --prefix     {params.prefix} \
            --primer_bed {params.primer_bed} \
            --tree       {params.tree_name} \
            --ref        {params.ref_name} \
            --clade_idx  {params.clade_idx} \
            --configfile {params.cfgfile} \
            --cores      {threads} \
            {params.se_flag}

        touch {output.run_txt}
        """
