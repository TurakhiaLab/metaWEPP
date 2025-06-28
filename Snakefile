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
DIR = config.get("DIR")
if "DIR" not in config:
    raise ValueError(
        "No DIR folder specified.\n"
        "Call snakemake with, e.g.,  --config DIR=Sample_"
    )

WEPP_ROOT        = "WEPP"
WEPP_WORKFLOW    = "WEPP/workflow/Snakefile"
WEPP_RESULTS_DIR = "WEPP/results"
WEPP_DATA_DIR    = "WEPP/data"
WEPP_CONFIG      = "WEPP/config/config.yaml"
WEPP_DATA_PREFIX = f"{WEPP_DATA_DIR}/{DIR}"

SIM_TOOL         = config.get("SIMULATION_TOOL", "none").upper() 

KRAKEN_DB = config.get("KRAKEN_DB")
if KRAKEN_DB is None:
    raise ValueError(
        "Please provide the path to the Kraken2 database via\n"
        "  --config KRAKEN_DB=<folder>"
    )
TAXID_MAP = os.path.join(config["KRAKEN_DB"], "seqid2taxid.map")
IS_SINGLE_END = config.get("SEQUENCING_TYPE", "p").lower() in {"s", "n"}

# Get FQ1 and FQ2
if SIM_TOOL == "MESS":
    FQ1 = "simulated_reads/fastq/merged_R1.fq.gz"
    FQ2 = "simulated_reads/fastq/merged_R2.fq.gz"

else:
    fq_dir = Path("data") / config["DIR"]
    if not fq_dir.exists():
        raise FileNotFoundError(f"Input folder {DIR} does not exist")

    # always look for R1
    r1_files = sorted(fq_dir.glob("*.fastq*"))
    if len(r1_files) < 0:
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
        return f"results/{DIR}"
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

SAMPLE_TAG   = config["DIR"]        
def tag(acc):  # AF013254.1  -->  AF013254.1_RSVA_test
    return f"{acc}_{SAMPLE_TAG}"

# ────────────────────────────────────────────────────────────────
# Helper function for wepp done
def wepp_done_files():
    return [
        f"{WEPP_RESULTS_DIR}/{acc}/{acc}_run.txt"
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

    # Find the genomes file use for simulation
    genome_dir = Path("data") / config["DIR"]
    fasta_candidates = sorted(genome_dir.glob("*.fa"))

    if not fasta_candidates:
        sys.exit(
            f"No *.fa files found in {genome_dir}\n"
            "Place at least one multi-FASTA there or fix the DIR setting."
        )

    # use the first one that was found
    GENOMES_FASTA = str(fasta_candidates[0])

    print(f"Using genome FASTA for MeSS: {GENOMES_FASTA}", file=sys.stderr)

    # 4) Split genomes to individual fasta files to input to MeSS (their multi genome doesn't work for some reason)
    checkpoint split_genomes:
        input:
            fasta = lambda wc: GENOMES_FASTA
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
        report     = "kraken_report.txt",
        kraken_out = "kraken_output.txt",
    threads: config.get("kraken_threads", 8)
    params:
        db        = KRAKEN_DB,
        mode_flag = "--single" if IS_SINGLE_END else "--paired",
        mate2     = (lambda wc: "" if IS_SINGLE_END else FQ2),
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
        kraken_out = "kraken_output.txt",
        mapping    = TAXID_MAP,
        r1         = FQ1,
        r2         = (lambda wc: [] if IS_SINGLE_END else FQ2),
    output:
        dir = directory(OUT_ROOT)
    params:
        script = "scripts/split_read.py",
        r2_arg = (lambda wc: "" if IS_SINGLE_END else f"--r2 {FQ2}"),
        ref_arg = "--ref-accessions " + ",".join(REF_ACCESSIONS),
        dir    =  OUT_ROOT,
    shell:
        r"""
        python {params.script} \
            --kraken-out {input.kraken_out} \
            --mapping    {input.mapping}    \
            --r1         {input.r1}         \
            {params.r2_arg} \
            {params.ref_arg} \
            --out-dir    {params.dir}
        """

# 7.5) Prepare wepp input directory
# For every accession (wildcard: {acc}) copy the reference *.fa and the matching 
# *.pb.gz / *.mat into the same OUT_ROOT/<acc>/ directory that already contains 
# the split reads.
# helper defined
TAG = config["DIR"]     # e.g. "RSVA_test"
def tag(acc):           #  AF013254.1  →  AF013254.1_RSVA_test   (DIR = RSVA_test)
    return f"{acc}_{TAG}"

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
            new_r1    = f"{WEPP_DATA_DIR}/{{acc}}_{TAG}/{{acc}}.fastq.gz",
            fasta_out = f"{WEPP_DATA_DIR}/{{acc}}_{TAG}/{{acc}}.fa",
            pb_out    = f"{WEPP_DATA_DIR}/{{acc}}_{TAG}/{{acc}}.pb.gz"
        params:
            data_dir  = lambda wc: f"{WEPP_DATA_DIR}/{tag(wc.acc)}"
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
            fastq_dir = OUT_ROOT,          # just a string
            r1 = lambda wc: optional_file(
                    f"{OUT_ROOT}/{wc.acc}/{wc.acc.split('.')[0]}_R1.fq.gz"),
            r2 = lambda wc: optional_file(
                    f"{OUT_ROOT}/{wc.acc}/{wc.acc.split('.')[0]}_R2.fq.gz"),
            fasta = lambda wc: ACC2FASTA[wc.acc],
            pb    = lambda wc: ACC2PB[wc.acc]
        output:
            new_r1    = f"{WEPP_DATA_DIR}/{{acc}}_{TAG}/{{acc}}_R1.fastq.gz",
            new_r2    = f"{WEPP_DATA_DIR}/{{acc}}_{TAG}/{{acc}}_R2.fastq.gz",
            fasta_out = f"{WEPP_DATA_DIR}/{{acc}}_{TAG}/{{acc}}.fa",
            pb_out    = f"{WEPP_DATA_DIR}/{{acc}}_{TAG}/{{acc}}.pb.gz"
        params:
            data_dir  = lambda wc: f"{WEPP_DATA_DIR}/{tag(wc.acc)}"
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
        # FASTQs are stored in   WEPP/data/<acc>_<TAG>/
        r1 = lambda wc:
            f"{WEPP_DATA_DIR}/{tag(wc.acc)}/{wc.acc}.fastq.gz"
            if IS_SINGLE_END else
            f"{WEPP_DATA_DIR}/{tag(wc.acc)}/{wc.acc}_R1.fastq.gz",
        r2 = lambda wc: [] if IS_SINGLE_END else
            f"{WEPP_DATA_DIR}/{tag(wc.acc)}/{wc.acc}_R2.fastq.gz",
        fasta_full = lambda wc:
            f"{WEPP_DATA_DIR}/{tag(wc.acc)}/{wc.acc}.fa",
        tree_full  = lambda wc:
            f"{WEPP_DATA_DIR}/{tag(wc.acc)}/{wc.acc}.pb.gz",
    output:
        run_txt = f"{WEPP_RESULTS_DIR}/{{acc}}/{{acc}}_run.txt",
    threads: config.get("wepp_threads", 32)
    params:
        snakefile   = str(WEPP_WORKFLOW),
        workdir     = str(WEPP_ROOT),
        primer_bed  = config["PRIMER_BED"],
        clade_idx   = config["CLADE_IDX"],
        cfgfile     = str(WEPP_CONFIG),
        resultsdir  = str(WEPP_RESULTS_DIR),
        prefix      = lambda wc: wc.acc.split('.')[0],
        ref_name    = lambda wc: f"{wc.acc}.fa",
        tree_name   = lambda wc: f"{wc.acc}.pb.gz",
        seq_type    = config.get("SEQUENCING_TYPE", "p").lower(),
        tag_dir     = lambda wc: tag(wc.acc),          # <acc>_<TAG>
    conda:
        "env/wepp.yaml"
    shell:
        r"""
        # ── skip WEPP when both FASTQs are empty ──────────────────────
        has_reads() {{
            local f=$1
            [ "$( (gzip -cd "$f" 2>/dev/null || cat "$f") | head -4 | wc -l )" -ge 4 ]
        }}

        if ! has_reads "{input.r1}" && ( [ -z "{input.r2}" ] || ! has_reads "{input.r2}" ); then
            echo "No reads for {wildcards.acc}; removing data folder and skipping WEPP." >&2
            rm -rf "{WEPP_DATA_DIR}/{params.tag_dir}"
            mkdir -p {params.resultsdir}/{wildcards.acc}
            touch {output.run_txt}
            exit 0
        fi

        # ── run inner WEPP workflow ───────────────────────────────────
        mkdir -p {params.resultsdir}/{wildcards.acc}

        python ./scripts/run_inner.py \
            --snakefile  {params.snakefile} \
            --workdir    {params.workdir} \
            --dir        {params.tag_dir} \
            --prefix     {params.prefix} \
            --primer_bed {params.primer_bed} \
            --tree       {params.tree_name} \
            --ref        {params.ref_name} \
            --clade_idx  {params.clade_idx} \
            --configfile {params.cfgfile} \
            --cores      {threads} \
            --sequencing_type {params.seq_type}

        touch {output.run_txt}
        """