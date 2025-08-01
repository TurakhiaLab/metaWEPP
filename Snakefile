# Snakefile Before keep pd name

# imports
import os
import csv, re, itertools

from pathlib import Path
import sys, glob

from collections import defaultdict
import itertools, json

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

# Original map file path
original_map = os.path.join(config["KRAKEN_DB"], "seqid2taxid.map")

# New map output path
TAXID_MAP = os.path.join("config", "accession2taxid.map")

# Ensure the output directory exists
os.makedirs("config", exist_ok=True)

# Convert the map
with open(original_map) as infile, open(TAXID_MAP, "w") as outfile:
    for line in infile:
        if line.strip():
            full_id, taxid = line.strip().split()
            accession = full_id.split("|")[-1]
            outfile.write(f"{accession}\t{taxid}\n")

IS_SINGLE_END = config.get("SEQUENCING_TYPE", "d").lower() in {"s", "n"}

# Get FQ1 and FQ2
if SIM_TOOL == "MESS":
    FQ1 = Path("data") / config["DIR"] / "merged_R1.fq.gz"
    FQ2 = Path("data") / config["DIR"] / "merged_R2.fq.gz"
    #"simulated_reads/fastq/merged_R1.fq.gz"

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
#   data/pathogens_for_wepp/pathogen*/<anything>.fa*
# ────────────────────────────────────────────────────────────────
PATHOGEN_ROOT = Path("data/pathogens_for_wepp")

# find any *.fa / *.fasta below the first directory level
FASTAS = sorted(PATHOGEN_ROOT.glob("*/*.[fF][aAn]*"))

if not FASTAS:
    print("\nWARNING: No pathogens given for variant analysis with WEPP!\n")

# build:  accession list,  accession to fasta,  accession to pb ----

ACC2FASTA, ACC2PB, ACC2DIR, ACC2CONFIG = {}, {}, {}, {}

for fa in FASTAS:                    #   data/pathogens_for_wepp/<dir>/<accession>.fa
    with open(fa) as fh:
        header = fh.readline().strip()
    path_dir = fa.parent 
    acc = header[1:].split()[0] 
    dir_name = fa.parent.name 
    ACC2FASTA[acc] = str(fa.resolve())  # Path to pathogen fa file 

    pb_file = next(itertools.chain(
        fa.parent.glob("*.pb.gz"),
        fa.parent.glob("*.pb")
    ), None)
    if pb_file is None:
        raise ValueError(f"No *.pb.gz / *.pb next to {fa}")
    ACC2PB[acc] = str(pb_file.resolve())

    # Find config.yaml for each pathogen
    config_file = path_dir / "config.yaml"
    ACC2CONFIG[acc] = str(config_file.resolve())

# Generate acc2dirname.json
for fasta in FASTAS:
    dir_name = fasta.parent.name
    with open(fasta) as f:
        header = f.readline().strip()
        if not header.startswith(">"):
            raise ValueError(f"{fasta} does not start with a FASTA header")
        acc = header[1:].split()[0]
        ACC2DIR[acc] = dir_name

# write JSON for use in split_read.py
Path("config").mkdir(exist_ok=True)
ACC2DIR_JSON = "config/acc2dirname.json"
with open(ACC2DIR_JSON, "w") as f:
    json.dump(ACC2DIR, f, indent=2)

# Extract accession IDs from .fa files' headers
REF_ACCESSIONS = []
for fasta in FASTAS:
    with open(fasta) as f:
        header = f.readline().strip()
        if not header.startswith(">"):
            raise ValueError(f"{fasta} does not start with a FASTA header")
        acc = header[1:].split()[0]
        REF_ACCESSIONS.append(acc)

print("Headers of reference fa file:", ", ".join(REF_ACCESSIONS))
print("Found accession → folder mappings:")
for acc, dirname in ACC2DIR.items():
    print(f"  {acc} -> {dirname}")
print(f"IS_SINGLE_END = {IS_SINGLE_END}", file=sys.stderr)

# ACC2CLASSIFIEDDIR
ACC2CLASSIFIEDDIR_JSON = "config/acc2classified_dir.json"


# dump mapping so split_read.py can use it
Path("config").mkdir(exist_ok=True)
ACC2DIR_JSON = "config/acc2dirname.json"
with open(ACC2DIR_JSON, "w") as fh:
    json.dump(ACC2DIR, fh, indent=2)

# ────────────────────────────────────────────────────────────────
# Helper function for splitting reads
def detect_out_root(fq1):
    sample = os.path.basename(os.path.dirname(fq1)) # e.g. data/sample_1/…
    return f"results/{sample}"

OUT_ROOT = detect_out_root(FQ1)


# ────────────────────────────────────────────────────────────────
# Helper function to find .pb.gz to each reference pathogen
ACC2PB = {}
for acc, fasta_path in ACC2FASTA.items():
    pdir = Path(fasta_path).parent
    # allow either *.pb.gz or *.pb – first match wins
    pb_candidates = list(itertools.chain(
        pdir.glob("*.pb.gz"),
        pdir.glob("*.pb")         
    ))
    if not pb_candidates:
        raise ValueError(f"No PB/MAT file found next to {fasta_path}")
    ACC2PB[acc] = str(pb_candidates[0])

# ────────────────────────────────────────────────────────────────
# Helper function for prepare wepp
def optional_file(path):
    """
    Return the path as a string if it exists,
    otherwise return an empty list 
    """
    return path if os.path.exists(path) else []

SAMPLE_TAG   = config["DIR"]        
def tag(acc):  
    return f"{ACC2DIR[acc]}_{SAMPLE_TAG}"


# ────────────────────────────────────────────────────────────────
# Helper function to get the path of tree file and genome file in WEPP/data
def WEPP_TREE(acc):
    """
    Return full path to the tree file for a given accession.
    Uses the output name in the form of WEPP/data/<dir_tag>/<tree_filename>.pb.gz
    """
    dir_tag = tag(acc)
    tree_filename = os.path.basename(ACC2PB[acc])  
    return f"{WEPP_DATA_DIR}/{dir_tag}/{tree_filename}"


def WEPP_REF(acc):
    """
    Return full path to the genome file for a given accession.
    Uses the output name in the form of WEPP/data/<dir_tag>/<filename>.pb.gz
    """
    dir_tag = tag(acc)
    tree_filename = os.path.basename(ACC2FASTA[acc])
                                                    
    return f"{WEPP_DATA_DIR}/{dir_tag}/{tree_filename}"

# ────────────────────────────────────────────────────────────────
ACC2CLASSIFIEDDIR = {}

def final_targets(_wc):
    ckpt = checkpoints.build_acc2classified_dir.get()  # wait for split
    files = []
    if Path(ACC2CLASSIFIEDDIR_JSON).exists():
        with open(ACC2CLASSIFIEDDIR_JSON) as f:
            ACC2CLASSIFIEDDIR = json.load(f)
    else:
        ACC2CLASSIFIEDDIR = {}          # fall-back (empty) so the Snakefile still parses
    for acc in REF_ACCESSIONS:
        dir_ = ACC2CLASSIFIEDDIR.get(acc)  # use .get() to avoid KeyError
        if dir_:
            dir_tag = f"{dir_}_{TAG}"
            files.append(f"{WEPP_RESULTS_DIR}/{dir_tag}/{acc}_run.txt")
    files.append(FQ1)
    if not IS_SINGLE_END:
        files.append(FQ2)
    return files

def final_results_files(wc):
    # ── wait until the checkpoint that builds acc2classified_dir.json is done ──
    ckpt = checkpoints.build_acc2classified_dir.get(**wc)   # ← blocks until file exists
    classified_json = ckpt.output[0]                        # same as ckpt.output.classified
    with open(classified_json) as fh:
        acc2dir = json.load(fh)                             # {acc : "NC_045512" …}

    wanted = []

    # one job per accession → per-accession FASTQs produced by split_per_accession
    for acc in REF_ACCESSIONS:
        out_dir = acc2dir.get(acc)          # directory name chosen by split_read.py
        if not out_dir:                     # accession not present → skip
            continue

        # R1 (and, if paired-end, R2) produced by split_per_accession
        wanted.append(f"{OUT_ROOT}/{out_dir}/{acc.split('.')[0]}_R1.fq.gz")
        if not IS_SINGLE_END:
            wanted.append(f"{OUT_ROOT}/{out_dir}/{acc.split('.')[0]}_R2.fq.gz")
    return wanted



# ────────────────────────────────────────────────────────────────
# Helper function for run wepp
def has_reads(fq):
    """
    Return shell code that exits 0 (true) when the FASTQ has ≥ 1 read.
    Works for .gz and plain fastq. 4 lines = 1 read.
    """
    return f'( (gzip -cd {fq} 2>/dev/null || cat {fq}) | head -4 | wc -l ) -ge 4'

def _classified_empty(p):
    """Return True only if file exists and its JSON is exactly {} (or file missing)."""
    try:
        with open(p) as fh:
            # Fast path: avoid loading huge JSON accidentally
            txt = fh.read().strip()
            return txt == "{}"
    except FileNotFoundError:
        return True


CLASSIFIED_EMPTY = _classified_empty(ACC2CLASSIFIEDDIR_JSON)

if CLASSIFIED_EMPTY:
    rule all:
        input:
            ACC2CLASSIFIEDDIR_JSON,
            f"{OUT_ROOT}/kraken_output.txt",
            f"{OUT_ROOT}/kraken_report.txt",
            f"{OUT_ROOT}/.split_done"
else: 
    rule all:
        input:
            final_targets,
            final_results_files,
            ACC2CLASSIFIEDDIR_JSON,
            f"{OUT_ROOT}/kraken_output.txt",
            f"{OUT_ROOT}/kraken_report.txt",

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

    def sim_reads_R1(wc):
        ckpt = checkpoints.split_genomes.get(**wc)           # wait for checkpoint
        genomes, = glob_wildcards(os.path.join(ckpt.output[0], "{genome}.fa"))
        return expand("simulated_reads/fastq/{genome}_R1.fq.gz", genome=genomes)

    def sim_reads_R2(wc):
        if IS_SINGLE_END:
            return []                                        # no second read file
        ckpt = checkpoints.split_genomes.get(**wc)
        genomes, = glob_wildcards(os.path.join(ckpt.output[0], "{genome}.fa"))
        return expand("simulated_reads/fastq/{genome}_R2.fq.gz", genome=genomes)


    rule merge_reads:
        input:
            sim_reads_1 = lambda wc: expand("simulated_reads/fastq/{genome}_R1.fq.gz", genome=get_genomes_from_checkpoint(wc)),
            sim_reads_2 = lambda wc: expand("simulated_reads/fastq/{genome}_R2.fq.gz", genome=get_genomes_from_checkpoint(wc))
        output:
            R1 = FQ1,
            R2 = FQ2
        params:
            FQ1_sim = "simulated_reads/fastq/merged_R1.fq.gz",
            FQ2_sim = "simulated_reads/fastq/merged_R2.fq.gz",
        shell:
            """
            zcat simulated_reads/fastq/*_R1.fq.gz | gzip > {params.FQ1_sim}
            zcat simulated_reads/fastq/*_R2.fq.gz | gzip > {params.FQ2_sim}
            
            mv {params.FQ1_sim} {output.R1}
            mv {params.FQ2_sim} {output.R2}

            rm -rf simulated_reads
            """


# ─── helper function for names ──────────────────────────────────────────
FQ1_GZ = FQ1 if str(FQ1).endswith(".gz") else str(FQ1) + ".gz"
FQ2_GZ = (
    "" if IS_SINGLE_END
    else (FQ2 if str(FQ2).endswith(".gz") else str(FQ2) + ".gz")
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

# 6.5) Kraken Visualization
rule kraken_visualization:
    input:
        report = "kraken_report.txt"
    output:
        dir = "classification_proportions.png"
    shell:
        """
        python scripts/kraken_data_visualization.py {input.report}
        """

checkpoint build_acc2classified_dir:
    input:
        kraken_out = f"kraken_output.txt",
        mapping    = TAXID_MAP,
        acc2dir    = ACC2DIR_JSON
    output:
        classified = ACC2CLASSIFIEDDIR_JSON
    params:
        script = "scripts/generate_acc2classified_dir.py"
    shell:
        "python {params.script} "
        "--kraken-out {input.kraken_out} "
        "--mapping {input.mapping} "
        "--acc2dir {input.acc2dir} "
        "-o {output.classified}"

# 7) Split Kraken output into per-taxid FASTQs 
# --- config / flags ---
ACC2CLASSIFIEDDIR_JSON = "config/acc2classified_dir.json"


# IS_SINGLE_END already defined elsewhere
# OUT_ROOT, TAXID_MAP, ACC2DIR_JSON, FQ1, FQ2, REF_ACCESSIONS defined elsewhere

# Common bits
SPLIT_SCRIPT = "scripts/split_read.py"
REF_ARG = "" if not REF_ACCESSIONS else "--ref-accessions " + ",".join(REF_ACCESSIONS)
DIR_ARG = f"--acc2dir {ACC2DIR_JSON}"

SHELL_BASE = r"""
    python {params.script} \
        --kraken-out {input.kraken_out} \
        --mapping    {input.mapping} \
        --r1         {input.r1} \
        {params.r2_arg} \
        {params.ref_arg} \
        {params.dir_arg} \
        --out-dir    {params.dir} \
        --pigz-threads {threads} \
        && touch {output.done}
"""

# ----------------- CASE 1: ACC2CLASSIFIEDDIR_JSON is EMPTY ------------------
if CLASSIFIED_EMPTY:
    if IS_SINGLE_END:
        rule split_per_accession:
            input:
                mapping    = TAXID_MAP,
                acc2dir    = ACC2DIR_JSON,
                r1         = FQ1,
                r2         = FQ2,
                kraken_out = "kraken_output.txt",
                classified = ACC2CLASSIFIEDDIR_JSON
            output:
                done = f"{OUT_ROOT}/.split_done"
            params:
                script  = "scripts/split_read.py",
                r2_arg  = f"--r2 {FQ2}",
                dir_arg = f"--acc2dir {ACC2DIR_JSON}",
                ref_arg = "" if not REF_ACCESSIONS else
                        "--ref-accessions " + ",".join(REF_ACCESSIONS),
                dir     = OUT_ROOT
            threads: 32
            shell:
                r"""
                python {params.script} \
                --kraken-out {input.kraken_out} \
                --mapping    {input.mapping} \
                --r1         {input.r1} \
                {params.ref_arg} \
                {params.dir_arg} \
                --out-dir    {params.dir} \
                --pigz-threads {threads} \
                && touch {output.done}
                """
    else:
        rule split_per_accession:
            input:
                mapping    = TAXID_MAP,
                acc2dir    = ACC2DIR_JSON,
                r1         = FQ1,
                r2         = FQ2,
                kraken_out = "kraken_output.txt",
                classified = ACC2CLASSIFIEDDIR_JSON
            output:
                done = f"{OUT_ROOT}/.split_done"
            params:
                script  = "scripts/split_read.py",
                r2_arg  = f"--r2 {FQ2}",
                dir_arg = f"--acc2dir {ACC2DIR_JSON}",
                ref_arg = "" if not REF_ACCESSIONS else
                        "--ref-accessions " + ",".join(REF_ACCESSIONS),
                dir     = OUT_ROOT
            threads: 32
            shell:
                r"""
                python {params.script} \
                --kraken-out {input.kraken_out} \
                --mapping    {input.mapping} \
                --r1         {input.r1} \
                {params.r2_arg} \
                {params.ref_arg} \
                {params.dir_arg} \
                --out-dir    {params.dir} \
                --pigz-threads {threads} \
                && touch {output.done}
                """

# ----------------- CASE 2: ACC2CLASSIFIEDDIR_JSON is NOT EMPTY --------------
else:
    if IS_SINGLE_END:
        rule split_per_accession:
            input:
                mapping    = TAXID_MAP,
                acc2dir    = ACC2DIR_JSON,
                r1         = FQ1,
                r2         = FQ2,
                kraken_out = "kraken_output.txt",
                classified = ACC2CLASSIFIEDDIR_JSON
            output:
                r1_out = f"{OUT_ROOT}/{{out_dir}}/{{acc}}_R1.fq.gz",
            params:
                script  = "scripts/split_read.py",
                r2_arg  = f"--r2 {FQ2}",
                dir_arg = f"--acc2dir {ACC2DIR_JSON}",
                ref_arg = "" if not REF_ACCESSIONS else
                        "--ref-accessions " + ",".join(REF_ACCESSIONS),
                dir     = OUT_ROOT,
                done = f"{OUT_ROOT}/.split_done"
            threads: 32
            shell:
                r"""
                python {params.script} \
                --kraken-out {input.kraken_out} \
                --mapping    {input.mapping} \
                --r1         {input.r1} \
                {params.ref_arg} \
                {params.dir_arg} \
                --out-dir    {params.dir} \
                --pigz-threads {threads}\
                && touch {params.done}            
                """
    else:
        rule split_per_accession:
            input:
                mapping    = TAXID_MAP,
                acc2dir    = ACC2DIR_JSON,
                r1         = FQ1,
                r2         = FQ2,
                kraken_out = "kraken_output.txt",
                classified = ACC2CLASSIFIEDDIR_JSON
            output:
                r1_out = f"{OUT_ROOT}/{{out_dir}}/{{acc}}_R1.fq.gz",
                r2_out = f"{OUT_ROOT}/{{out_dir}}/{{acc}}_R2.fq.gz",
            params:
                script  = "scripts/split_read.py",
                r2_arg  = f"--r2 {FQ2}",
                dir_arg = f"--acc2dir {ACC2DIR_JSON}",
                ref_arg = "" if not REF_ACCESSIONS else
                        "--ref-accessions " + ",".join(REF_ACCESSIONS),
                dir     = OUT_ROOT,
                done = f"{OUT_ROOT}/.split_done"
            threads: 32
            shell:
                r"""
                python {params.script} \
                --kraken-out {input.kraken_out} \
                --mapping    {input.mapping} \
                --r1         {input.r1} \
                {params.r2_arg} \
                {params.ref_arg} \
                {params.dir_arg} \
                --out-dir    {params.dir} \
                --pigz-threads {threads}\
                && touch {params.done}            
                """


rule move_shared_files:
    input:
        acc_reads = expand(f"{OUT_ROOT}/{{dir}}/{{acc}}_R1.fq.gz", zip,
                           dir=[ACC2CLASSIFIEDDIR[acc] for acc in REF_ACCESSIONS if acc in ACC2CLASSIFIEDDIR],
                           acc=[acc for acc in REF_ACCESSIONS if acc in ACC2CLASSIFIEDDIR]),
        kraken_out = "kraken_output.txt",
        kraken_report = "kraken_report.txt",
        viz = "classification_proportions.png",
        classified = ACC2CLASSIFIEDDIR_JSON
    output:
        new_kraken_out = f"{OUT_ROOT}/kraken_output.txt",
        new_kraken_report = f"{OUT_ROOT}/kraken_report.txt",
        new_classified = f"{OUT_ROOT}/classification_proportions.png"
    shell:
        r"""
        mv {input.kraken_out} {output.new_kraken_out}
        mv {input.kraken_report} {output.new_kraken_report}
        mv {input.viz} {output.new_classified}
        """


# 7.5) Prepare wepp input directory
# For every accession (wildcard: {acc}) copy the reference *.fa and the matching 
# *.pb.gz / *.pb into the same OUT_ROOT/<acc>/ directory that already contains 
# the split reads.
# helper defined
TAG = config["DIR"]
def tag(acc):
    return f"{ACC2DIR[acc]}_{TAG}"

if IS_SINGLE_END:
    # ───────── single-end ────────────────────────────────────────────────────
    rule prepare_wepp_inputs:
        input:
            r1          = lambda wc:
                f"{OUT_ROOT}/{wc.dir_tag.replace('_' + TAG, '')}/"
                f"{wc.acc.split('.')[0]}_R1.fq.gz",
            fasta       = lambda wc: ACC2FASTA[wc.acc],
            pb          = lambda wc: ACC2PB[wc.acc],
            viz         = f"{OUT_ROOT}/classification_proportions.png", 
            classified  = ACC2CLASSIFIEDDIR_JSON,
            kraken_out  = f"{OUT_ROOT}/kraken_output.txt", 
            kraken_report = f"{OUT_ROOT}/kraken_report.txt",
        output:
            new_r1      = f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}.fastq.gz",
        params:
            data_dir    = lambda wc: f"{WEPP_DATA_DIR}/{wc.dir_tag}",
            tree_dest   = lambda wc: WEPP_TREE(wc.acc),
            fasta_dest  = lambda wc: WEPP_REF(wc.acc),
        shell:
            r"""
            mkdir -p {params.data_dir}
            cp {input.fasta} {params.fasta_dest}
            cp {input.pb}    {params.tree_dest}

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
            r1          = lambda wc:
                f"{OUT_ROOT}/{wc.dir_tag.replace('_' + TAG, '')}/"
                f"{wc.acc.split('.')[0]}_R1.fq.gz",
            r2          = lambda wc:
                f"{OUT_ROOT}/{wc.dir_tag.replace('_' + TAG, '')}/"
                f"{wc.acc.split('.')[0]}_R2.fq.gz",
            fasta       = lambda wc: ACC2FASTA[wc.acc],
            pb          = lambda wc: ACC2PB[wc.acc],
            viz         = f"{OUT_ROOT}/classification_proportions.png",   # ← f-string
            classified  = ACC2CLASSIFIEDDIR_JSON,
            kraken_out  = f"{OUT_ROOT}/kraken_output.txt",
            kraken_report = f"{OUT_ROOT}/kraken_report.txt",
        output:
            new_r1      = f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}_R1.fastq.gz",
            new_r2      = f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}_R2.fastq.gz",
        params:
            data_dir    = lambda wc: f"{WEPP_DATA_DIR}/{wc.dir_tag}",
            tree_dest   = lambda wc: WEPP_TREE(wc.acc),
            fasta_dest  = lambda wc: WEPP_REF(wc.acc),
        shell:
            r"""
            mkdir -p {params.data_dir}
            cp {input.fasta} {params.fasta_dest}
            cp {input.pb}    {params.tree_dest}

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
        # FASTQs are stored in WEPP/data/<acc>_<TAG>/
        r1         = f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}.fastq.gz"
                      if IS_SINGLE_END else
                      f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}_R1.fastq.gz",
        r2         = [] if IS_SINGLE_END else
                      f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}_R2.fastq.gz",
    output:
        run_txt    = f"{WEPP_RESULTS_DIR}/{{dir_tag}}/{{acc}}_run.txt",
    threads: config.get("wepp_threads", 32)
    params:
        snakefile   = str(WEPP_WORKFLOW),
        workdir     = str(WEPP_ROOT),
        seq_type    = config.get("SEQUENCING_TYPE", "d").lower(),
        primer_bed  = config["PRIMER_BED"],
        min_af      = config["MIN_AF"],
        min_q       = config["MIN_Q"],
        max_reads   = config["MAX_READS"],
        clade_list  = config["CLADE_LIST"],
        clade_idx   = config["CLADE_IDX"],
        cfgfile     = str(WEPP_CONFIG),
        resultsdir  = str(WEPP_RESULTS_DIR),
        prefix      = lambda wc: wc.acc.split('.')[0],
        ref_name    = lambda wc: f"{wc.acc}.fa",
        tag_dir     = lambda wc: tag(wc.acc),  
        tree_name   = lambda wc: os.path.basename(WEPP_TREE(wc.acc)),
        tree_full   = lambda wc: WEPP_TREE(wc.acc),
        fasta_name  = lambda wc: os.path.basename(WEPP_REF(wc.acc)),
        fasta_full  = lambda wc: WEPP_REF(wc.acc),
        customconfig = lambda wc: ACC2CONFIG.get(wc.acc, "")
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
            echo "No reads for {wildcards.acc}; removing data folder and skipping WEPP."
            rm -rf "{WEPP_DATA_DIR}/{params.tag_dir}"
            mkdir -p {params.resultsdir}/{wildcards.acc}
            touch {output.run_txt}
            exit 0
        fi

        if [ -n "{params.customconfig}" ] && [ -f "{params.customconfig}" ]; then
            custom_arg="--customconfig {params.customconfig}"
        else
            custom_arg=""
        fi

        python ./scripts/run_inner.py \
            --dir        {params.tag_dir} \
            --prefix     {params.prefix} \
            --primer_bed {params.primer_bed} \
            --min_af     {params.min_af} \
            --min_q      {params.min_q} \
            --max_reads  {params.max_reads} \
            --tree       {params.tree_name} \
            --ref        {params.fasta_name} \
            --clade_idx  {params.clade_idx} \
            --snakefile  {params.snakefile} \
            --workdir    {params.workdir} \
            --clade_list {params.clade_list} \
            $custom_arg \
            --cores      {threads} \
            --sequencing_type {params.seq_type}

        touch {output.run_txt}
        """
