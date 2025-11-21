# Snakefile

# imports
import os
import csv, re, itertools
import subprocess

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
        "Call snakemake with, e.g.,  --config DIR=TEST_DIR"
    )

rule all:
    input:
        f"results/{DIR}/classification_proportions.png",
        f"results/{DIR}/kraken_output.txt"

KRAKEN_DB = config.get("KRAKEN_DB")
if KRAKEN_DB is None:
    raise ValueError(
        "Please provide the path to the Kraken2 database via\n"
        "  --config KRAKEN_DB=<folder path>"
    )


WEPP_ROOT        = "WEPP"
WEPP_WORKFLOW    = "WEPP/workflow/Snakefile"
WEPP_RESULTS_DIR = "WEPP/results"
WEPP_DATA_DIR    = "WEPP/data"
WEPP_CMD_LOG     = "WEPP/cmd_log"

CONFIG      = "config/config.yaml"

# Dashboard variable from config
DASHBOARD_ENABLED = config.get("DASHBOARD_ENABLED", False)
MIN_DEPTH = float(config.get("MIN_DEPTH_FOR_WEPP", 0))

## Get comma-separated pathogens string
#pathogen_str = config.get("PATHOGENS_FOR_DASHBOARD", "")
#PATHOGENS_FOR_DASHBOARD = [p.strip() for p in pathogen_str.split(",") if p.strip()]
#print(f"Dashboard = {DASHBOARD_ENABLED}", file=sys.stderr)
#print("PATHOGENS_FOR_DASHBOARD:", ", ".join(PATHOGENS_FOR_DASHBOARD))

# File that stores the mapping from accession to taxid
NCBI_TAXID_MAP = os.path.join(config["KRAKEN_DB"], "seqid2taxid.map")

## New map output path
#TAXID_MAP = os.path.join("config", "accession2taxid.map")
## Convert the map
#with open(NCBI_TAXID_MAP) as infile, open(TAXID_MAP, "w") as outfile:
#    for line in infile:
#        if line.strip():
#            full_id, taxid = line.strip().split()
#            accession = full_id.split("|")[-1]
#            outfile.write(f"{accession}\t{taxid}\n")

IS_SINGLE_END = config.get("SEQUENCING_TYPE", "d").lower() in {"s", "n"}

fq_dir = Path("data") / config["DIR"]
if not fq_dir.exists():
    raise FileNotFoundError(f"Input folder {DIR} does not exist")

# look for input FASTQs
def gzip_if_needed(file_path):
    """Gzips a file if it's not already gzipped, returns new path."""
    if not str(file_path).endswith(".gz"):
        print(f"Compressing {file_path}...", file=sys.stderr)
        try:
            subprocess.run(["gzip", str(file_path)], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error compressing {file_path}: {e}", file=sys.stderr)
            raise
        return Path(str(file_path) + ".gz")
    return file_path

if IS_SINGLE_END:
    r1_files = sorted(fq_dir.glob("*.fastq*"))
    if len(r1_files) != 1:
        raise RuntimeError(
            f"{fq_dir} must contain exactly one .fastq file for single-end mode, "
            f"but found {len(r1_files)}."
        )
    FQ1 = str(gzip_if_needed(r1_files[0]))
    FQ2 = ""
else: # Paired-end
    r1_files = sorted(fq_dir.glob("*_R1.fastq*"))
    if len(r1_files) != 1:
        raise RuntimeError(
            f"{fq_dir} must contain exactly one *_R1.fastq* file for paired-end mode, "
            f"but found {len(r1_files)}."
        )

    r2_files = sorted(fq_dir.glob("*_R2.fastq*"))
    if len(r2_files) != 1:
        raise RuntimeError(
            f"{fq_dir} must contain exactly one *_R2.fastq* file for paired-end mode, "
            f"but found {len(r2_files)}."
        )
    FQ1 = str(gzip_if_needed(r1_files[0]))
    FQ2 = str(gzip_if_needed(r2_files[0]))

print(f"Input FASTQs chosen:\n  FQ1 = {FQ1}\n  FQ2 = {FQ2}", file=sys.stderr)


# ────────────────────────────────────────────────────────────────
# Adds new pathogens with interative script to data/pathogens_for_wepp
# ────────────────────────────────────────────────────────────────

PATHOGEN_ROOT = Path("data/pathogens_for_wepp")
haplotype_pathogens = []

rule add_pathogen:
    output:
        f"results/{DIR}/.add_pathogen.done"
    params:
        db=KRAKEN_DB
    threads: workflow.cores
    run:
        subprocess.run(
            ["python", "scripts/add_ref_mat.py", "--db", params.db, "--threads", str(threads)],
            check=True
        )

        if PATHOGEN_ROOT.exists() and PATHOGEN_ROOT.is_dir():
            for pathogen_dir in filter(Path.is_dir, PATHOGEN_ROOT.iterdir()):
                fasta_files = list(pathogen_dir.glob("*.[fF][aAn]*"))
                pb_files = list(pathogen_dir.glob("*.pb*"))

                num_fasta = len(fasta_files)
                num_pb = len(pb_files)

                if num_fasta == 1 and num_pb == 1:
                    haplotype_pathogens.append(pathogen_dir)
                    continue

                if num_fasta != 1:
                    print(f"\n[WARN]: Expected 1 FASTA file in {pathogen_dir}, but found {num_fasta}.\n", file=sys.stderr)
                
                if num_pb != 1:
                    print(f"\n[WARN]: Expected 1 .pb or .pb.gz file in {pathogen_dir}, but found {num_pb}.\n", file=sys.stderr)
        
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        Path(output[0]).touch()

rule kraken:
    input:
        r1 = FQ1,
        r2 = lambda wc: [] if IS_SINGLE_END else [FQ2],
        add_pathogen_done = rules.add_pathogen.output,
    output:
        kraken_report = f"results/{DIR}/kraken_report.txt",
        kraken_out    = f"results/{DIR}/kraken_output.txt",
    threads: workflow.cores
    params:
        db        = KRAKEN_DB,
        mode_flag = "" if IS_SINGLE_END else "--paired",
    shell:
        r"""
        kraken2 --db {params.db} --threads {threads} {params.mode_flag} \
                {input.r1} {input.r2} \
                --report {output.kraken_report} \
                --output {output.kraken_out}
        """

rule kraken_visualization:
    input:
        kraken_report = rules.kraken.output.kraken_report,
    output:
        png = f"results/{DIR}/classification_proportions.png"
    shell:
        r"""
        python scripts/kraken_data_visualization.py \
            {input.kraken_report} {output.png} {float(config.get("MIN_PROP_FOR_WEPP", 0.01))}
        """
#
#rule split_read:
#    input:
#        kraken_out = rules.kraken.output.kraken_out,
#        kraken_report = rules.kraken.output.kraken_report,
#        png = rules.kraken_visualization.output.png,
#        r1 = FQ1,
#        r2 = lambda wc: [] if IS_SINGLE_END else [FQ2],
#    output:
#        "results/{DIR}/.split_read.done"
#    params:
#        script = "scripts/split_read.py",
#    threads: workflow.cores
#    run:
#        subprocess.run(
#            ["python", params.script, "-k", input.kraken_out, "-r", input.kraken_out, "--r1", input.r1, "--r2", input.r2, "-o results/{DIR}", "-t", threads],
#            check=True
#        )
#        Path(output[0]).touch()


## 6.5) Kraken Visualization
#
#checkpoint build_acc2classified_dir:
#    input:
#        kraken_out = f"{OUT_ROOT}/kraken_output.txt",
#        mapping    = TAXID_MAP,
#        acc2dir    = ACC2DIR_JSON
#    output:
#        classified = ACC2CLASSIFIEDDIR_JSON
#    params:
#        script = "scripts/generate_acc2classified_dir.py"
#    shell:
#        "python {params.script} "
#        "--kraken-out {input.kraken_out} "
#        "--mapping {input.mapping} "
#        "--acc2dir {input.acc2dir} "
#        "-o {output.classified}"
#
## 7) Split Kraken output into per-taxid FASTQs 
## --- config ---
#ACC2CLASSIFIEDDIR_JSON = "config/acc2classified_dir.json"
#
## canonical paths used for split_per_accession rule
#KRAKEN_OUT     = f"{OUT_ROOT}/kraken_output.txt"
#KRAKEN_REPORT  = f"{OUT_ROOT}/kraken_report.txt"
#SPLIT_SENTINEL = f"{OUT_ROOT}/.split_done"
#
#SPLIT_INPUTS = {
#    "mapping":       TAXID_MAP,
#    "acc2dir":       ACC2DIR_JSON,
#    "r1":            FQ1,
#    "kraken_out":    KRAKEN_OUT,
#    "kraken_report": KRAKEN_REPORT,
#    "classified":    ACC2CLASSIFIEDDIR_JSON,
#}
#if not IS_SINGLE_END:
#    SPLIT_INPUTS["r2"] = FQ2
#
#SPLIT_PARAMS = {
#    "script":    "scripts/split_read.py",
#    "dir_arg":   f"--acc2dir {ACC2DIR_JSON}",
#    "ref_arg":   "" if not REF_ACCESSIONS else "--ref-accessions " + ",".join(REF_ACCESSIONS),
#    "dir":       OUT_ROOT
#}
#
#r2_arg = "--r2 {input.r2} \\" if not IS_SINGLE_END else ""
#
#SPLIT_SHELL = r"""
#python {params.script} \
#--kraken-out {input.kraken_out} \
#--kraken-report {input.kraken_report} \
#--mapping {input.mapping} \
#--r1 {input.r1} \
#""" + r2_arg + r"""
#{params.ref_arg} \
#{params.dir_arg} \
#--out-dir {params.dir} \
#--pigz-threads {threads} \
#"""
#
## ----------------- CASE 1: ACC2CLASSIFIEDDIR_JSON is EMPTY ------------------
#if CLASSIFIED_EMPTY:
#    rule split_per_accession:
#        input:  **SPLIT_INPUTS
#        output:
#            done = f"{OUT_ROOT}/.split_done"
#        params: **SPLIT_PARAMS,
#        threads: workflow.cores
#        shell:
#            SPLIT_SHELL + r" && touch {output.done}"
#            
## ----------------- CASE 2: ACC2CLASSIFIEDDIR_JSON is NOT EMPTY --------------
#else:
#    if IS_SINGLE_END:
#        rule split_per_accession:
#            input: **SPLIT_INPUTS
#            output:
#                r1_out = f"{OUT_ROOT}/{{out_dir}}/{{acc}}_R1.fq.gz",
#            params: 
#                **SPLIT_PARAMS,
#                done = SPLIT_SENTINEL
#            threads: workflow.cores
#            shell:
#                SPLIT_SHELL + r" && touch {params.done}"
#    else:
#        rule split_per_accession:
#            input: **SPLIT_INPUTS
#            output:
#                r1_out = f"{OUT_ROOT}/{{out_dir}}/{{acc}}_R1.fq.gz",
#                r2_out = f"{OUT_ROOT}/{{out_dir}}/{{acc}}_R2.fq.gz",
#            params: 
#                **SPLIT_PARAMS,
#                done = SPLIT_SENTINEL
#            threads: workflow.cores
#            shell:
#                SPLIT_SHELL + r" --r2 {input.r2} && touch {params.done}"
#
#checkpoint coverage_calculate:
#    input:
#        classified = ACC2CLASSIFIEDDIR_JSON, 
#        r1s        = split_fastqs_for_coverage 
#    output:
#        coverage = ACC2COVERED_JSON
#    params:
#        out_root      = OUT_ROOT,
#        pathogen_root = str(PATHOGEN_ROOT),
#        min_depth     = MIN_DEPTH,
#        script        = "scripts/calc_coverage_json.py"
#    shell:
#        r"""
#        python {params.script} \
#          --acc2classified {input.classified} \
#          --out-root {params.out_root} \
#          --pathogen-root {params.pathogen_root} \
#          --min-depth {params.min_depth} \
#          --out-json {output.coverage}
#        """
#
#
#
## 7.5) Prepare wepp input directory
## For every accession (wildcard: {acc}) copy the reference *.fa and the matching 
## *.pb.gz / *.pb into the same OUT_ROOT/<acc>/ directory that already contains 
## the split reads.
## helper defined
#TAG = config["DIR"]
#def tag(acc):
#    return f"{ACC2DIR[acc]}_{TAG}"
#
#if IS_SINGLE_END:
#    # ───────── single-end ────────────────────────────────────────────────────
#    rule prepare_wepp_inputs:
#        input:
#            r1          = lambda wc:
#                f"{OUT_ROOT}/{wc.dir_tag.replace('_' + TAG, '')}/"
#                f"{wc.acc.split('.')[0]}_R1.fq.gz",
#            fasta       = lambda wc: ACC2FASTA[wc.acc],
#            pb          = lambda wc: ACC2PB[wc.acc],
#            viz         = f"{OUT_ROOT}/classification_proportions.png", 
#            classified  = ACC2CLASSIFIEDDIR_JSON,
#            kraken_out  = f"{OUT_ROOT}/kraken_output.txt", 
#            kraken_report = f"{OUT_ROOT}/kraken_report.txt",
#            coverage = ACC2COVERED_JSON
#        output:
#            new_r1      = f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}.fastq.gz",
#        params:
#            data_dir    = lambda wc: f"{WEPP_DATA_DIR}/{wc.dir_tag}",
#            tree_dest   = lambda wc: WEPP_TREE(wc.acc),
#            fasta_dest  = lambda wc: WEPP_REF(wc.acc),
#            jsonl_src   = lambda wc: ACC2JSONL.get(wc.acc, ""),
#            jsonl_dest  = lambda wc: WEPP_JSONL(wc.acc)
#        shell:
#            r"""
#            mkdir -p {params.data_dir}
#            cp {input.fasta} {params.fasta_dest}
#            cp {input.pb}    {params.tree_dest}
#
#            if [ -f "{input.r1}" ]; then
#                cp {input.r1} {output.new_r1}
#            else
#                echo | gzip -c > {output.new_r1}
#            fi
#
#            if [ -n "{params.jsonl_src}" ]; then
#                cp "{params.jsonl_src}" "{params.jsonl_dest}"
#            fi
#            """
#
#else:
#    # ───────── paired-end ───────────────────────────────────────────────────
#    rule prepare_wepp_inputs:
#        input:
#            r1          = lambda wc:
#                f"{OUT_ROOT}/{wc.dir_tag.replace('_' + TAG, '')}/"
#                f"{wc.acc.split('.')[0]}_R1.fq.gz",
#            r2          = lambda wc:
#                f"{OUT_ROOT}/{wc.dir_tag.replace('_' + TAG, '')}/"
#                f"{wc.acc.split('.')[0]}_R2.fq.gz",
#            fasta       = lambda wc: ACC2FASTA[wc.acc],
#            pb          = lambda wc: ACC2PB[wc.acc],
#            viz         = f"{OUT_ROOT}/classification_proportions.png",
#            classified  = ACC2CLASSIFIEDDIR_JSON,
#            kraken_out  = f"{OUT_ROOT}/kraken_output.txt",
#            kraken_report = f"{OUT_ROOT}/kraken_report.txt",
#            coverage = ACC2COVERED_JSON
#        output:
#            new_r1      = f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}_R1.fastq.gz",
#            new_r2      = f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}_R2.fastq.gz",
#        params:
#            data_dir    = lambda wc: f"{WEPP_DATA_DIR}/{wc.dir_tag}",
#            tree_dest   = lambda wc: WEPP_TREE(wc.acc),
#            fasta_dest  = lambda wc: WEPP_REF(wc.acc),
#            jsonl_src   = lambda wc: ACC2JSONL.get(wc.acc, ""),
#            jsonl_dest  = lambda wc: WEPP_JSONL(wc.acc)
#        shell:
#            r"""
#            mkdir -p {params.data_dir}
#            cp {input.fasta} {params.fasta_dest}
#            cp {input.pb}    {params.tree_dest}
#
#            if [ -f "{input.r1}" ]; then
#                cp {input.r1} {output.new_r1}
#            else
#                echo | gzip -c > {output.new_r1}
#            fi
#
#            if [ -f "{input.r2}" ]; then
#                cp {input.r2} {output.new_r2}
#            else
#                echo | gzip -c > {output.new_r2}
#            fi
#
#            if [ -n "{params.jsonl_src}" ]; then
#                cp "{params.jsonl_src}" "{params.jsonl_dest}"
#            fi
#            """
#
## 8) Invoke WEPP’s Snakefile for each taxid
## Run WEPP
#rule run_wepp:
#    input:
#        # FASTQs are stored in WEPP/data/<acc>_<TAG>/
#        r1         = f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}.fastq.gz"
#                      if IS_SINGLE_END else
#                      f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}_R1.fastq.gz",
#        r2         = [] if IS_SINGLE_END else
#                      f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}_R2.fastq.gz",
#    output:
#        run_txt    = f"{WEPP_RESULTS_DIR}/{{dir_tag}}/{{acc}}_run.txt",
#    threads: workflow.cores
#    params:
#        snakefile   = str(WEPP_WORKFLOW),
#        workdir     = str(WEPP_ROOT),
#        seq_type    = config.get("SEQUENCING_TYPE", "d").lower(),
#        primer_bed  = config["PRIMER_BED"],
#        min_af      = config["MIN_AF"],
#        min_q       = config["MIN_Q"],
#        max_reads   = config["MAX_READS"],
#        clade_list  = config["CLADE_LIST"],
#        clade_idx   = config["CLADE_IDX"],
#        min_prop    = config["MIN_PROP"],
#        min_len    = config["MIN_LEN"],
#        cfgfile     = str(CONFIG),
#        resultsdir  = str(WEPP_RESULTS_DIR),
#        prefix      = lambda wc: wc.acc.split('.')[0],
#        ref_name    = lambda wc: f"{wc.acc}.fa",
#        tag_dir     = lambda wc: tag(wc.acc),  
#        tree_name   = lambda wc: os.path.basename(WEPP_TREE(wc.acc)),
#        tree_full   = lambda wc: WEPP_TREE(wc.acc),
#        fasta_name  = lambda wc: os.path.basename(WEPP_REF(wc.acc)),
#        fasta_full  = lambda wc: WEPP_REF(wc.acc),
#        cmd_log     = f"{WEPP_CMD_LOG}/{TAG}_dashboard_run.txt",
#        pathogens_name = lambda wc: dir(wc.acc)
#    conda:
#        "env/wepp.yaml"
#    shell:
#        r"""
#        # ── skip WEPP when both FASTQs are empty ──────────────────────
#        min_reads=100
#        has_reads() {{
#            local f=$1
#            local num_reads
#            num_reads=$( (gzip -cd "$f" 2>/dev/null || cat "$f") | wc -l )
#            [ "$((num_reads / 4))" -ge "$min_reads" ]
#        }}
#
#        if ! has_reads "{input.r1}" && ( [ -z "{input.r2}" ] || ! has_reads "{input.r2}" ); then
#            echo "Fewer than $min_reads for {wildcards.acc}; removing data folder and skipping WEPP."
#            rm -rf "{WEPP_DATA_DIR}/{params.tag_dir}"
#            mkdir -p {params.resultsdir}/{wildcards.acc}
#            touch {output.run_txt}
#            exit 0
#        fi
#
#        mkdir -p {WEPP_CMD_LOG}
#        [ -f {params.cmd_log} ] || touch {params.cmd_log}
#
#        python ./scripts/run_inner.py \
#            --dir        {params.tag_dir} \
#            --prefix     {params.prefix} \
#            --primer_bed {params.primer_bed} \
#            --min_af     {params.min_af} \
#            --min_q      {params.min_q} \
#            --max_reads  {params.max_reads} \
#            --tree       {params.tree_name} \
#            --ref        {params.fasta_name} \
#            --snakefile  {params.snakefile} \
#            --workdir    {params.workdir} \
#            --cfgfile    {params.cfgfile} \
#            --cores      {threads} \
#            --sequencing_type {params.seq_type} \
#            --pathogens_name {params.pathogens_name} \
#            --cmd_log {params.cmd_log}
#
#        touch {output.run_txt} 
#        """
#
#rule emit_dashboard_instructions:
#    input:
#        runs    = run_txts_for_all
#    output:
#        howto   = f"{WEPP_CMD_LOG}/{TAG}_dashboard_howto.txt"
#    params:
#        script = "scripts/emit_dashboard_howto.py",
#        cmd_log = f"{WEPP_CMD_LOG}/{TAG}_dashboard_run.txt"
#    shell:
#        "python scripts/emit_dashboard_howto.py --cmd-log {params.cmd_log} --out {output.howto}"
#
#
#rule run_wepp_dashboard:
#    input:
#        run_txt = f"{WEPP_RESULTS_DIR}/{{dir_tag}}/{{acc}}_run.txt",
#    output:
#        dash_txt = f"{WEPP_RESULTS_DIR}/{{dir_tag}}/{{acc}}_dashboard_run.txt",
#    threads: workflow.cores
#    params:
#        snakefile    = str(WEPP_WORKFLOW),
#        workdir      = str(WEPP_ROOT),
#        seq_type     = config.get("SEQUENCING_TYPE", "d").lower(),
#        primer_bed   = config["PRIMER_BED"],
#        min_af       = config["MIN_AF"],
#        min_q        = config["MIN_Q"],
#        max_reads    = config["MAX_READS"],
#        clade_list   = config["CLADE_LIST"],
#        clade_idx    = config["CLADE_IDX"],
#        cfgfile      = str(CONFIG),
#        resultsdir   = str(WEPP_RESULTS_DIR),
#        prefix       = lambda wc: wc.acc.split('.')[0],
#        ref_name     = lambda wc: f"{wc.acc}.fa",
#        tag_dir      = lambda wc: tag(wc.acc),
#        tree_name    = lambda wc: os.path.basename(WEPP_TREE(wc.acc)),
#        tree_full    = lambda wc: WEPP_TREE(wc.acc),
#        fasta_name   = lambda wc: os.path.basename(WEPP_REF(wc.acc)),
#        fasta_full   = lambda wc: WEPP_REF(wc.acc),
#        cmd_log      = f"{WEPP_CMD_LOG}/{TAG}_dashboard_run.txt",
#        pathogens_name = lambda wc: dir(wc.acc),
#        taxonium_file  = lambda wc: (WEPP_JSONL(wc.acc) if ACC2JSONL.get(wc.acc) else "")
#    conda:
#        "env/wepp.yaml"
#    shell:
#        r"""
#        if [ -n "{params.taxonium_file}" ]; then
#            tax_arg="--taxonium_file {params.taxonium_file}"
#        else
#            tax_arg=""
#        fi
#
#        python ./scripts/run_inner.py \
#            --dir        {params.tag_dir} \
#            --prefix     {params.prefix} \
#            --primer_bed {params.primer_bed} \
#            --min_af     {params.min_af} \
#            --min_q      {params.min_q} \
#            --max_reads  {params.max_reads} \
#            --tree       {params.tree_name} \
#            --ref        {params.fasta_name} \
#            --snakefile  {params.snakefile} \
#            --workdir    {params.workdir} \
#            --cfgfile    {params.cfgfile} \
#            --cores      {threads} \
#            --sequencing_type {params.seq_type} \
#            --pathogens_name {params.pathogens_name} \
#            --cmd_log {params.cmd_log}
#            $tax_arg
#
#        touch {output.dash_txt}
#        """
#
#
#