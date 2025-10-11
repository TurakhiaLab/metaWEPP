ctx = workflow.globals["ctx"]

IS_SINGLE_END = ctx.is_single_end
FQ1 = str(ctx.fq1)
FQ2 = [] if IS_SINGLE_END or ctx.fq2 is None else str(ctx.fq2)
OUT_ROOT = str(ctx.out_root)
KRAKEN_DB = str(ctx.kraken_db)
TAXID_MAP = str(ctx.taxid_map_path)
ACC2DIR_JSON = str(ctx.acc2dir_json_path)
ACC2CLASSIFIEDDIR_JSON = str(ctx.acc2classifieddir_json_path)
ACC2COVERED_JSON = str(ctx.acc2covered_json_path)
ADD_REF_SENTINEL = str(ctx.add_ref_sentinel)
SPLIT_SENTINEL = str(ctx.split_sentinel)
KRAKEN_OUT = str(ctx.kraken_out)
KRAKEN_REPORT = str(ctx.kraken_report)
VISUALIZATION = str(ctx.visualization)
EXCLUDE_TAXIDS_STR = ctx.exclude_taxids_str
PATHOGEN_ROOT = str(ctx.pathogen_root)
MIN_DEPTH = ctx.min_depth
SPECIES_REPORT = f"{OUT_ROOT}/species_over_1pct.txt"
REBUILD_SENTINEL = f"{OUT_ROOT}/.kraken_db_rebuilt"
REBUILD_REQUIRED = f"{OUT_ROOT}/.kraken_db_rebuild_required"
ADD_REF_INITIAL_SENTINEL = str(ctx.add_ref_initial_sentinel)

rule add_ref_mat_initial:
    output:
        done = temp(ADD_REF_INITIAL_SENTINEL)
    params:
        script = "scripts/add_ref_mat.py",
        db = KRAKEN_DB
    resources:
        serial = 1
    shell:
        r"""
        python {params.script} --db {params.db} --rebuild-sentinel {REBUILD_REQUIRED}
        touch {output.done}
        """

# 6) Kraken2 classification
rule kraken:
    input:
        r1 = FQ1,
        r2 = FQ2,
        prep = ADD_REF_INITIAL_SENTINEL,
    output:
        report     = f"{OUT_ROOT}/kraken_report.txt",
        kraken_out = f"{OUT_ROOT}/kraken_output.txt",
    threads: workflow.cores
    params:
        db        = KRAKEN_DB,
        mode_flag = "" if IS_SINGLE_END else "--paired",
        mate2     = "" if IS_SINGLE_END else str(ctx.fq2),
    resources:
        serial = 1
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
        report = f"{OUT_ROOT}/kraken_report.txt"
    output:
        png = f"{OUT_ROOT}/classification_proportions.png"
    params:
        exclude = EXCLUDE_TAXIDS_STR
    resources:
        serial = 1
    shell:
        r"""
        python scripts/kraken_data_visualization.py \
            {input.report} {output.png} "{params.exclude}"
        """

rule display_species:
    input:
        report = f"{OUT_ROOT}/kraken_report.txt"
    output:
        summary = SPECIES_REPORT
    params:
        script = "scripts/display_species.py",
        threshold = 1.0
    resources:
        serial = 1
    shell:
        "python {params.script} --report {input.report} --out {output.summary} --threshold {params.threshold}"

checkpoint build_acc2classified_dir:
    input:
        kraken_out = f"{OUT_ROOT}/kraken_output.txt",
        mapping    = TAXID_MAP,
        acc2dir    = ACC2DIR_JSON,
        ref_mat_ready = ADD_REF_SENTINEL
    output:
        classified = ACC2CLASSIFIEDDIR_JSON
    params:
        script = "scripts/generate_acc2classified_dir.py"
    resources:
        serial = 1
    shell:
        "python {params.script} "
        "--kraken-out {input.kraken_out} "
        "--mapping {input.mapping} "
        "--acc2dir {input.acc2dir} "
        "-o {output.classified}"

rule add_ref_mat:
    input:
        species = SPECIES_REPORT
    output:
        done = temp(ADD_REF_SENTINEL)
    params:
        script = "scripts/add_ref_mat.py",
        db     = KRAKEN_DB
    resources:
        serial = 1
    shell:
        r"""
        python {params.script} --db {params.db} --species-summary {input.species} --rebuild-sentinel {REBUILD_REQUIRED}
        touch {output.done}
        """

SPLIT_INPUTS = {
    "mapping":       TAXID_MAP,
    "acc2dir":       ACC2DIR_JSON,
    "r1":            FQ1,
    "kraken_out":    KRAKEN_OUT,
    "kraken_report": KRAKEN_REPORT,
    "classified":    ACC2CLASSIFIEDDIR_JSON,
    "ref_mat_ready": ADD_REF_SENTINEL,
    "db_rebuilt":    REBUILD_SENTINEL,
}
if not IS_SINGLE_END:
    SPLIT_INPUTS["r2"] = str(ctx.fq2)

SPLIT_PARAMS = {
    "script": "scripts/split_read.py",
    "dir_arg": f"--acc2dir {ACC2DIR_JSON}",
    "ref_arg": lambda wildcards: ctx.build_ref_arg(),
    "dir": OUT_ROOT,
}

_r2_segment = "" if IS_SINGLE_END else "--r2 {input.r2} "

SPLIT_SHELL_BASE = (
    "python {params.script} "
    "--kraken-out {input.kraken_out} "
    "--kraken-report {input.kraken_report} "
    "--mapping {input.mapping} "
    "--r1 {input.r1} "
    + _r2_segment +
    "{params.dir_arg} "
    "--out-dir {params.dir} "
    "--pigz-threads {threads}"
    "{params.ref_arg}"
)


rule rebuild_kraken_db:
    input:
        ref_mat_ready = ADD_REF_SENTINEL
    output:
        stamp = REBUILD_SENTINEL
    params:
        db = KRAKEN_DB,
        flag = REBUILD_REQUIRED
    threads: workflow.cores
    resources:
        serial = 1
    shell:
        r"""
        if [ -f {params.flag} ]; then
            cd {params.db} && k2 build --db . --threads {threads} && rm -f {params.flag}
        fi
        touch {output.stamp}
        """

if ctx.classified_json_empty():
    rule split_per_accession:
        input:  **SPLIT_INPUTS
        output:
            done = f"{OUT_ROOT}/.split_done"
        params: **SPLIT_PARAMS
        threads: workflow.cores
        resources:
            serial = 1
        shell:
            SPLIT_SHELL_BASE + r" && touch {output.done}"
else:
    if IS_SINGLE_END:
        rule split_per_accession:
            input: **SPLIT_INPUTS
            output:
                r1_out = f"{OUT_ROOT}/{{out_dir}}/{{acc}}_R1.fq.gz",
            params:
                **SPLIT_PARAMS,
                done = SPLIT_SENTINEL
            threads: workflow.cores
            resources:
                serial = 1
            shell:
                SPLIT_SHELL_BASE + r" && touch {params.done}"
    else:
        rule split_per_accession:
            input: **SPLIT_INPUTS
            output:
                r1_out = f"{OUT_ROOT}/{{out_dir}}/{{acc}}_R1.fq.gz",
                r2_out = f"{OUT_ROOT}/{{out_dir}}/{{acc}}_R2.fq.gz",
            params:
                **SPLIT_PARAMS,
                done = SPLIT_SENTINEL
            threads: workflow.cores
            resources:
                serial = 1
            shell:
                SPLIT_SHELL_BASE + r" && touch {params.done}"

checkpoint coverage_calculate:
    input:
        classified = ACC2CLASSIFIEDDIR_JSON,
        r1s        = split_fastqs_for_coverage
    output:
        coverage = ACC2COVERED_JSON
    params:
        out_root      = OUT_ROOT,
        pathogen_root = PATHOGEN_ROOT,
        min_depth     = MIN_DEPTH,
        script        = "scripts/calc_coverage_json.py"
    resources:
        serial = 1
    shell:
        r"""
        python {params.script} \
          --acc2classified {input.classified} \
          --out-root {params.out_root} \
          --pathogen-root {params.pathogen_root} \
          --min-depth {params.min_depth} \
          --out-json {output.coverage}
        """
