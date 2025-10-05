import os

ctx = workflow.globals["ctx"]

IS_SINGLE_END = ctx.is_single_end
OUT_ROOT = str(ctx.out_root)
WEPP_DATA_DIR = str(ctx.wepp_data_dir)
WEPP_RESULTS_DIR = str(ctx.wepp_results_dir)
WEPP_CMD_LOG = str(ctx.wepp_cmd_log)
WEPP_WORKFLOW = str(ctx.wepp_workflow)
WEPP_ROOT = str(ctx.wepp_root)
WEPP_CONFIG = str(ctx.wepp_config)
TAG = ctx.sample_tag
if IS_SINGLE_END:
    rule prepare_wepp_inputs:
        input:
            r1 = lambda wc: f"{OUT_ROOT}/{ctx.dir_for_acc(wc.acc)}/{wc.acc.split('.')[0]}_R1.fq.gz",
            fasta = lambda wc: ctx.fasta_for(wc.acc),
            pb = lambda wc: ctx.pb_for(wc.acc),
            viz = f"{OUT_ROOT}/classification_proportions.png",
            classified = str(ctx.acc2classifieddir_json_path),
            kraken_out = str(ctx.kraken_out),
            kraken_report = str(ctx.kraken_report),
            coverage = str(ctx.acc2covered_json_path)
        output:
            new_r1 = f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}.fastq.gz",
        params:
            data_dir = lambda wc: f"{WEPP_DATA_DIR}/{wc.dir_tag}",
            tree_dest = lambda wc: ctx.wepp_tree(wc.acc),
            fasta_dest = lambda wc: ctx.wepp_ref(wc.acc),
            jsonl_src = lambda wc: ctx.jsonl_for(wc.acc),
            jsonl_dest = lambda wc: ctx.wepp_jsonl(wc.acc)
        resources:
            serial = 1
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

            if [ -n "{params.jsonl_src}" ]; then
                cp "{params.jsonl_src}" "{params.jsonl_dest}"
            fi
            """
else:
    rule prepare_wepp_inputs:
        input:
            r1 = lambda wc: f"{OUT_ROOT}/{ctx.dir_for_acc(wc.acc)}/{wc.acc.split('.')[0]}_R1.fq.gz",
            r2 = lambda wc: f"{OUT_ROOT}/{ctx.dir_for_acc(wc.acc)}/{wc.acc.split('.')[0]}_R2.fq.gz",
            fasta = lambda wc: ctx.fasta_for(wc.acc),
            pb = lambda wc: ctx.pb_for(wc.acc),
            viz = f"{OUT_ROOT}/classification_proportions.png",
            classified = str(ctx.acc2classifieddir_json_path),
            kraken_out = str(ctx.kraken_out),
            kraken_report = str(ctx.kraken_report),
            coverage = str(ctx.acc2covered_json_path)
        output:
            new_r1 = f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}_R1.fastq.gz",
            new_r2 = f"{WEPP_DATA_DIR}/{{dir_tag}}/{{acc}}_R2.fastq.gz",
        params:
            data_dir = lambda wc: f"{WEPP_DATA_DIR}/{wc.dir_tag}",
            tree_dest = lambda wc: ctx.wepp_tree(wc.acc),
            fasta_dest = lambda wc: ctx.wepp_ref(wc.acc),
            jsonl_src = lambda wc: ctx.jsonl_for(wc.acc),
            jsonl_dest = lambda wc: ctx.wepp_jsonl(wc.acc)
        resources:
            serial = 1
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

            if [ -n "{params.jsonl_src}" ]; then
                cp "{params.jsonl_src}" "{params.jsonl_dest}"
            fi
            """


rule run_wepp:
    input:
        r1 = lambda wc: f"{WEPP_DATA_DIR}/{wc.dir_tag}/{wc.acc}.fastq.gz",
        r2 = lambda wc: "" if IS_SINGLE_END else f"{WEPP_DATA_DIR}/{wc.dir_tag}/{wc.acc}_R2.fastq.gz",
    output:
        run_txt = f"{WEPP_RESULTS_DIR}/{{dir_tag}}/{{acc}}_run.txt"
    threads: workflow.cores
    params:
        snakefile  = str(ctx.wepp_workflow),
        workdir    = str(ctx.wepp_root),
        seq_type   = str(config.get("SEQUENCING_TYPE", "d")).lower(),
        primer_bed = config["PRIMER_BED"],
        min_af     = config["MIN_AF"],
        min_q      = config["MIN_Q"],
        max_reads  = config["MAX_READS"],
        min_len    = config["MIN_LEN"],
        cfgfile    = WEPP_CONFIG,
        resultsdir = WEPP_RESULTS_DIR,
        prefix     = lambda wc: wc.acc.split('.')[0],
        ref_name   = lambda wc: f"{wc.acc}.fa",
        tag_dir    = lambda wc: ctx.tag(wc.acc),
        tree_name  = lambda wc: os.path.basename(ctx.wepp_tree(wc.acc)),
        tree_full  = lambda wc: ctx.wepp_tree(wc.acc),
        fasta_name = lambda wc: os.path.basename(ctx.wepp_ref(wc.acc)),
        fasta_full = lambda wc: ctx.wepp_ref(wc.acc),
        cmd_log    = f"{WEPP_CMD_LOG}/{TAG}_dashboard_run.txt",
        pathogens_name = lambda wc: ctx.dir_for_acc(wc.acc)
    conda:
        "env/wepp.yaml"
    resources:
        serial = 1
    shell:
        r"""
        min_reads=100
        has_reads() {{
            local f=$1
            local num_reads
            num_reads=$( (gzip -cd "$f" 2>/dev/null || cat "$f") | wc -l )
            [ "$((num_reads / 4))" -ge "$min_reads" ]
        }}

        if ! has_reads "{input.r1}" && ( [ -z "{input.r2}" ] || ! has_reads "{input.r2}" ); then
            echo "Fewer than $min_reads for {wildcards.acc}; removing data folder and skipping WEPP."
            rm -rf "{WEPP_DATA_DIR}/{params.tag_dir}"
            mkdir -p {params.resultsdir}/{wildcards.acc}
            touch {output.run_txt}
            exit 0
        fi

        mkdir -p {WEPP_CMD_LOG}
        [ -f {params.cmd_log} ] || touch {params.cmd_log}

        python ./scripts/run_inner.py \
            --dir        {params.tag_dir} \
            --prefix     {params.prefix} \
            --primer_bed {params.primer_bed} \
            --min_af     {params.min_af} \
            --min_q      {params.min_q} \
            --max_reads  {params.max_reads} \
            --tree       {params.tree_name} \
            --ref        {params.fasta_name} \
            --snakefile  {params.snakefile} \
            --workdir    {params.workdir} \
            --cfgfile    {params.cfgfile} \
            --cores      {threads} \
            --sequencing_type {params.seq_type} \
            --pathogens_name {params.pathogens_name} \
            --cmd_log {params.cmd_log}

        touch {output.run_txt}
        """


rule emit_dashboard_instructions:
    input:
        runs = run_txts_for_all
    output:
        howto = f"{WEPP_CMD_LOG}/{TAG}_dashboard_howto.txt"
    params:
        script = "scripts/emit_dashboard_howto.py",
        cmd_log = f"{WEPP_CMD_LOG}/{TAG}_dashboard_run.txt"
    resources:
        serial = 1
    shell:
        "python scripts/emit_dashboard_howto.py --cmd-log {params.cmd_log} --out {output.howto}"


rule run_wepp_dashboard:
    input:
        run_txt = f"{WEPP_RESULTS_DIR}/{{dir_tag}}/{{acc}}_run.txt"
    output:
        dash_txt = f"{WEPP_RESULTS_DIR}/{{dir_tag}}/{{acc}}_dashboard_run.txt"
    threads: workflow.cores
    params:
        snakefile    = str(ctx.wepp_workflow),
        workdir      = str(ctx.wepp_root),
        seq_type     = str(config.get("SEQUENCING_TYPE", "d")).lower(),
        primer_bed   = config["PRIMER_BED"],
        min_af       = config["MIN_AF"],
        min_q        = config["MIN_Q"],
        max_reads    = config["MAX_READS"],
        clade_list   = config["CLADE_LIST"],
        clade_idx    = config["CLADE_IDX"],
        cfgfile      = WEPP_CONFIG,
        resultsdir   = WEPP_RESULTS_DIR,
        prefix       = lambda wc: wc.acc.split('.')[0],
        ref_name     = lambda wc: f"{wc.acc}.fa",
        tag_dir      = lambda wc: ctx.tag(wc.acc),
        tree_name    = lambda wc: os.path.basename(ctx.wepp_tree(wc.acc)),
        tree_full    = lambda wc: ctx.wepp_tree(wc.acc),
        fasta_name   = lambda wc: os.path.basename(ctx.wepp_ref(wc.acc)),
        fasta_full   = lambda wc: ctx.wepp_ref(wc.acc),
        cmd_log      = f"{WEPP_CMD_LOG}/{TAG}_dashboard_run.txt",
        pathogens_name = lambda wc: ctx.dir_for_acc(wc.acc),
        taxonium_file  = lambda wc: ctx.wepp_jsonl(wc.acc)
    conda:
        "env/wepp.yaml"
    resources:
        serial = 1
    shell:
        r"""
        if [ -n "{params.taxonium_file}" ]; then
            tax_arg="--taxonium_file {params.taxonium_file}"
        else
            tax_arg=""
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
            --snakefile  {params.snakefile} \
            --workdir    {params.workdir} \
            --cfgfile    {params.cfgfile} \
            --cores      {threads} \
            --sequencing_type {params.seq_type} \
            --pathogens_name {params.pathogens_name} \
            --cmd_log {params.cmd_log}
            $tax_arg

        touch {output.dash_txt}
        """
