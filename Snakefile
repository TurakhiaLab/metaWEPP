# Snakefile
import os
configfile: "config.yaml"
include: config["wepp_workflow"]


# module wepp:
#     snakefile: config["wepp_workflow"]
#     config: config["wepp_config"]

# use rule * from wepp as wepp_*


# 1) Load main pipeline config



# 2) Constants from config
DATASET_MAP = config["taxid_to_dataset"]     #  { '2697049':'NC_045512v2', … }
WEPP_DIR         = config["wepp_results_dir"]
CLASSIFIED_PREF  = config["classified_prefix"]
SIM_TOOL         = config.get("simulation_tool", "none").upper()
#WEPP_SNAKEFILE   = config["wepp_workflow"]   # e.g. /home/.../wepp/SARS2-WBE/workflow/Snakefile
WEPP_CONFIGFILE  = config["wepp_config"]     # e.g. /home/.../wepp/SARS2-WBE/config/config.yaml
TARGET_TAXIDS   = config["target_taxids"]
WEPP_SMK        = config["wepp_workflow"]
#WEPP_CFG        = config["wepp_config"]
WEPP_DIR_TPL    = config["wepp_results_dir"]
DATA_DIR     = "/home/jseangmany@AD.UCSD.EDU/wepp/WEPP/data"
WEPP_OUT   = config["wepp_results_dir"]   # "/…/results/{taxid}"
TAXIDS     = config["target_taxids"]

# 3) Top‐level: expect one WEPP/run.txt per taxid under WEPP_DIR
rule all:
    input:
        expand(
            "{wepp_dir}/{dataset}_run.txt",
            wepp_dir = config["wepp_results_dir"],
            dataset  = [ DATASET_MAP[t] for t in config["target_taxids"] ]
        ),
        expand(
            WEPP_DIR + "/{taxid}_run.txt",
            taxid=TARGET_TAXIDS
        ),
        expand(WEPP_OUT + "/{taxid}_run.txt", taxid=TAXIDS)

# 4) Optional ART simulation
if SIM_TOOL == "ART":
    rule simulate_reads_art:
        input:  "split_genomes.done"
        output: touch(config["simulate_reads_art_done"])
        shell:
            """
            mkdir -p art
            for fasta in individual_genomes/*.fa; do
                genome=$(basename "$fasta" .fa)
                art_illumina --paired --rndSeed 0 --noALN --maskN 0 --seqSys MSv3 \
                             --len 150 --rcount 800000 -m 200 -s 10 \
                             --in "$fasta" --out art/"$genome"
            done
            touch {output}
            """

    rule simulate_reads:
        input:  config["simulate_reads_art_done"]
        output: touch(config["simulate_reads_done"])
        shell: "touch {output}"

    rule merge_reads:
        input:  config["simulate_reads_done"]
        output:
            R1   = config["fq1"],
            R2   = config["fq2"],
            done = touch(config["merge_reads_done"])
        shell:
            """
            mkdir -p $(dirname {output.R1})
            cat art/*1.fq > {output.R1}
            cat art/*2.fq > {output.R2}
            touch {output.done}
            """

# 5) Optional MeSS simulation
elif SIM_TOOL == "MESS":
    rule simulate_reads_mess:
        input:  "mess_inputs/tsv_generation.done"
        output: touch(config["simulate_reads_mess_done"])
        shell:
            """
            mkdir -p simulated_reads
            for tsv in individual_genomes/*.tsv; do
                mess simulate --input "$tsv" \
                              --tech illumina \
                              --output simulated_reads/ \
                              --fasta individual_genomes/
            done
            touch {output}
            """

    rule simulate_reads:
        input:  config["simulate_reads_mess_done"]
        output: touch(config["simulate_reads_done"])
        shell: "touch {output}"

    rule merge_reads:
        input:  config["simulate_reads_done"]
        output:
            R1   = config["fq1"],
            R2   = config["fq2"],
            done = touch(config["merge_reads_done"])
        shell:
            """
            mkdir -p $(dirname {output.R1})
            cat simulated_reads/*_R1.fq > {output.R1}
            cat simulated_reads/*_R2.fq > {output.R2}
            touch {output.done}
            """

# 6) Kraken2 classification (conditionally depends on merge_reads)
if SIM_TOOL in ("ART","MESS"):
    kraken_input = dict(merge_done=config["merge_reads_done"],
                        r1=config["fq1"], r2=config["fq2"])
else:
    kraken_input = dict(r1=config["fq1"], r2=config["fq2"])

rule kraken:
    input:
        **kraken_input
    output:
        report     = config["kraken_report"],
        kraken_out = config["kraken_output"],
        done       = touch(config["kraken_done"])
    threads: config.get("kraken_threads",4)
    params:
        db = config["kraken_db"]
    shell:
        """
        mkdir -p $(dirname {output.report})
        kraken2 --db {params.db} --threads {threads} \
                --paired {input.r1} {input.r2} \
                --report {output.report} \
                --output {output.kraken_out} \
                --classified-out '/home/jseangmany@AD.UCSD.EDU/art/snakemake_art/kraken_reports/wepp_reads/classified_run2_R#'
        touch {output.done}
        """

# 7) Split Kraken output into per‐taxid FASTQs
rule split_per_taxid:
    input:
        kraken_out = config["kraken_output"],
        r1         = config["fq1"],
        r2         = config["fq2"]
    output:
        # one pair per taxid
        expand("/home/jseangmany@AD.UCSD.EDU/wepp/WEPP/data/{taxid}/{taxid}_R1.fq", taxid=config["target_taxids"]),
        expand("/home/jseangmany@AD.UCSD.EDU/wepp/WEPP/data/{taxid}/{taxid}_R2.fq", taxid=config["target_taxids"]),
        done = touch("split.done")
    params:
        script = "/home/jseangmany@AD.UCSD.EDU/art/snakemake_art/scripts/split_classified_reads.py",
        outdir = "/home/jseangmany@AD.UCSD.EDU/wepp/WEPP/data",
       # taxids = ",".join(map(str, config["target_taxids"]))
    shell:
        """
        python {params.script} \
          --kraken-out {input.kraken_out} \
          --r1 {input.r1} \
          --r2 {input.r2} \
          --out-dir {params.outdir} 
        """

# 8) Invoke WEPP’s Snakefile for each taxid


# pipeline config
WEPP_WD       = os.path.dirname(WEPP_SMK)   # WEPP’s working directory
WEPP_OUT_TPL  = config["wepp_results_dir"]  # e.g. "/home/.../wepp/SARS2-WBE/results/{taxid}"

rule run_wepp:
    input:
        r1      = config["wepp_data_dir"]   + "/{taxid}/{taxid}_R1.fq",
        r2      = config["wepp_data_dir"]   + "/{taxid}/{taxid}_R2.fq"
    output:
        run_txt = config["wepp_results_dir"] + "/{taxid}/{taxid}_run.txt"
    threads: config.get("wepp_threads", 4)
    shell:
        """
        snakemake \
          --snakefile {config[wepp_workflow]} \
          --directory {config[wepp_root]} \
          --cores {threads} \
          results/{wildcards.taxid}/{wildcards.taxid}_run.txt \
          --configfile {config[wepp_config]} \
          --config taxid={wildcards.taxid} \
                     wepp_data_dir={config[wepp_data_dir]} \
                     wepp_results_dir={config[wepp_results_dir]}
        """

# rule run_wepp:
#     input:
#         r1 = config["wepp_data_dir"] + "/{taxid}/{taxid}_R1.fq",
#         r2 = config["wepp_data_dir"] + "/{taxid}/{taxid}_R2.fq"
#     output:
#         run_txt = config["wepp_results_dir"] + "/{taxid}/{taxid}_run.txt"
#     threads: config.get("wepp_threads", 4)



# rule run_wepp:
#     input:
#         r1 = DATA_DIR + "/{taxid}/{taxid}_R1.fq.gz",
#         r2 = DATA_DIR + "/{taxid}/{taxid}_R2.fq.gz"
#     output:
#         run_txt = WEPP_OUT_TPL + "/{taxid}_run.txt"
#     params:
#         wepp_smk = WEPP_SMK,
#         wepp_cfg = WEPP_CFG,
#         wepp_wd  = WEPP_WD,
#         config_path = config["config_path"]  
#     threads: config.get("wepp_threads", 4)
#     shell:
#         r"""
#         mkdir -p $(dirname {output.run_txt})

#         cd {params.wepp_wd}

#         snakemake \
#           --snakefile {params.wepp_smk} \
#           --configfile {params.wepp_cfg} \
#           --config config_path={params.config_path} \
#           --cores {threads} \
#           --use-conda \
#           {output.run_txt} \
#           --config fq1={input.r1} fq2={input.r2} taxid={wildcards.taxid}
#         """
