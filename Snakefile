# metaWEPP Snakefile (refactored)

configfile: "config/config.yaml"

import subprocess
import sys

from workflow.context import build_context, ContextError


# ---------------------------------------------------------------------------
# CLI helpers and config validation
# ---------------------------------------------------------------------------
argv = sys.argv[1:]
requested_rules = {arg for arg in argv if arg and not arg.startswith("-")}
is_dry_run = any(arg in {"-n", "--dry-run", "--dryrun"} for arg in argv)

if "test" not in requested_rules:
    missing_keys = [key for key in ("DIR", "KRAKEN_DB") if not config.get(key)]
    if missing_keys:
        raise ValueError(
            "Missing required config entries: " + ", ".join(missing_keys)
        )

    if not is_dry_run:
        print("\nLaunching pathogen selection via add_ref_mat.py...", file=sys.stderr)
        subprocess.run(
            [
                "python",
                "scripts/add_ref_mat.py",
                "--db",
                str(config["KRAKEN_DB"]),
            ],
            check=True,
        )


# ---------------------------------------------------------------------------
# Build shared workflow context and expose it to rule modules
# ---------------------------------------------------------------------------
try:
    ctx = build_context(config)
except ContextError as exc:
    raise ValueError(str(exc))

workflow.globals["ctx"] = ctx
workflow.global_resources.update(serial=1)

if "test" not in requested_rules:
    print(
        "\n\n"
        f"        DIR={ctx.dir_name}, OUT_ROOT={ctx.out_root}\n"
        f"        KRAKEN_DB={ctx.kraken_db}\n"
        f"        SEQUENCING_TYPE={config.get('SEQUENCING_TYPE')}\n"
        f"        MIN_DEPTH_FOR_WEPP={config.get('MIN_DEPTH_FOR_WEPP')}\n"
        f"        DASHBOARD_ENABLED={config.get('DASHBOARD_ENABLED')}\n"
        "\n"
    )


# ---------------------------------------------------------------------------
# Shared Snakemake helper callables (wrapping context methods)
# ---------------------------------------------------------------------------

def split_fastqs_for_coverage(wc):
    return ctx.split_fastqs_for_coverage(checkpoints, wc)


def final_targets(wc):
    return ctx.final_targets(checkpoints, wc)


def final_targets_dashboard(wc):
    return ctx.final_targets_dashboard(checkpoints, wc)


def final_results_files(wc):
    return ctx.final_results_files(checkpoints, wc)


def run_txts_for_all(wc):
    return ctx.run_txts_for_all(checkpoints, wc)


def dashboard_howto_path(wc):
    return ctx.dashboard_howto_path(wc)


# ---------------------------------------------------------------------------
# Determine rule all inputs based on presence of classification output
# ---------------------------------------------------------------------------
if ctx.classified_json_empty():
    ALL_INPUTS = [
        str(ctx.acc2classifieddir_json_path),
        str(ctx.kraken_out),
        str(ctx.kraken_report),
        str(ctx.split_sentinel),
        str(ctx.visualization),
    ]
elif ctx.dashboard_enabled:
    ALL_INPUTS = (
        final_targets_dashboard,
        final_results_files,
        [
            str(ctx.acc2classifieddir_json_path),
            str(ctx.acc2covered_json_path),
            str(ctx.kraken_out),
            str(ctx.kraken_report),
            str(ctx.visualization),
        ],
    )
else:
    ALL_INPUTS = (
        final_targets,
        final_results_files,
        [
            str(ctx.acc2classifieddir_json_path),
            str(ctx.acc2covered_json_path),
            str(ctx.kraken_out),
            str(ctx.kraken_report),
            str(ctx.visualization),
            dashboard_howto_path,
        ],
    )


# ---------------------------------------------------------------------------
# Include rule definitions from modular .smk files
# ---------------------------------------------------------------------------
include: "workflow/rules/kraken.smk"
include: "workflow/rules/wepp.smk"


# ---------------------------------------------------------------------------
# Aggregate rule targets
# ---------------------------------------------------------------------------
rule all:
    input: ALL_INPUTS

rule test:
    run:
        pass