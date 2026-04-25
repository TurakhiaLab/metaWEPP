import os
import csv, re, itertools
import subprocess
import pandas as pd
import shutil
import sys, glob
import time
from pathlib import Path
from collections import defaultdict
from itertools import zip_longest
import itertools, json

# ────────────────────────────────────────────────────────────────
# 1. DEFINE BASE DIRECTORY & LOAD CONFIG
# ────────────────────────────────────────────────────────────────

# This gets the directory where the Snakefile is located
BASE_DIR = Path(workflow.basedir)

# Load default config from the package location
configfile: str(BASE_DIR / "config/config.yaml")

# 2) Constants from config
DIR = config.get("DIR")
KRAKEN_DB = config.get("KRAKEN_DB")

requested_rules = set()
for arg in sys.argv[1:]:
    if not arg.startswith("-"):
        requested_rules.add(arg)

# ────────────────────────────────────────────────────────────────
# 2. DEFINE PATHS FOR INTERNAL RESOURCES
# ────────────────────────────────────────────────────────────────

RUNNING_HELP  = "help" in requested_rules

PATHOGEN_ROOT = Path("data/pathogens_for_wepp")
ADDED_TAXONS  = PATHOGEN_ROOT / "added_taxons.csv"

# Setting variables based on conda installation or docker image
env_root = Path(sys.prefix)
wepp_conda_path = env_root / "WEPP"
WEPP_DATA = Path.cwd() / "WEPP"

if wepp_conda_path.exists():
    WEPP_ROOT = wepp_conda_path
    WEPP_EXECUTABLE = env_root / "bin" / "run-wepp"
    if not WEPP_EXECUTABLE.exists():
        WEPP_EXECUTABLE = WEPP_ROOT / "run-wepp"
    WEPP_DATA.mkdir(parents=True, exist_ok=True)
else:
    WEPP_ROOT = WEPP_DATA
    WEPP_EXECUTABLE = WEPP_ROOT / "run-wepp"
    if not WEPP_ROOT.exists() and not RUNNING_HELP:
        print(f"Error: WEPP not installed at {WEPP_ROOT}", file=sys.stderr)
        sys.exit(1)

WEPP_SNAKEFILE = WEPP_ROOT / "workflow" / "Snakefile"

# User inputs are relative to where the user runs the command
IS_SINGLE_END = config.get("SEQUENCING_TYPE", "d").lower() in {"s", "n"}

rule all:
    input:
        f"results/{DIR}/classification_proportions.png",
        f"results/{DIR}/.split_read.done",
        f"results/{DIR}/run_wepp.txt",
        f"results/{DIR}/.wepp.done",
        f"results/{DIR}/.wepp_dashboard.done"

# ────────────────────────────────────────────────────────────────
# 3. HELPER FUNCTIONS
# ────────────────────────────────────────────────────────────────
PATHOGENS = [x.strip() for x in str(config.get("PATHOGENS") or "").split(",") if x.strip()] or ["default"]


def _parse_per_pathogen(raw, pathogens, default, caster=str, name="option"):
    """Parse a per-pathogen config value into a {pathogen: value} dict.

    Accepted forms:
      - None / empty string          → every pathogen uses ``default``.
      - Single value ``"x"``         → broadcast ``x`` to every pathogen.
      - ``N`` comma-separated values → 1-to-1 mapping with ``PATHOGENS``.

    Each option supplies its own sentinel via ``default``:
      - ``CLADE_LIST``         → ``""``   (no clade annotations; this is
        the literal per-species value, matching the original behaviour –
        e.g. the leading slot in ``",nextstrain:pango,..."`` becomes an
        empty string for the ``default`` pathogen)
      - ``CLADE_IDX``          → ``-1``   (no lineage annotations)
      - ``PRIMER_BED``         → ``""``   (no primer trimming for this
        species; absolute paths that only split on commas, never on slashes)
      - ``CORES_PER_PATHOGEN`` → ``"auto"`` sentinel, resolved later by
        the checkpoint so it can share leftover ``--cores`` across the
        species that asked for auto allocation.

    A mismatch between the number of values and ``len(PATHOGENS)`` raises
    a clear error pointing at the offending option.
    """
    if raw is None or str(raw).strip() == "":
        return {p: default for p in pathogens}

    parts = [x.strip() for x in str(raw).split(",")]

    def cast(token):
        if token == "":
            return default
        return caster(token)

    if len(parts) == 1:
        value = cast(parts[0])
        return {p: value for p in pathogens}

    if len(parts) != len(pathogens):
        raise ValueError(
            f"ERROR: '{name}' has {len(parts)} value(s) but PATHOGENS has "
            f"{len(pathogens)}. Provide either 1 value (broadcast to every "
            f"species) or exactly {len(pathogens)} comma-separated values."
        )

    return {p: cast(v) for p, v in zip(pathogens, parts)}


def _cast_cores(token):
    """Cast a CORES_PER_PATHOGEN token to either an int or the 'auto' sentinel."""
    s = str(token).strip().lower()
    if s in {"auto", ""}:
        return "auto"
    return max(1, int(s))


CLADE_LIST_MAP = _parse_per_pathogen(
    config.get("CLADE_LIST"), PATHOGENS, default="", name="CLADE_LIST"
)
CLADE_IDX_MAP = _parse_per_pathogen(
    config.get("CLADE_IDX"), PATHOGENS, default=-1, caster=int, name="CLADE_IDX"
)
PRIMER_BED_MAP = _parse_per_pathogen(
    config.get("PRIMER_BED"), PATHOGENS, default="", name="PRIMER_BED"
)
CORES_PER_PATHOGEN_MAP = _parse_per_pathogen(
    config.get("CORES_PER_PATHOGEN"),
    PATHOGENS,
    default="auto",
    caster=_cast_cores,
    name="CORES_PER_PATHOGEN",
)


def build_wepp_metadata(pathogens, clade_list_map, clade_idx_map, primer_bed_map, cores_map):
    """Build pathogen -> {clade_list, clade_idx, primer_bed, cores} metadata."""
    if not pathogens:
        raise ValueError("ERROR: PATHOGENS is required and must include 'default'.")
    if "default" not in pathogens:
        raise ValueError("ERROR: 'default' must be present in PATHOGENS.")

    meta = {}
    for pathogen in pathogens:
        meta[pathogen] = {
            "clade_list": str(clade_list_map.get(pathogen, "")).strip(),
            "clade_idx": int(clade_idx_map.get(pathogen, -1)),
            "primer_bed": str(primer_bed_map.get(pathogen, "")).strip(),
            "cores": cores_map.get(pathogen, "auto"),
        }
    return meta

if not RUNNING_HELP and os.environ.get("METAWEPP_CONFIG_PRINTED") != "1":
    print("\nmetaWEPP config arguments:", file=sys.stderr)
    for key in sorted(config.keys()):
        print(f"  {key}={config.get(key)}", file=sys.stderr)
    print("", file=sys.stderr)
    os.environ["METAWEPP_CONFIG_PRINTED"] = "1"

FQ1 = ""
FQ2 = ""

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

def check_input_files():
    if not DIR:
        raise ValueError("No DIR folder specified. Call snakemake with --config DIR=TEST_DIR")
    if not KRAKEN_DB:
        raise ValueError("Please provide KRAKEN_DB via --config KRAKEN_DB=<path>")

    # User data is expected in the CURRENT working directory's "data/" folder
    fq_dir = Path("data") / config["DIR"]
    if not fq_dir.exists():
        raise FileNotFoundError(f"Input folder {DIR} does not exist")

    if IS_SINGLE_END:
        r1_files = sorted(fq_dir.glob("*.fastq*"))
        if len(r1_files) != 1: raise RuntimeError("Single-end requires one file")
        fq1 = str(gzip_if_needed(r1_files[0]))
        fq2 = ""
    else:
        r1_files = sorted(fq_dir.glob("*_R1.fastq*"))
        r2_files = sorted(fq_dir.glob("*_R2.fastq*"))
        if len(r1_files) != 1 or len(r2_files) != 1: raise RuntimeError("Paired-end requires exactly one *_R1.fastq and *_R2.fastq")
        fq1 = str(gzip_if_needed(r1_files[0]))
        fq2 = str(gzip_if_needed(r2_files[0]))
    return fq1, fq2

if not RUNNING_HELP:
    FQ1, FQ2 = check_input_files()
    
# ────────────────────────────────────────────────────────────────
# 4. Auto-update added_taxons.csv using Kraken2 on ref FASTAs
# ────────────────────────────────────────────────────────────────

rule check_pathogen:
    output:
        ADDED_TAXONS
    params:
        pathogen_root = PATHOGEN_ROOT, 
        kraken_db = KRAKEN_DB
    run:
        # Logic remains mostly the same, just using absolute params.pathogen_root
        pathogen_root = Path(params.pathogen_root)
        added_taxons = Path(output[0])
        pathogen_root.mkdir(parents=True, exist_ok=True)
        
        # Load existing added_taxons.csv (if present)
        csv_dirs = set()
        if added_taxons.exists():
            with open(added_taxons, "r") as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split(",", 1)
                        if len(parts) == 2: csv_dirs.add(parts[1].strip())
        else:
            added_taxons.touch()

        # Scan pathogen directories
        new_lines = []
        for pathogen_dir in sorted(pathogen_root.iterdir()):
            if not pathogen_dir.is_dir():
                continue

            folder_name = pathogen_dir.name
            if folder_name in csv_dirs:
                continue

            fasta_files = (
                list(pathogen_dir.glob("*.fa")) +
                list(pathogen_dir.glob("*.fasta")) +
                list(pathogen_dir.glob("*.fna"))
            )

            if not fasta_files:
                continue

            if not (
                any(pathogen_dir.glob("*.pb")) or
                any(pathogen_dir.glob("*.pb.gz"))
            ):
                continue

            ref_fasta = fasta_files[0]

            # Run kraken2
            result = subprocess.run(
                ["kraken2", "--db", params.kraken_db, str(ref_fasta)],
                check=True,
                capture_output=True,
                text=True,
            )

            taxid = None
            for line in result.stdout.splitlines():
                parts = line.split()
                if len(parts) >= 3 and parts[0] == "C":
                    taxid = parts[2]
                    break

            if taxid is None:
                print(
                    f"[WARN] Could not find taxon id for '{ref_fasta}' "
                    f"using DB '{params.kraken_db}'. Skipping.",
                    file=sys.stderr,
                )
                continue

            new_lines.append((taxid, folder_name))
            csv_dirs.add(folder_name)

        # Append new entries to added_taxons.csv
        if new_lines:
            with open(added_taxons, "a") as f:
                for taxid, folder_name in new_lines:
                    f.write(f"{taxid},{folder_name}\n")

rule add_pathogen:
    input:
        rules.check_pathogen.output
    output:
        f"results/{DIR}/.add_pathogen.done"
    params:
        db = KRAKEN_DB,
        script = BASE_DIR / "scripts/add_ref_mat.py",
        add_species = config.get("ADD_SPECIES_RUNTIME", True)
    threads: workflow.cores
    run:
        if params.add_species:
            subprocess.run(
                ["python", str(params.script), "--db", params.db, "--threads", str(threads)],
                check=True
            )
        
        # Using the absolute PATHOGEN_ROOT defined at top
        if PATHOGEN_ROOT.exists() and PATHOGEN_ROOT.is_dir():
            for pathogen_dir in filter(Path.is_dir, PATHOGEN_ROOT.iterdir()):
                fasta_files = list(pathogen_dir.glob("*.[fF][aAn]*"))
                pb_files = list(pathogen_dir.glob("*.pb*"))
                num_fasta = len(fasta_files)
                num_pb = len(pb_files)
                if num_fasta == 1 and num_pb == 1:
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
        db = KRAKEN_DB,
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
    params:
        min_prop = float(config.get("MIN_PROP_FOR_WEPP", 0.01)),
        script = BASE_DIR / "scripts/kraken_data_visualization.py",
        add_species_runtime = config.get("ADD_SPECIES_RUNTIME", True)
    shell:
        r"""
        python {params.script} \
            {input.kraken_report} {output.png} {params.min_prop} {params.add_species_runtime}
        """

rule split_read:
    input:
        kraken_out = rules.kraken.output.kraken_out,
        kraken_report = rules.kraken.output.kraken_report,
        png = rules.kraken_visualization.output.png,
        r1 = FQ1,
        r2 = lambda wc: [] if IS_SINGLE_END else [FQ2],
    output:
        "results/{DIR}/pathogen_coverage.tsv",
        "results/{DIR}/.split_read.done"
    params:
        script = BASE_DIR / "scripts/split_read.py",
        sequencing_type = config.get("SEQUENCING_TYPE"),
    threads: workflow.cores
    run:
        output_dir = Path(output[0]).parent
        cmd = [
            "python",
            str(params.script),
            "-k", input.kraken_out,
            "-r", input.kraken_report,
            "--r1", input.r1,
            "-o", str(output_dir),
            "-t", str(threads),
            "-s", params.sequencing_type,
        ]
        if input.r2:
            cmd.extend(["--r2", input.r2[0]])

        subprocess.run(cmd, check=True)
        
        tsv_path = Path(output[0])
        if not tsv_path.exists(): tsv_path.write_text("")
        Path(output[1]).touch()

# ────────────────────────────────────────────────────────────────
# 5. RUN WEPP (PARALLEL ACROSS SPECIES)
#
# `prepare_for_wepp` is a checkpoint so the DAG can be expanded once
# species gets its own command file under results/{DIR}/wepp_cmds/ and
# is dispatched as an independent `run_wepp_single` job.
# ────────────────────────────────────────────────────────────────

checkpoint prepare_for_wepp:
    input:
        "results/{DIR}/.split_read.done",
        coverage = "results/{DIR}/pathogen_coverage.tsv"
    output:
        run_wepp_txt = "results/{DIR}/run_wepp.txt",
        cmd_dir = directory("results/{DIR}/wepp_cmds")
    run:
        cmd_dir = Path(output.cmd_dir)
        cmd_dir.mkdir(parents=True, exist_ok=True)

        try:
            coverage_df = pd.read_csv(input.coverage, sep="\t", header=None, names=["pathogen", "depth"])
        except pd.errors.EmptyDataError:
            coverage_df = pd.DataFrame(columns=["pathogen", "depth"])

        coverage_df = coverage_df.dropna()
        # If coverage_df is empty → create empty run_wepp.txt (no species to run)
        if coverage_df.empty:
            Path(output.run_wepp_txt).write_text("")
            return

        # Keep Pathogens meeting depth requirement
        MIN_DEPTH = float(config.get("MIN_DEPTH_FOR_WEPP", 0))
        selected_pathogens = coverage_df[coverage_df["depth"] >= MIN_DEPTH]["pathogen"].tolist()
        if not selected_pathogens:
            Path(output.run_wepp_txt).write_text("")
            return

        WEPP_METADATA = build_wepp_metadata(
            PATHOGENS,
            CLADE_LIST_MAP,
            CLADE_IDX_MAP,
            PRIMER_BED_MAP,
            CORES_PER_PATHOGEN_MAP,
        )
        DEFAULT_META = WEPP_METADATA["default"]

        # ── Resolve CORES_PER_PATHOGEN for every selected species ──
        # Explicit integer values are honoured verbatim and take precedence
        # shares whatever budget remains after subtracting the explicit allocations:
        #     remaining = --cores − Σ(explicit)
        #     auto_share = ceil(remaining / n_auto)   (clamped to ≥ 1)
        selected_requests = {
            p: WEPP_METADATA.get(p, DEFAULT_META).get("cores", "auto")
            for p in selected_pathogens
        }
        explicit_total = sum(
            int(v) for v in selected_requests.values() if v != "auto"
        )
        auto_pathogens = [p for p, v in selected_requests.items() if v == "auto"]
        if auto_pathogens:
            remaining = int(workflow.cores) - explicit_total
            auto_share = -(-remaining // len(auto_pathogens))  # ceil division
            if auto_share < 1:
                auto_share = 1
        else:
            auto_share = 1  # unused when no auto slots exist

        def resolve_cores(pathogen):
            requested = selected_requests[pathogen]
            if requested == "auto":
                return auto_share
            return max(1, int(requested))

        with open(output.run_wepp_txt, "w") as run_wepp_f:
            for pathogen in selected_pathogens:
                raw_meta = WEPP_METADATA.get(pathogen, DEFAULT_META)
                cl = (raw_meta["clade_list"] or "").replace(":", ",")
                clade_list_arg = f"CLADE_LIST={cl} " if cl else ""
                cl_idx = f"{raw_meta['clade_idx']}"
                primer_bed = raw_meta["primer_bed"]
                pathogen_cores = resolve_cores(pathogen)

                # Source directory (shared read-only pathogen assets)
                pathogen_dir = PATHOGEN_ROOT / pathogen
                
                # Destination Directory inside WEPP
                dest_dir = WEPP_DATA / "data" / f"{DIR}_{pathogen}"
                dest_dir.mkdir(parents=True, exist_ok=True)
                
                # 1. Copy files from data/pathogens_for_wepp/{pathogen}
                for file in pathogen_dir.glob("*"):
                    if file.is_file(): shutil.copy2(file, dest_dir)
                
                # 2. Copy files from results/{pathogen} if present
                result_dir = Path(f"results/{DIR}/{pathogen}")
                if result_dir.exists():
                    for file in result_dir.glob("*"):
                        if file.is_file(): shutil.copy2(file, dest_dir)

                tree_file = list(dest_dir.glob("*.pb")) + list(dest_dir.glob("*.pb.gz"))
                ref_file = list(dest_dir.glob("*.fa")) + list(dest_dir.glob("*.fasta")) + list(dest_dir.glob("*.fna"))

                if not tree_file: raise FileNotFoundError(f"No tree file in {pathogen_dir}")
                if not ref_file: raise FileNotFoundError(f"No ref file in {pathogen_dir}")

                tree_name = tree_file[0].name
                ref_name = ref_file[0].name

                taxonium_files = list(dest_dir.glob("*.jsonl")) + list(dest_dir.glob("*.jsonl.gz"))
                taxonium_arg = f"TAXONIUM_FILE={taxonium_files[0].name} " if taxonium_files else ""

                wepp_args = (
                    f"-s {WEPP_SNAKEFILE} "
                    f"--directory {WEPP_DATA} "
                    f"--use-conda "
                    f"--config DIR={DIR}_{pathogen} FILE_PREFIX=metaWEPP_run "
                    f"TREE={tree_name} REF={ref_name} "
                    f"SEQUENCING_TYPE={config.get('SEQUENCING_TYPE', 'd')} "
                    f"PRIMER_BED={primer_bed} "
                    f"MIN_AF={config.get('MIN_AF')} "
                    f"MIN_DEPTH={config.get('MIN_DEPTH')} "
                    f"MIN_Q={config.get('MIN_Q')} "
                    f"MIN_PROP={config.get('MIN_PROP')} "
                    f"MIN_LEN={config.get('MIN_LEN')} "
                    f"MAX_READS={config.get('MAX_READS')} "
                    f"{clade_list_arg}CLADE_IDX={cl_idx} {taxonium_arg}"
                    f"DASHBOARD_ENABLED={config.get('DASHBOARD_ENABLED')}"
                )

                # Full command (used for run_wepp.txt / dashboard instructions)
                full_cmd = f"{WEPP_EXECUTABLE} --cores {pathogen_cores} {wepp_args}\n"
                run_wepp_f.write(full_cmd)

                # Per-species command dispatched in parallel by Snakemake.
                # `--nolock` is required because concurrent WEPP invocations
                # share the same --directory (WEPP).
                per_cmd = (
                    f"{WEPP_EXECUTABLE} --nolock --cores {pathogen_cores} {wepp_args}\n"
                )
                (cmd_dir / f"{pathogen}.cmd").write_text(per_cmd)


def _cores_per_pathogen_threads(wildcards):
    """Schedule-side thread count for ``run_wepp_single``.

    Reads the ``--cores N`` that was baked into the per-species command file
    by the checkpoint. Snakemake requires ``threads <= workflow.cores``, so
    we cap the value here for scheduling purposes only — the WEPP
    subinvocation still receives the *requested* ``--cores N`` (unclamped),
    giving ``CORES_PER_PATHOGEN`` precedence over the parent ``--cores``.
    """
    cmd_path = Path(
        f"results/{wildcards.DIR}/wepp_cmds/{wildcards.pathogen}.cmd"
    )
    try:
        m = re.search(r"--cores\s+(\d+)", cmd_path.read_text())
        if m:
            return max(1, min(int(workflow.cores), int(m.group(1))))
    except Exception:
        pass
    return max(1, int(workflow.cores))


def _aggregate_wepp_done(wildcards):
    """Collect per-species .done files after the checkpoint expands the DAG."""
    cmd_dir = Path(
        checkpoints.prepare_for_wepp.get(DIR=wildcards.DIR).output.cmd_dir
    )
    pathogens = [p.stem for p in cmd_dir.glob("*.cmd")]
    return expand(
        "results/{DIR}/wepp_done/{pathogen}.done",
        DIR=wildcards.DIR,
        pathogen=pathogens,
    )


def _run_under_pty(cmd, log_path, env):
    """Run *cmd* in a child process attached to a pseudo-TTY and append
    every byte the child writes (stdout + stderr, in chronological order)
    to *log_path*. Returns the child's exit status.
    """
    import pty
    import termios

    master_fd, slave_fd = pty.openpty()
    try:
        attrs = termios.tcgetattr(slave_fd)
        attrs[1] &= ~termios.OPOST   # output flags
        attrs[3] &= ~termios.ECHO    # local flags
        termios.tcsetattr(slave_fd, termios.TCSANOW, attrs)
    except (termios.error, OSError):
        pass

    try:
        proc = subprocess.Popen(
            ["/bin/bash", "-c", cmd],
            stdin=slave_fd,
            stdout=slave_fd,
            stderr=slave_fd,
            env=env,
            close_fds=True,
        )
        os.close(slave_fd)
        slave_fd = -1

        with open(log_path, "ab") as logf:
            while True:
                try:
                    data = os.read(master_fd, 4096)
                except OSError:
                    break
                if not data:
                    break
                logf.write(data)
                logf.flush()
        return proc.wait()
    finally:
        if slave_fd >= 0:
            try:
                os.close(slave_fd)
            except OSError:
                pass
        try:
            os.close(master_fd)
        except OSError:
            pass


rule run_wepp_single:
    """Run WEPP for a single species. Multiple instances of this rule run
    concurrently; Snakemake schedules them based on ``threads`` vs
    ``--cores``, while the per-species ``--cores N`` baked into the command
    file is what WEPP actually uses (may exceed ``--cores``)."""
    input:
        cmd = "results/{DIR}/wepp_cmds/{pathogen}.cmd"
    output:
        done = "results/{DIR}/wepp_done/{pathogen}.done"
    threads: _cores_per_pathogen_threads
    run:
        base_cmd = Path(input.cmd).read_text().strip()
        if not base_cmd:
            Path(output.done).touch()
            return
        m = re.search(r"--cores\s+(\d+)", base_cmd)
        if m:
            requested_cores = int(m.group(1))
        else:
            requested_cores = max(1, int(workflow.cores))
            base_cmd = f"{base_cmd} --cores {requested_cores}"

        # Each parallel WEPP run gets its own log file so the parent terminal stays readable
        log_path = Path(f"results/{wildcards.DIR}/wepp_logs/{wildcards.pathogen}.log")
        log_path.parent.mkdir(parents=True, exist_ok=True)

        def fmt_elapsed(seconds):
            seconds = int(seconds)
            h, rem = divmod(seconds, 3600)
            m_, s = divmod(rem, 60)
            return f"{h}h{m_:02d}m{s:02d}s" if h else f"{m_}m{s:02d}s"

        print(
            f"[metaWEPP] Launching WEPP for '{wildcards.pathogen}' with "
            f"{requested_cores} cores (scheduler slot: {threads}) "
            f"-> {log_path}",
            file=sys.stderr,
        )
        sys.stderr.flush()

        # Header written before the PTY-attached child appends.
        log_path.write_text(
            f"# metaWEPP log for pathogen='{wildcards.pathogen}'\n"
            f"# requested cores: {requested_cores}\n"
            f"# command:\n{base_cmd}\n\n"
        )

        # Run the inner WEPP under a pseudo-TTY so tools that gate output
        # on isatty() (Snakemake, conda, tqdm, ...) emit the full chatter they would on a real terminal. 
        env = dict(os.environ)
        env.setdefault("TERM", "dumb")
        env.setdefault("NO_COLOR", "1")
        env.setdefault("PYTHONUNBUFFERED", "1")

        start_ts = time.time()
        try:
            rc = _run_under_pty(base_cmd, log_path, env)
            if rc != 0:
                raise subprocess.CalledProcessError(rc, base_cmd)
        except subprocess.CalledProcessError as exc:
            elapsed = fmt_elapsed(time.time() - start_ts)
            sep = "-" * 78
            try:
                full_log = log_path.read_text(errors="replace")
            except Exception:
                full_log = "<unable to read log file>"
            print(
                f"\n[metaWEPP] '{wildcards.pathogen}' WEPP FAILED after "
                f"{elapsed} (exit {exc.returncode}) -- full log: {log_path}\n"
                f"{sep}\n{full_log}\n{sep}\n",
                file=sys.stderr,
            )
            sys.stderr.flush()
            raise

        elapsed = fmt_elapsed(time.time() - start_ts)
        print(
            f"[metaWEPP] '{wildcards.pathogen}' WEPP complete in {elapsed} "
            f"(log: {log_path})",
            file=sys.stderr,
        )
        sys.stderr.flush()
        Path(output.done).parent.mkdir(parents=True, exist_ok=True)
        Path(output.done).touch()


rule run_wepp:
    """Aggregator: waits for every per-species WEPP job to finish."""
    input:
        _aggregate_wepp_done
    output:
        "results/{DIR}/.wepp.done"
    run:
        Path(output[0]).touch()

rule print_dashboard_instructions:
    input:
        f"results/{DIR}/.wepp.done"
    output:
        f"results/{DIR}/.wepp_dashboard.done"
    run:
        if not config.get("DASHBOARD_ENABLED"):
            print("\n\n\n\nTo view the dashboard, run the following commands:\n")
            with open(f"results/{DIR}/run_wepp.txt", "r") as f:
                for line in f:
                    line = line.strip()
                    if "DASHBOARD_ENABLED" in line:
                        modified_line = line.replace("DASHBOARD_ENABLED=False", "DASHBOARD_ENABLED=True")
                    else:
                        modified_line = line + " DASHBOARD_ENABLED=True"
                    print(modified_line, "\n")
            print("\n")
        
        Path(output[0]).touch()

rule help:
    message: "Printing metaWEPP configuration help"
    params:
        script = BASE_DIR / "scripts/metawepp_help.py"
    shell:
        "python {params.script}"