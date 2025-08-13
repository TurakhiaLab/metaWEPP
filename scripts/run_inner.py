#!/usr/bin/env python3
import argparse
import subprocess
import os
import re
import shlex
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--dir",         required=True)
parser.add_argument("--prefix",      required=True)
parser.add_argument("--primer_bed",  required=True)
parser.add_argument("--min_af",      required=True)
parser.add_argument("--min_q",       required=True)
parser.add_argument("--max_reads",   required=True)
parser.add_argument("--tree",        required=True)
parser.add_argument("--ref",         required=True)
parser.add_argument("--clade_idx",   required=True)
parser.add_argument("--snakefile",   required=True)
parser.add_argument("--workdir",     required=True)
parser.add_argument("--customconfig",required=False)
parser.add_argument("--cores",       default="32")
parser.add_argument("--sequencing_type", required=True,
                    help="Value to forward as SEQUENCING_TYPE in the inner config")
parser.add_argument("--clade_list", required=False, default="")
parser.add_argument("--cmd_log",     required=True)
parser.add_argument("--pathogens_name", required=True)
parser.add_argument("--dashboard_enabled", action="store_true",
                    help="If set, reuse the logged command for pathogens_name and run it with DASHBOARD_ENABLED=True")

args = parser.parse_args()

cmd = [
    "snakemake",
    "--snakefile",  args.snakefile,
    "--directory",  args.workdir,
    "--cores",      args.cores,
    "--use-conda",
    "--config",
    f"DIR={args.dir}",
    f"FILE_PREFIX={args.prefix}",
    f"TREE={args.tree}",
    f"REF={args.ref}",
]

def upsert_cmd_log(path, name, cmd_str):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    if path.exists():
        with open(path) as f:
            for line in f:
                parts = line.rstrip("\n").split(" : ", 1)
                if len(parts) == 2 and parts[0] == name:
                    continue  # drop old entry for this name
                lines.append(line.rstrip("\n"))
    lines.append(f"{name} : {cmd_str}")
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as f:
        f.write("\n".join(lines) + "\n")
    os.replace(tmp, path)

def parse_config_file(path):
    config_dict = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            match = re.match(r'^([^:#\s]+)\s*:\s*["\']?(.*?)["\']?\s*$', line)
            if match:
                key, value = match.groups()
                config_dict[key] = value
            else:
                raise ValueError(f"Invalid line in config file: {line}")
    return config_dict

if args.customconfig:
    config = parse_config_file(args.customconfig)
    cmd.append(f"PRIMER_BED={config.get('PRIMER_BED', args.primer_bed)}")
    cmd.append(f"MIN_AF={config.get('MIN_AF', args.min_af)}")
    cmd.append(f"MIN_Q={config.get('MIN_Q', args.min_q)}")
    cmd.append(f"MAX_READS={config.get('MAX_READS', args.max_reads)}")
    cmd.append(f"CLADE_IDX={config.get('CLADE_IDX', args.clade_idx)}")
    cmd.append(f"SEQUENCING_TYPE={config.get('SEQUENCING_TYPE', args.sequencing_type)}")
    if "CLADE_LIST" in config or args.clade_list:
        cmd.append(f"CLADE_LIST={config.get('CLADE_LIST', args.clade_list)}")
else:
    cmd.extend([
        f"PRIMER_BED={args.primer_bed}",
        f"MIN_AF={args.min_af}",
        f"MIN_Q={args.min_q}",
        f"MAX_READS={args.max_reads}",
        f"CLADE_IDX={args.clade_idx}",
        f"SEQUENCING_TYPE={args.sequencing_type}",
    ])
    if args.clade_list:
        cmd.append(f"CLADE_LIST={args.clade_list}")

cmd_log_path = Path(args.cmd_log)

if not args.dashboard_enabled:
    upsert_cmd_log(cmd_log_path, args.pathogens_name, " ".join(cmd))
    print("Running command:", " ".join(cmd))
    subprocess.run(cmd, check=True)
else:
    if not cmd_log_path.exists():
        raise SystemExit(f"cmd_log not found: {cmd_log_path}") # TODO: need to fix later

    matched_command_str = None
    with open(cmd_log_path, "r") as f:
        lines = f.readlines()
    # search last match
    for line in reversed(lines):
        line = line.rstrip("\n")
        parts = line.split(" : ", 1)
        if len(parts) == 2 and parts[0] == args.pathogens_name:
            matched_command_str = parts[1]
            break

    if matched_command_str is None:
        raise SystemExit(f"No command found in cmd_log for PATHOGENS_FOR_DASHBOARD '{args.pathogens_name}'")

    tokens = shlex.split(matched_command_str)

    # Put DASHBOARD_ENABLED=True into the --config block
    try:
        idx_config = tokens.index("--config")
        j = idx_config + 1
        k = j
        while k < len(tokens) and not tokens[k].startswith("-"):
            k += 1
        # replace if present
        replaced = False
        for m in range(j, k):
            if tokens[m].startswith("DASHBOARD_ENABLED="):
                tokens[m] = "DASHBOARD_ENABLED=True"
                replaced = True
                break
        if not replaced:
            tokens.insert(k, "DASHBOARD_ENABLED=True")
    except ValueError:
        tokens.extend(["--config", "DASHBOARD_ENABLED=True"])
    tokens.extend(["--forcerun", "dashboard_serve"])

    print("Running command (dashboard enabled):", " ".join(tokens))
    subprocess.run(tokens, check=True)
