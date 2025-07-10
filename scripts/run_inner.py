#!/usr/bin/env python3
import argparse
import subprocess
import os
import re

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
parser.add_argument("--clade_list", required=False, default="")  # Optional extra

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

# Optional values, using fallback if missing from config
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

print("Running command:", " ".join(cmd))
subprocess.run(cmd, check=True)
