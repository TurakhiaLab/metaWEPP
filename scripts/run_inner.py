#!/usr/bin/env python3
import argparse
import subprocess
import os

parser = argparse.ArgumentParser()
parser.add_argument("--dir",         required=True)
parser.add_argument("--prefix",      required=True)
parser.add_argument("--primer_bed",  required=True)
parser.add_argument("--tree",        required=True)
parser.add_argument("--ref",         required=True)
parser.add_argument("--clade_idx",   required=True)
parser.add_argument("--snakefile",   required=True)
parser.add_argument("--workdir",     required=True)
parser.add_argument("--configfile",  required=False)
parser.add_argument("--cores",       default="4")
parser.add_argument("--is_single_end", action="store_true",
                    help="If set, tell the inner workflow we are single-end")

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
    f"PRIMER_BED={args.primer_bed}",
    f"TREE={args.tree}",
    f"REF={args.ref}",
    f"CLADE_IDX={args.clade_idx}",
]

# add sequencing-type flag if single-end
if args.is_single_end:
    cmd.append("SEQUENCING_TYPE=s") 

if args.configfile:
    cmd += ["--configfile", args.configfile]

print("Running command:", " ".join(cmd))
subprocess.run(cmd, check=True)
