import argparse
import subprocess
import os

parser = argparse.ArgumentParser()
parser.add_argument("--dir", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--primer_bed", required=True)
parser.add_argument("--tree", required=True)
parser.add_argument("--ref", required=True)
parser.add_argument("--clade_idx", required=True)
parser.add_argument("--snakefile", required=True)
parser.add_argument("--workdir", required=True)
parser.add_argument("--configfile", required=False)
parser.add_argument("--cores", default="4")

args = parser.parse_args()

cmd = [
    "snakemake",
    "--snakefile", args.snakefile,
    "--directory", args.workdir,
    "--cores", args.cores,
    "--use-conda",
    "--config",
    f'DIR={args.dir}',
    f'FILE_PREFIX={args.prefix}',
    f'PRIMER_BED={args.primer_bed}',
    f'TREE={args.tree}',
    f'REF={args.ref}',
    f'CLADE_IDX={args.clade_idx}'
]

if args.configfile:
    cmd += ["--configfile", args.configfile]

print("Running command:", " ".join(cmd))
subprocess.run(cmd, check=True)