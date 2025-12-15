#!/usr/bin/env python3

"""Print metaWEPP configuration help similar to the WEPP CLI help output."""

from textwrap import TextWrapper
from typing import List
import sys


OPTIONS = [
    ("DIR", "Folder containing the metagenomic reads."),
    ("KRAKEN_DB", "Folder containing the Kraken2 database."),
    (
        "SEQUENCING_TYPE",
        "Sequencing read type (s: Illumina single-ended, d: Illumina double-ended, n: ONT long reads).",
    ),
    ("PRIMER_BED", "BED file for primers, expected under WEPP/primers."),
    (
        "MIN_AF",
        "Alleles with an allele frequency below this threshold are masked (Illumina: 0.5%, Ion Torrent: 1.5%, ONT: 2%).",
    ),
    ("MIN_DEPTH", "Sites with read depth below this threshold are masked by WEPP."),
    ("MIN_Q", "Alleles with a Phred score below this threshold are masked by WEPP."),
    (
        "MIN_PROP",
        "Minimum proportion of haplotypes detected by WEPP (Wastewater samples: 0.5%, Clinical samples: 5%).",
    ),
    ("MIN_LEN", "Minimum read length to consider after ivar trim (default: 80)."),
    (
        "MAX_READS",
        "Maximum number of reads considered by WEPP from the sample (useful for reducing runtime).",
    ),
    ("DASHBOARD_ENABLED", "Enables the WEPP dashboard for visualizing haplotype results."),
    (
        "PATHOGENS",
        "List of pathogens with custom WEPP settings; species not listed use default settings.",
    ),
    (
        "CLADE_LIST",
        "Comma-separated clade annotation schemes in the MAT file, ordered to match PATHOGENS; leave blank for species without clade annotations.",
    ),
    (
        "CLADE_IDX",
        "Comma-separated clade indices for each pathogen; use -1 for species without lineage annotations, ordered to match PATHOGENS.",
    ),
    ("MIN_DEPTH_FOR_WEPP", "Minimum read coverage required to run WEPP for any pathogen species."),
    (
        "MIN_PROP_FOR_WEPP",
        "Minimum relative abundance before metaWEPP prompts adding a species for haplotype-level analysis.",
    ),
]


def print_help() -> None:
    print("\nmetaWEPP configuration options\n")
    wrapper = TextWrapper(width=96, subsequent_indent=" " * 20)
    for name, description in OPTIONS:
        lines = wrapper.wrap(description)
        if not lines:
            continue
        print(f"{name:<18} {lines[0]}")
        for line in lines[1:]:
            print(f"{'':<18} {line}")
    print()


def main(argv: List[str]) -> int:
    if len(argv) == 0 or argv[0] in {"help", "--help", "-h"}:
        print_help()
        return 0

    print("Invalid command. Use 'help' to see available configuration options.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
