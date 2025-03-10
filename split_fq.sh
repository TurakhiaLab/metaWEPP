#!/bin/bash
# separate_reads.sh
#
# This script processes two FASTQ files (R1 and R2) and extracts any records
# whose header starts with one of the target prefixes. For each such record,
# it appends the corresponding taxid to the header and writes the record into
# a new FASTQ file named according to the prefix.
#
# Define input FASTQ files
FASTQ_FILES=("output/combined_R1.fq" "output/combined_R2.fq")

# Process each FASTQ file using awk.
# The awk script reads the FASTQ file record by record (4 lines per record),
# checks if the header (line 1) starts with one of the defined prefixes, and if so:
#   - Appends " taxid:<taxid>" to the header.
#   - Outputs the record to an output file named <prefix>.fq.
awk '
BEGIN {
    # Define mapping: read prefix -> taxid.
    target["NC_000913.3"] = "511145"
    target["AF013254.1"] = "11250"
    target["ON811098.1"] = "2697049"
    target["OQ557947.1"] = "10244"
}
{
  # FASTQ records come in blocks of 4 lines.
  if ((NR - 1) % 4 == 0) {
    header = $0
    # Remove initial "@" and any trailing /1 or /2 for matching.
    id = header
    sub(/^@/, "", id)
    sub(/\/[12]$/, "", id)
    matched = 0
    # Loop over target prefixes.
    for (p in target) {
      if (index(id, p) == 1) {
         matched = 1
         tax = target[p]
         # Determine output file name based on prefix.
         outFile = p ".fq"
         # Append taxid to header.
         header = header " taxid:" tax
         # Print the record (4 lines) to the proper file.
         print header >> outFile
         getline; print >> outFile    # sequence line
         getline; print >> outFile    # plus line
         getline; print >> outFile    # quality line
         break
      }
    }
    # If not matched, skip the next three lines.
    if (!matched) {
      getline; getline; getline;
    }
  }
}
' "${FASTQ_FILES[@]}"
