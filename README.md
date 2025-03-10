# jaden's branch
### What I did:
I've worked on creating an automated Snakemake workflow that splits mixed genome fasta files, runs Art on each genome, combines all generated reads, inputs them into Kraken for classification, and gets the accuracy of the classification.
### Kraken accuracy usage:
Calculating kraken accuracy currently isn't part of the Snakemake workflow, but there is a separate Python script `evaluate_kraken.py` that calculates it. The way that this script works is that it parses the expected taxids, and counts them in the kraken report. It then compares that to the total amount of classified reads as a ratio to get the classification accuracy. Here are some usage guidelines:
* The input must be the generated Kraken report from running Kraken, and the accession numbers that are mapped to the taxids. For example:
```
GROUND_TRUTH_TAXIDS = {
    "NC_000913.3": "511145",  # E. coli K-12 substr. MG1655
    "AF013254.1": "11250",    # genome name
    "ON811098.1": "2697049",  # genome name
    "OQ557947.1": "10244"     # genome name
}
```
For all simulated reads that were inputted into Kraken (and that you want to be checked for accuracy), you must fill in the taxid mapping.
* File directory should go into line 48, `kraken_counts, total_reads = load_kraken_report("<directory>")`.
### `split_fq.sh` usage:
Splits combined reads into separated reads, each with additional information of the taxid. The name of the file is of the following format: `<name_of_read>.fq`.
