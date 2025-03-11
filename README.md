# jaden's branch
### What I did:
I've worked on creating an automated Snakemake workflow that splits mixed genome fasta files, runs Art on each genome, combines all generated reads, inputs them into Kraken for classification, and gets the accuracy of the classification.
## Snakefile
This automated workflow does the same workflow mentioned above: Splits input mixed genome fasta file, runs Art on each genome, combines all reads, inputs them into classification, gets the accuracy of the classification, and splits the merged reads into separate fastq files. Inside the `Snakefile` file, there are notes about what is the input to each rule, and what should be the output. Here are some **important** usage guidelines.
* Directories must be written inside the `Snakefile` manually. Future steps will be to accept directories as an argument.
* Directories must follow the input output flow. For example, some directories will be created, like the `output` directory. Within the `output` directory is where the merged reads will be stored. The next rule runs Kraken on those merged reads, and it takes in as input the reads **inside** the `output` directory. So if you were to change where the reads will be stored, it will cause the other rules to not work properly. For full effectiveness and understanding, please read each comment before each rule in `Snakefile`, as there are explicit information on which directories are used as input, and which directories are used as output.
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
Future steps will be to automate this process and get these taxid mapping automatically. For now, it must be written inside the `evaluate_kraken.py` file manually.
For all simulated reads that were inputted into Kraken (and that you want to be checked for accuracy), you must fill in the taxid mapping.
* File directory should go into line 48, `kraken_counts, total_reads = load_kraken_report("<directory>")`.
### `split_fq.sh` usage:
Splits combined reads into separated reads, each with additional information of the taxid. The name of the file is of the following format: `<name_of_read>.fq`. However, it does need the following mapping inputted manually. Future steps will be to automate this process and get these mappings. Here is an example:
```
    # Define mapping: read prefix -> taxid.
    target["NC_000913.3"] = "511145"
    target["AF013254.1"] = "11250"
    target["ON811098.1"] = "2697049"
    target["OQ557947.1"] = "10244"
```
It must be implemented manually within the `split_fq.sh` file.
