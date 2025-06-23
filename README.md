# META-WEPP: Metagenomic Wastewater-Based Epidemiology using Phylogenetic Placements

## Table of Contents
- [Introduction](#intro)
- [Installation](#install)
  - [Option-1: Install via Shell Commands](#shell)
- [Quick Start](#example)
  - [Example-1: Simulated Data](#mess)
  - [Example-2: Real World Data](#real-world) 
- [Usage Guide](#usage)
- [Building Kraken Databases](#build-database)

<br>


## <a name="intro"></a> Introduction

META-WEPP is a Snakemake-based bioinformatics pipeline designed to enable rapid classification and haplotype-level analysis of mixed-pathogen metagenomic samples. Developed for flexible, high-throughput use in public health surveillance, META-WEPP integrates [Kraken2](https://github.com/DerrickWood/kraken2) for taxonomic classification and routes identified pathogen reads into [WEPP](https://github.com/TurakhiaLab/WEPP) for phylogenetic placement and haplotype inference. The pipeline automates the end-to-end workflow—from raw mixed reads to lineage-level analysis—with optional support for simulated read generation using [MeSS](https://github.com/metagenlab/MeSS). 

<div align="center">
    <img src="docs/images/metawepp-figure.png" width="600">
    <div><b>Figure 1: META-WEPP Pipeline Visual</b></div>
</div>


## <a name="install"></a> Installation

### <a name="shell"></a> Option-1: Install via Shell Commands.

**Step 1:** Clone the repository.
```
git clone https://github.com/TurakhiaLab/metagenomic-WBE.git
cd metagenomic-WBE
```
**Step 2:** Install Conda (if your system does not have it already).
```
wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/download/24.11.3-2/Miniforge3-24.11.3-2-Linux-x86_64.sh"
bash Miniforge3.sh -b -p "${HOME}/conda"

source "${HOME}/conda/etc/profile.d/conda.sh"
source "${HOME}/conda/etc/profile.d/mamba.sh"
```

⚠️ Make sure that these two channels, `conda-forge` and `bioconda` is installed.

```
conda config --add channels conda-forge
conda config --add channels bioconda
```

**Step 3:** Install Kraken.
Replace `$KRAKEN2_DIR` with the directory in which you would like to install Kraken2's scripts. The following commands install kraken and also update the `$PATH` variable for easily running the tool.
```
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh $KRAKEN2_DIR
echo -e '\nexport PATH="$KRAKEN2_DIR:$PATH"' >> ~/.bashrc 
source ~/.bashrc
cd ..
```

**Step 4:** Install MeSS.

```
conda create -n mess mess=0.10.0
conda activate mess
```

**Step 5:** Install WEPP.

Follow the WEPP installation guide starting from option 3 on the [WEPP repo](https://github.com/TurakhiaLab/WEPP/tree/main?tab=readme-ov-file#-option-3-install-via-shell-commands-requires-sudo-access).

---

##  <a name="example"></a> Quick Start

### <a name="mess"></a> Example - 1 SARS-CoV-2 Dataset: Run the pipeline with MeSS simulated data

This example will simulate reads from our `filtered_genomes.fa` mixed metagenomic fasta file using MeSS's `illumina` simulator.

**Step 1:** Download the SARS-CoV-2 MAT and Reference FASTA File:
```
wget https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2021/12/05/public-2021-12-05.all.masked.pb.gz
```
Note that we have already provided the reference fasta file located in the `genomes` directory.

**Step 2:** Build the Kraken database
```
mkdir test_kraken_DB
kraken2-build --download-taxonomy --db test_kraken_DB
k2 add-to-library --db test_kraken_DB --file genomes/*.fa 
kraken2-build --build --db test_kraken_DB
```
⚠️ Note that you must add the reference genome (in this example, `NC_045512v2.fa`) into the custom database for the pipeline to work.

**Step 3:**  Run the pipeline
```
snakemake --config target_taxids=2697049 SIMULATE_TOOL=MeSS METAGENOMIC_REF=genomes/filtered_genomes.fa KRAKEN_DB=test_kraken_DB TREE=public-2021-12-05.all.masked.pb.gz PRIMER_BED=nimagenV2.bed CLADE_IDX=1 --resources mess_slots=1 --cores 32
```

**Step 4:**  Analyze Results

All results can be found in the `WEPP/results/2697049` directory. This taxid is mapped to SARS-CoV-2, so analysis for this example is done on SARS-CoV-2. 

### <a name="real-world"></a> Example - 2 SARS-CoV-2 Dataset: Run the pipeline with real world (non simulated) data

This example will take our own metagenomic wastewater reads and use them as input for our analysis.

⚠️ Note that if you've already done Example 1, you may skip steps 1 and 2 for this example.

**Step 1:** Download the SARS-CoV-2 MAT, Reference FASTA File, and wastewater metagenomic reads:
```
wget https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2021/12/05/public-2021-12-05.all.masked.pb.gz
```
Note that we have already provided the reference fasta file located in the `genomes` directory, and also our wastewater metagenomic reads located in the `example_metagenomic_reads` directory.



**Step 2:** Build the Kraken database
```
mkdir test_kraken_DB
kraken2-build --download-taxonomy --db test_kraken_DB
k2 add-to-library --db test_kraken_DB --file genomes/*.fa 
kraken2-build --build --db test_kraken_DB
```
⚠️ Note that you must add the reference genome (in this example, `NC_045512v2.fa`) into the custom database for the pipeline to work.

**Step 3:**  Run the pipeline
```
snakemake --config target_taxids=2697049 SIMULATE_TOOL=none FQ1=example_metagenomic_reads/mixed_reads_R1.fastq.gz FQ2=example_metagenomic_reads/mixed_reads_R2.fastq.gz KRAKEN_DB=test_kraken_DB TREE=public-2021-12-05.all.masked.pb.gz PRIMER_BED=nimagenV2.bed CLADE_IDX=1 --resources mess_slots=1 --cores 32
```

**Step 4:**  Analyze Results

All results can be found in the `WEPP/results/2697049` directory. This taxid is mapped to SARS-CoV-2, so analysis for this example is done on SARS-CoV-2. 

---
## <a name="usage"></a> Usage Guide:

META-WEPP requires `kraken_db` and `target_taxids` as config arguments through the command line, while the remaining ones can be taken from the config file. It also requires --cores from the command line, which is the number of threads used by the workflow.

Example 1:
```
snakemake --config kraken_db=test_kraken_DB target_taxids=2697049 TREE=public-2021-12-05.all.masked.pb.gz --resources mess_slots=1 --cores 32
```
This will run the full pipeline and run WEPP for the taxid `2697049`, and uses the provided MAT, `public-2021-12-05.all.masked.pb.gz`.

The `config.yaml` file has the following arguments:


1. `KRAKEN_DB` - Name of the Kraken database.
2. `KRAKEN_REPORT` - Name of the Kraken report. (This tells you a report of what has been classified by Kraken)
3. `KRAKEN_OUTPUT` - Name of the Kraken output. (This tells you which reads were mapped to the corresponding genome)
4. `SIMULATION_TOOL` - Input `"MeSS"` to simulate reads with MeSS, or input `"None"` to provide your own reads.
5. `COVERAGE` - MESS's genomic coverage - Learn more about MESS's coverage calculation [here](https://metagenlab.github.io/MeSS/guide/simulate/coverage/).
6. `REF` - The reference genome in fasta.
7. `TREE` - The Mutation-Annotated Tree.
8. `TARGET_TAXIDS` - The taxids to be analyzed.
9. `METAGENOMIC_REF` - Reference mixed fasta file if simulating with MeSS.
10. `CLADE_IDX` - Clade index for inferring lineages from MAT: Generally '1' for SARS-CoV-2 MAT and '0' for other MATs.
11. `PRIMER_BED` - BED file for primers. These are located in the `WEPP/primers` directory.
12. `FQ1` - R1 reads in `readname_R1.fastq.gz` format.
13. `FQ2` - R2 reads in `readname_R2.fastq.gz` format.

⚠️ If you are providing your own metagenomic wastewater reads, you must provide reference genomes (in the example above, `NC_045512v2.fa`) and a MAT.

⚠️ If you are simulating with MeSS, along with the reference genomes and MAT, you must also provide a reference mixed genome in fasta. In the quick start, the `filtered_genomes.fa` is the reference fasta file, and it must be a mixed genome sample.

---
##  <a name="build-database"></a> Building Kraken Databases
If you would like more information on building a Kraken database, see below:

### How to build a custom Kraken Database:

**Step 1:** Install the taxonomy. This is necessary for building Kraken databases.
```
kraken2-build --download-taxonomy --db $DBNAME
```
(Replace "$DBNAME" above with your preferred database name/location. The database will use approximately 100 GB of disk space during creation. )

**Step 2:** Add sequence to the database's genomic library using the --add-to-library switch, e.g.:
```
kraken2-build --add-to-library /path/to/chr1.fa --db $DBNAME
kraken2-build --add-to-library /path/to/chr2.fa --db $DBNAME
```

Add a list of files to the database's genomic library (all the .fa files in your current working directory)
```
k2 add-to-library --db test_kraken_DB --file *.fa
```
You can also add a multi fasta file (metagenomic fasta file) in the genomic library.

⚠️ For this to work, the FASTA sequence headers must include either the NCBI accession numbers or the text `kraken:taxid` followed by the taxonomy ID for the genome. For example: `>sequence100|kraken:taxid|9606|`
```
kraken2-build --add-to-library /path/to/multi_fasta.fa --db $DBNAME
```

(Optional) Install one or more [reference libraries](https://github.com/DerrickWood/kraken2/wiki/Manual#standard-kraken-2-database).
```
kraken2-build --download-library bacteria --db $DBNAME
```
Installing the reference libraries can help if you do not have custom genomes you would like to input into the Kraken database.

**Step 3:** Build the database 
```
kraken2-build --build --db $DBNAME
```
Customize kmer with `--kmer-len` and `--minimizer-len` option if needed. This may improve classification rate if tuned properly.

**Step 4:** (Optional) Remove intermediate files after building a custom database which helps to free disk space.
```
kraken2-build --clean --db $DBNAME
```
View more information at the official [Kraken2 documentation](https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases).
