# META-WEPP: Metagenomic Wastewater-Based Epidemiology using Phylogenetic Placements

## Table of Contents
- [Introduction](#intro)
- [Installation](#install)
  - [Option-1: Install via Shell Commands](#shell)
- [Quick Start](#example)
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

conda config --add channels conda-forge
conda config --add channels bioconda
```
**Step 3:** Install Kraken.
```
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh $KRAKEN2_DIR
cp $KRAKEN2_DIR/kraken2{,-build,-inspect} $HOME/bin
```
Replace `$KRAKEN2_DIR` with the directory in which you would like to install Kraken2's scripts.

**Step 4:** Install MeSS.
```
conda create -n mess mess=0.10.0
```
**Step 5:** Install WEPP.

Follow the WEPP installation guide starting from option 3 step 2 on the [WEPP repo](https://github.com/TurakhiaLab/WEPP/tree/main?tab=readme-ov-file#-option-3-install-via-shell-commands-requires-sudo-access).

---

##  <a name="example"></a> Quick Start

### <a name="MeSS"></a> Example - 1 SARS-CoV-2 Dataset: Run the pipeline with MeSS simulated data

This example will simulate reads from our `filtered_genomes.fa` mixed metagenomic fasta file using MeSS's `illumina` simulator.

**Step 1:** Download the SARS-CoV-2 MAT and Reference FASTA File:
```
wget https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2021/12/05/public-2021-12-05.all.masked.pb.gz
cp WEPP/NC_045512v2.fa ./genomes
```

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
snakemake --config kraken_db=test_kraken_DB target_taxids=2697049 --resources mess_slots=1 --cores 32
```

**Step 4:**  Analyze Results

All results can be found in the `WEPP/results/2697049` directory. This taxid is mapped to SARS-CoV-2, so analysis for this example is done on SARS-CoV-2. 

---
## <a name="usage"></a> Usage Guide:

META-WEPP requires `kraken_db` and `target_taxids` as config arguments through the command line, while the remaining ones can be taken from the config file. It also requires --cores from the command line, which is the number of threads used by the workflow.

Example 1:
```
snakemake --config kraken_db=test_kraken_DB target_taxids=2697049 --resources mess_slots=1 --cores 32
```
This will run the full pipeline and run WEPP for the taxid `2697049`.

The `config.yaml` file has the following arguments:


1. `kraken_db` - Name of the Kraken database.
2. `kraken_report` - Name of the Kraken report. (This tells you a report of what has been classified by Kraken)
3. `kraken_output` - Name of the Kraken output. (This which reads were mapped to the corresponding genome)
4. `simulation_tool` - Input "MeSS" to simulate reads with MeSS, or input "None" to provide your own reads.
5. `REF` - The reference genome in fasta.
6. `TREE` - Mutation-Annotated Tree
7. `target_taxids` - The taxids to be analyzed.
8. `mixed_genomes_fasta` - Reference mixed fasta file if simulating with MeSS

⚠️ If you are providing your own metagenomic wastewater reads, you must provide reference genomes (in the example above, `NC_045512v2.fa`) and a MAT.

⚠️ If you are simulating with MeSS, along with the reference genomes and MAT, you must also provide a reference mixed genome in fasta. In the quick start, the `filtered_genomes.fa` is the reference fasta file, and it must be a mixed genome sample.

---
##  <a name="build-database"></a> Building Kraken Databases
If you would like more information on building a Kraken database, see below:

### How to build a custom Kraken Database:

**Step 1:** Install a taxonomy. 
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
You can also add a multi fasta file in the genomic library.

⚠️ For this to work, the FASTA sequence headers must include either the NCBI accession numbers or the text `kraken:taxid` followed by the taxonomy ID for the genome. For example: `>sequence100|kraken:taxid|9606|`
```
kraken2-build --add-to-library /path/to/multi_fasta.fa --db $DBNAME
```

(Optional) Install one or more reference libraries: https://github.com/DerrickWood/kraken2/wiki/Manual#standard-kraken-2-database
```
kraken2-build --download-library bacteria --db $DBNAME
```
Installing the reference libraries can help if you do not have custom genomes you would like to input into the Kraken database.

**Step 3:** Build the database 
```
kraken2-build --build --db $DBNAME
```
Customize kmer with `--kmer-len` and `--minimizer-len` option if needed.

**Step 4:** (Optional) Remove intermediate files after building a custom database which helps to free disk space.
```
kraken2-build --clean --db $DBNAME
```
More information in https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases
