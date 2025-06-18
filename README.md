# META-WEPP: Metagenomic Wastewater-based Epidemiology using Phylogenetic Placement

## Table of Contents
- [Introduction](#intro)
- [Installation](#install)
  - [Option-1: Install via Shell Commands](#shell)
- [Quick Start](#example)

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
git clone --recurse-submodules https://github.com/TurakhiaLab/metagenomic-WBE.git
cd metagenomic-WBE
```
**Step 2:** Install Conda (if your system does not have it already).
```
wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/download/24.11.3-2/Miniforge3-24.11.3-2-Linux-x86_64.sh"
bash Miniforge3.sh -b -p "${HOME}/conda"

source "${HOME}/conda/etc/profile.d/conda.sh"
source "${HOME}/conda/etc/profile.d/mamba.sh"
```
**Step 3:** Install Kraken.
```
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh $KRAKEN2_DIR
cp $KRAKEN2_DIR/kraken2{,-build,-inspect} $HOME/bin
```
Replace $KRAKEN2_DIR with the directory in which you would like to install Kraken2's scripts.

**Step 4:** Install MeSS.
```
conda create -n mess mess
```
**Step 5:** Install WEPP.

Follow the WEPP installation guide starting from option 3 step 2 on the [WEPP repo](https://github.com/TurakhiaLab/WEPP/tree/main?tab=readme-ov-file#-option-3-install-via-shell-commands-requires-sudo-access).

---

##  <a name="example"></a> Quick Start=

### <a name="MeSS"></a> Example - 1: Run the pipeline with MeSS simulated data
**Step 1:** Download the MAT:
```
wget https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2021/12/05/public-2021-12-05.all.masked.pb.gz
```
**Step 2:** Build the Kraken database
```
mkdir test_kraken_DB
kraken2-build --download-taxonomy --db test_kraken_DB
kraken2-build --add-to-library filtered_genomes.fa --db test_kraken_DB
kraken2-build --build --db test_kraken_DB
```

**Step 3:**  Run the pipeline
```
snakemake --cores 32
```

**Step 4:**  Analyze Results

All results can be found in the `WEPP/results/2697049` directory. This taxid is mapped to SARS-CoV-2, so analysis for this example is done on SARS-CoV-2.

---
## Usage Guide:

The entire pipeline can be ran through the following command:
```
snakemake --cores 32
```
This will run the full pipeline and run WEPP for each taxid in `target_taxids` in the `config.yaml` file.

The `config.yaml` file has the following arguments:


1. kraken_db - Name of the Kraken database.
2. kraken_report - Name of the Kraken report. (This tells you a report of what has been classified by Kraken)
3. kraken_output - Name of the Kraken output. (This which reads were mapped to the corresponding genome)
4. simulation_tool - Input "MeSS" to simulate reads with MeSS, or input "None" to provide your own reads.
5. REF - The reference mixed genome in fasta.
6. TREE - Mutation-Annotated Tree
7. target_taxids - The taxids to be analyzed.


If you would like more information on building a Kraken database, see below:
**Step 0:**  Build a Kraken database
**1.** Install a taxonomy. 
```
kraken2-build --download-taxonomy --db $DBNAME
```
(Replace "$DBNAME" above with your preferred database name/location. The database will use approximately 100 GB of disk space during creation. )

**2.** Add sequence to the database's genomic library using the --add-to-library switch, e.g.:
```
kraken2-build --add-to-library /path/to/chr1.fa --db $DBNAME
kraken2-build --add-to-library /path/to/chr2.fa --db $DBNAME
```

Add a list of files to the database's genomic library
```
for file in /path/to/chr*.fa
do
    kraken2-build --add-to-library $file --db $DBNAME
done
```

(Optional) Install one or more reference libraries: https://github.com/DerrickWood/kraken2/wiki/Manual#standard-kraken-2-database
```
kraken2-build --download-library bacteria --db $DBNAME
```

**3.** Build the database 
```
kraken2-build --build --db $DBNAME
```
Customized kmer with `--kmer-len` and `--minimizer-len` option if needed.

**4.** (Optional) Remove intermediate files after building a custom database which helps to free disk space.
```
kraken2-build --clean --db $DBNAME
```
More information in https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases
