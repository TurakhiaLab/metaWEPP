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
git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git 
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
The following commands install kraken and also update the `$PATH` variable for easily running the tool.
```
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh .
echo -e "\nexport PATH=\"$(pwd):\$PATH\"" >> ~/.bashrc
source ~/.bashrc
cd ..
```

**Step 4:** Install MeSS.

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n mess mess
conda activate mess
conda install -c conda-forge singularity
conda deactivate

echo 'export SINGULARITY_TMPDIR="$PWD/.singularity/tmp"' >> ~/.bashrc
echo 'export SINGULARITY_CACHEDIR="$PWD/.singularity/cache"' >> ~/.bashrc
echo 'mkdir -p "$SINGULARITY_TMPDIR" "$SINGULARITY_CACHEDIR"' >> ~/.bashrc
source ~/.bashrc
```

**Step 5:** Install WEPP.

Follow the WEPP installation guide starting from option 3 on the [WEPP repo](https://github.com/TurakhiaLab/WEPP/tree/main?tab=readme-ov-file#-option-3-install-via-shell-commands-requires-sudo-access).

---

##  <a name="example"></a> Quick Start

### <a name="mess"></a> Example - 1 SARS-CoV-2 Dataset: Run the pipeline with MeSS simulated data

**Step 1:** Download the RSVA MAT and Reference FASTA File.
```
wget https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2021/12/05/public-2021-12-05.all.masked.pb.gz

mkdir -p data/pathogens_for_wepp/sars_cov_2
mkdir data/simulated_metagenomic_sample

mv public-2021-12-05.all.masked.pb.gz data/pathogens_for_wepp/sars_cov_2
cp metagenomic_references/NC_045512.2.fasta data/pathogens_for_wepp/sars_cov_2
cp metagenomic_example.fa data/simulated_metagenomic_sample
```
**Step 2:** Build the Kraken database.
```
mkdir test_kraken_DB
kraken2-build --download-taxonomy --db test_kraken_DB
k2 add-to-library --db test_kraken_DB --file metagenomic_references/*
kraken2-build --build --db test_kraken_DB
```

**Step 3:**  Run the pipeline.
```
snakemake --config DIR=simulated_metagenomic_sample SIMULATION_TOOL=MESS KRAKEN_DB=test_kraken_DB CLADE_IDX=1 --resources mess_slots=1 --cores 32
```

**Step 4:**  Analyze Results.

All results can be found in the `WEPP/results/sars_cov_2` directory. This accession name is mapped to SARS-CoV-2, so analysis for this example is done on SARS-CoV-2. 

### <a name="real-world"></a> Example - 2: Real World Data

This example will take our own metagenomic wastewater reads and use them as input for our analysis, running on our RSVA dataset.

**Step 1:** Download the RSVA MAT and Reference FASTA File
```
mkdir -p data/pathogens_for_wepp/rsv_a
mkdir data/real_metagenomic_sample

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/011/ERR14763711/ERR14763711_*.fastq.gz https://hgdownload.gi.ucsc.edu/hubs/GCF/002/815/475/GCF_002815475.1/UShER_RSV-A/2025/04/25/rsvA.2025-04-25.pb.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/475/GCF_002815475.1_ASM281547v1/GCF_002815475.1_ASM281547v1_genomic.fna.gz 

gunzip GCF_002815475.1_ASM281547v1_genomic.fna.gz 
mv ERR14763711_1.fastq.gz data/real_metagenomic_sample/ERR14763711_R1.fastq.gz
mv ERR14763711_2.fastq.gz data/real_metagenomic_sample/ERR14763711_R2.fastq.gz
mv rsvA.2025-04-25.pb.gz data/pathogens_for_wepp/rsv_a
mv GCF_002815475.1_ASM281547v1_genomic.fna data/pathogens_for_wepp/rsv_a/GCF_002815475.1_ASM281547v1_genomic.fa
```

**Step 2:** Build the Kraken database
```
mkdir test_kraken_DB
kraken2-build --download-taxonomy --db test_kraken_DB
k2 add-to-library --db test_kraken_DB --file metagenomic_references/* 
kraken2-build --build --db test_kraken_DB
```
⚠️ Note that you must add the reference genome (in this example, `NC_045512v2.fa`) into the custom database for the pipeline to work. This is done for us in the third line of Step 2.


**Step 3:**  Run the pipeline
```
snakemake --config DIR=real_metagenomic_sample KRAKEN_DB=test_kraken_DB CLADE_IDX=1 --resources mess_slots=1 --cores 32
```

**Step 4:**  Analyze Results

All results can be found in the `WEPP/results/rsva_a` directory. 


## <a name="usage"></a> Usage Guide:

### Data Organization

Visualization of META-WEPP's workflow directories
```
📁 META-WEPP             
└───📁data                                       # [User Created] Contains data to analyze 
     ├───📁pathogens_for_wepp                    # [User Created] Pathogens for Variant Analysis
          ├───📁SARS_COV_2_real                   
               ├───ON811098.fa                   # SARS COV 2 reference genome
               ├───ON811098.pb.gz                # SARS COV 2 mat

          ├───📁RSV_A_real                   
               ├───NC_045512.fa                  # RSV A reference genome
               ├───NC_045512.pb.gz               # RSV A mat

     ├───📁real_metagenomic_sample               # [User Created] Folder containing wastewater reads
          ├───metagenomic_reads_R1.fastq.gz      
          ├───metagenomic_reads_R2.fastq.gz

     ├───📁simulated_metagenomic_sample          # [User Created] Folder containing metagenomic fasta file for simulating reads with MeSS               
          ├───metagenomic_reference.fa           

     ├───📁simulated_reads                       # [META-WEPP Generated] Reads generated by MeSS
          ├───📁fastq   
              ├───merged_R1.fastq.gz 
              ├───merged_R2.fastq.gz

└───📁results                                    # [META-WEPP Generated] Contains pathogens specific reads found by META-WEPP
      ├───📁real_metagenomic_sample                 
           ├───📁ON811098
                ├───ON811098_R1.fastq.gz    
                ├───ON811098_R2.fastq.gz

           ├───📁NC_045512
                ├───NC_045512.fastq.gz         
                ├───NC_045512.fastq.gz

      ├───📁simulated_metagenomic_sample                        
           ├───📁ON811098
                ├───ON811098_R1.fastq.gz    
                ├───ON811098_R2.fastq.gz

           ├───📁NC_045512
                ├───NC_045512.fastq.gz         
                ├───NC_045512.fastq.gz
     
```

### Run Command

META-WEPP requires `KRAKEN_DB` and `DIR` as config arguments through the command line, while the remaining ones can be taken from the config file. It requires `--cores` from the command line, which is the number of threads used by the workflow, and also requires `--resources mess_slots=1` to prevent MeSS running in parallel which causes some issues.

Using all parameters from the config file:
```
snakemake --config KRAKEN_DB=test_kraken_DB DIR=simulated_metagenomic_sample --resources mess_slots=1 --cores 32
```
Overriding `CLADE_IDX`, `PRIMER_BED`, and `SIMULATE_TOOL`:
```
snakemake --config SIMULATE_TOOL=MESS KRAKEN_DB=test_kraken_DB DIR=simulated_metagenomic_sample CLADE_IDX=1 PRIMER_BED=none.bed --resources mess_slots=1 --cores 32
```

### Arguments

META-WEPP requries the following arguments:

1. `KRAKEN_DB` - Name of the Kraken database.
2. `SIMULATION_TOOL` - Input `"MESS"` to simulate reads with MeSS, or don't include the command in the command-line argument to provide your own real reads.
3. `DIR` - Directory of the metagenomic reference fasta file (simulation) or real metagenomic reads (providing your own reads).
4. `COVERAGE` - MESS's genomic coverage - Learn more about MESS's coverage calculation [here](https://metagenlab.github.io/MeSS/guide/simulate/coverage/).
5. `METAGENOMIC_REF` - Reference mixed fasta file if simulating with MeSS.
6. `CLADE_IDX` - Clade index for inferring lineages from MAT: Generally '1' for SARS-CoV-2 MAT and '0' for other MATs.
7. `PRIMER_BED` - BED file for primers. These are located in the `WEPP/primers` directory.
8. `SEQUENCING_TYPE` - Sequencing read type (s:Illumina single-ended, d:Illumina double-ended, or n:ONT long reads)

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

⚠️ If you would like to save disk memory, perform the following commands:
```
rm -rf test_kraken_DB/taxonomy 
rm -rf test_kraken_DB/library  
```
View more information at the official [Kraken2 documentation](https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases).
