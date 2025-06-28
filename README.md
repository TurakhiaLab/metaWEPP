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

**Step 3:** Install Kraken.
The following commands install kraken and also update the `$PATH` variable for easily running the tool.
```
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh .
echo -e '\nexport PATH="$(pwd):$PATH"' >> ~/.bashrc
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

### <a name="mess"></a> Example - 1 RSV Dataset: Run the pipeline with MeSS simulated data

**Step 1:** Download the RSVA MAT, Reference FASTA File
```
wget https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2021/12/05/public-2021-12-05.all.masked.pb.gz

mkdir -p data/pathogens_for_detailed_analysis
cd data/pathogens_for_detailed_analysis
wget https://hgdownload.gi.ucsc.edu/hubs/GCF/002/815/475/GCF_002815475.1/UShER_RSV-A/2025/04/25/rsvA.2025-04-25.pb.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/475/GCF_002815475.1_ASM281547v1/GCF_002815475.1_ASM281547v1_genomic.fna.gz
gunzip GCF_002815475.1_ASM281547v1_genomic.fna.gz
mv GCF_002815475.1_ASM281547v1_genomic.fna NC_038235.fa
cd ../..
```

**Step 2:** Prepare simulated genomes
```
mkdir genomes
cd genomes
mv ../data/pathogens_for_detailed_analysis/NC_038235.fa .
```

**Step 3:** Build the Kraken database
```
mkdir test_kraken_DB
kraken2-build --download-taxonomy --db test_kraken_DB
k2 add-to-library --db test_kraken_DB --file /data/pathogens_for_detailed_analysis/NC_038235.fa 
kraken2-build --build --db test_kraken_DB
rm -rf test_kraken_DB/taxonomy #  To save disk memory
rm -rf test_kraken_DB/library  #  To save disk memory
```

**Step 4:**  Run the pipeline
```
snakemake --config SIMULATION_TOOL=MESS METAGENOMIC_REF=genomes/NC_038235.fa KRAKEN_DB=test_kraken_DB CLADE_IDX=1 --resources mess_slots=1 --cores 32
```

**Step 5:**  Analyze Results

All results can be found in the `WEPP/results/2697049` directory. This taxid is mapped to SARS-CoV-2, so analysis for this example is done on SARS-CoV-2. 

### <a name="real-world"></a> Example - 2: Real World Data

This example will take our own metagenomic wastewater reads and use them as input for our analysis, running on our RSVA dataset.

⚠️ Note that if you've already done Example 1, you may skip steps 1 and 2 for this example.

**Step 1:** Download the RSVA MAT, Reference FASTA File
```
wget https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2021/12/05/public-2021-12-05.all.masked.pb.gz

mkdir -p data/pathogens_for_detailed_analysis
cd data/pathogens_for_detailed_analysis
wget https://hgdownload.gi.ucsc.edu/hubs/GCF/002/815/475/GCF_002815475.1/UShER_RSV-A/2025/04/25/rsvA.2025-04-25.pb.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/475/GCF_002815475.1_ASM281547v1/GCF_002815475.1_ASM281547v1_genomic.fna.gz
gunzip GCF_002815475.1_ASM281547v1_genomic.fna.gz
mv GCF_002815475.1_ASM281547v1_genomic.fna NC_038235.fa
cd ../..
```

**Step 2:** Download astewater metagenomic reads:
```
mkdir -p data/RSVA_real
cd data/RSVA_real
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/011/ERR14763711/ERR14763711_*.fastq.gz https://hgdownload.gi.ucsc.edu/hubs/GCF/002/815/475/GCF_002815475.1/UShER_RSV-A/2025/04/25/rsvA.2025-04-25.pb.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/475/GCF_002815475.1_ASM281547v1/GCF_002815475.1_ASM281547v1_genomic.fna.gz
gunzip GCF_002815475.1_ASM281547v1_genomic.fna.gz
mv GCF_002815475.1_ASM281547v1_genomic.fna NC_038235.fa
mv ERR14763711_1.fastq.gz ERR14763711_R1.fastq.gz
mv ERR14763711_2.fastq.gz ERR14763711_R2.fastq.gz
cd ../..
```

**Step 3:** Build the Kraken database
```
mkdir test_kraken_DB
kraken2-build --download-taxonomy --db test_kraken_DB
k2 add-to-library --db test_kraken_DB --file data/RSVA_real/NC_038235.fa 
kraken2-build --build --db test_kraken_DB
rm -rf test_kraken_DB/taxonomy #  To save disk memory
rm -rf test_kraken_DB/library  #  To save disk memory
```
⚠️ Note that you must add the reference genome (in this example, `NC_038235.fa`) into the custom database for the pipeline to work.


**Step 4:**  Run the pipeline
```
snakemake --config FQ_DIR=RSVA_real SIMULATION_TOOL=none  KRAKEN_DB=test_kraken_DB CLADE_IDX=1 --resources mess_slots=1 --cores 32
```

**Step 5:**  Analyze Results

All results can be found in the `WEPP/results/NC_038235.1` directory. 


## <a name="usage"></a> Usage Guide:

### Data Organization

Visualization of META-WEPP's workflow directories
```text
📁 META-WEPP             
└───📁data                                      # [User Created] Contains data to analyze 
     ├───📁pathogens_for_wepp                   # [User Created] Pathogens selected for analysis
          ├───📁SARS_COV_2_real                   
               ├───sars_cov_2_reference.fa   
               ├───sars_cov_2_mat.pb.gz
          ├───📁RSV_A_real                   
               ├───rsv_a_2_reference.fa   
               ├───rsv_a_mat.pb.gz

     ├───📁metagenomic_sample_1                  # [User Created] Sample input reads (if providing reads)
          ├───metagenomic_reads_R1.fastq.gz      # Paired-ended reads
          ├───metagenomic_reads_R2.fastq.gz
     ├───📁metagenomic_sample_2                  
          ├───metagenomic_reference.fa           # [User Created] Sample input fasta file (if simulating reads)
          ├───metagenomic_reads_R1.fastq.gz      # [META-WEPP Generated] These are created after MeSS simulation
          ├───metagenomic_reads_R2.fastq.gz

└───📁results                                    # [META-WEPP Generated] Contains final META-WEPP results
      ├───📁meta_genomic_sample_1                # [META-WEPP Generated] Contains split reads
           ├───📁SARS_COV_2
                ├───sars_cov_2_reads_R1.fastq.gz    
                ├───sars_cov_2_reads_R2.fastq.gz
           ├───📁RSV_A
                ├───rsv_a_reads_R1.fastq.gz         
                ├───rsv_a_reads_R2.fastq.gz
      ├───📁meta_genomic_sample_2                        
           ├───📁SARS_COV_2
                ├───sars_cov_2_reads_R1.fastq.gz    
                ├───sars_cov_2_reads_R2.fastq.gz
           ├───📁RSV_A
                ├───rsv_a_reads_R1.fastq.gz         
                ├───rsv_a_reads_R2.fastq.gz
     
```

### Run Command

META-WEPP requires `KRAKEN_DB`, `TARGET_TAXIDS`, and `SIMULATE_TOOL` as config arguments through the command line, while the remaining ones can be taken from the config file. If you are setting `SIMULATE_TOOL=none`, then META-WEPP also requires `FQ1` and `FQ2` through the command line. It requires `--cores` from the command line, which is the number of threads used by the workflow, and also requires `--resources mess_slots=1` to prevent MeSS running in parallel which causes some issues.

Using all parameters from the config file:
```
snakemake --config SIMULATE_TOOL=MESS KRAKEN_DB=test_kraken_DB TARGET_TAXIDS=2697049 --resources mess_slots=1 --cores 32
```
Overriding `CLADE_IDX` and `PRIMER_BED`:
```
snakemake --config SIMULATE_TOOL=MESS KRAKEN_DB=test_kraken_DB TARGET_TAXIDS=2697049 CLADE_IDX=1 PRIMER_BED=none.bed --resources mess_slots=1 --cores 32
```

This will run the full pipeline and run WEPP for the taxid `2697049`, and uses the provided MAT and REF genome file.

### Arguments

META-WEPP has the following arguments:

1. `KRAKEN_DB` - Name of the Kraken database.
2. `SIMULATION_TOOL` - Input `"MESS"` to simulate reads with MeSS, or leave it blank, `""`, to provide your own reads.
3. `COVERAGE` - MESS's genomic coverage - Learn more about MESS's coverage calculation [here](https://metagenlab.github.io/MeSS/guide/simulate/coverage/).
4. `METAGENOMIC_REF` - Reference mixed fasta file if simulating with MeSS.
5. `CLADE_IDX` - Clade index for inferring lineages from MAT: Generally '1' for SARS-CoV-2 MAT and '0' for other MATs.
6. `PRIMER_BED` - BED file for primers. These are located in the `WEPP/primers` directory.
7. `SEQUENCING_TYPE` - `"d"` for paired end reads, `"s"` for single ended reads.

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
