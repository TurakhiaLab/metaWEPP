<div align="center">

# metaWEPP: Metagenomic Wastewater-Based Epidemiology using Phylogenetic Placements

[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://github.com/TurakhiaLab/WEPP/blob/main/LICENSE

[![License][license-badge]][license-link]
[<img src="https://img.shields.io/badge/Build with-CMake-green.svg?logo=snakemake">](https://cmake.org)
[<img src="https://img.shields.io/badge/Made with-Snakemake-aquamarine.svg?logo=snakemake">](https://snakemake.readthedocs.io/en/v7.19.1/index.html)

</div>

## Table of Contents
- [Introduction](#intro)
- [Installation](#install)
  - [Option-1: Install via Dockerfile](#docker)
  - [Option-2: Install via Shell Commands](#shell)
- [Quick Start](#example)
  - [Example-1: Simulated Data](#mess)
  - [Example-2: Real World Data](#real-world) 
- [User Guide](#guide)
  - [Data Organization](#data)
  - [Arguments](#arg)
  - [Run Command](#run) 
  - [MAT Download](#mat) 
- [Building Kraken Database](#build-database)

<br>


## <a name="intro"></a> Introduction

metaWEPP is a Snakemake-based bioinformatics pipeline that achieves haplotype-level resolution from metagenomic sequencing data. As illustrated in Figure 1, metaWEPP can analyze mixed sequencing read samples originating from diverse sources, including environmental metagenomes and clinical samples containing co-infections. It leverages [Kraken2](https://github.com/DerrickWood/kraken2) for taxonomic classification and then employs [WEPP](https://github.com/TurakhiaLab/WEPP) to determine haplotypes for each identified pathogen, enabling high-resolution variant analysis of metagenomic samples. The pipeline also integrates WEPP’s interactive visualization dashboard, which allows in-depth exploration of detected haplotypes for each pathogen — offering unprecedented fine-grained insight into metagenomic samples. For benchmarking and testing, metaWEPP incorporates the metagenomic simulation tool, [MeSS](https://github.com/metagenlab/MeSS). 

<div align="center">
    <img src="metaWEPP.pdf" width="600">
    <div><b>Figure 1: metaWEPP Pipeline</b></div>
</div>


## <a name="install"></a> Installation

### <a name="docker"></a> Option-1: Install via Dockerfile.

**Step 1:** Clone the repository.
```
git clone https://github.com/TurakhiaLab/metaWEPP.git
cd metaWEPP
```

**Step 2:** Build a Docker Image.
```
cd docker
docker build -t metawepp . 
cd ..
```
**Step 3:** Start and run Docker container. The command below will take you inside the Docker container with metaWEPP already installed. 
```
docker run -it metawepp
```
⚠️ Note: MeSS cannot be used to simulate reads within the Docker container, as it internally relies on Singularity, which causes issues with the Docker container.


### <a name="shell"></a> Option-2: Install via Shell Commands (requires sudo access).

**Step 1:** Clone the repository.
```
git clone https://github.com/TurakhiaLab/metaWEPP.git
cd metaWEPP
```

**Step 2:** Install Kraken.
The following commands install kraken and also update the `$PATH` variable for running the tool easily.
```
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh .
echo -e "\nexport PATH=\"$(pwd):\$PATH\"" >> ~/.bashrc
source ~/.bashrc
cd ..
```
**Step 3:** Install WEPP.

```
git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git
```
View the WEPP installation guid starting from option 3 on the [WEPP repo](https://github.com/TurakhiaLab/WEPP/tree/main?tab=readme-ov-file#-option-3-install-via-shell-commands-requires-sudo-access).

**Step 4:** Install MeSS.

Follow the MeSS installation guide on the [MeSS Quick Start](https://github.com/metagenlab/MeSS?tab=readme-ov-file#zap-quick-start).

---

##  <a name="example"></a> Quick Start

### <a name="mess"></a> Example - 1: MeSS simulated metagenomic sample

**Step 1:** Download the SARS-CoV-2 MAT and the corresponding reference file for its variant analysis with WEPP.
```
wget https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2021/12/05/public-2021-12-05.all.masked.pb.gz

mkdir -p data/pathogens_for_wepp/sars_cov_2
mkdir data/simulated_metagenomic_sample

mv public-2021-12-05.all.masked.pb.gz data/pathogens_for_wepp/sars_cov_2
cp metagenomic_references/NC_045512.2.fasta data/pathogens_for_wepp/sars_cov_2
cp metagenomic_example.fa data/simulated_metagenomic_sample
```
**Step 2:** Prepare the config.yaml for SARS-CoV-2 variant analysis.
```
cat <<EOF > data/pathogens_for_wepp/sars_cov_2/config.yaml
PRIMER_BED: none.bed
CLADE_IDX: 1
SEQUENCING_TYPE: d
CLADE_LIST: nextstrain,pango
EOF
```

**Step 3:** Build a custom Kraken database for taxonomic classification.
```
mkdir test_kraken_DB
kraken2-build --download-taxonomy --db test_kraken_DB
k2 add-to-library --db test_kraken_DB --file metagenomic_references/*
kraken2-build --build --db test_kraken_DB
```
⚠️ Note that the reference genome of any pathogen being analyzed with WEPP must be included in the Kraken database.


**Step 4:**  Run the pipeline.
```
snakemake --config DIR=simulated_metagenomic_sample SIMULATION_TOOL=MESS KRAKEN_DB=test_kraken_DB --resources mess_slots=1 --cores 32
```

**Step 5:**  Analyze Results.

Taxonomic classification results can be viewed at `results/simulated_metagenomic_sample/classification_proportions.png`, while all results generated by WEPP for SARS-CoV-2 variant analysis are present in `WEPP/results/sars_cov_2_simulated_metagenomic_sample`. 

Expected metagenomic classification proportions:

- 54.10% Severe acute respiratory syndrome coronavirus 2
- 25.25% Human respiratory syncytial virus
- 19.13% Dengue virus type 4
- 1.53% Unclassified

Expected SARS-CoV-2 lineage abundance results, which can be viewed at `WEPP/results/sars_cov_2_simulated_metagenomic_sample/NC_045512_lineage_abundance.csv`:

- AY.103,1.000000

### <a name="real-world"></a> Example - 2: Real world metagenomic sample

**Step 1:** Download the RSV-A MAT and the corresponding reference file for its variant analysis with WEPP.
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

**Step 2:** Prepare the config.yaml for RSV-A variant analysis.
```
cat <<EOF > data/pathogens_for_wepp/rsv_a/config.yaml
PRIMER_BED: RSVA_all_primers_best_hits.bed
CLADE_IDX: 0
SEQUENCING_TYPE: d
CLADE_LIST: annotation_1
EOF
```

**Step 3:** Build a custom Kraken database for taxonomic classification.
```
mkdir test_kraken_DB
kraken2-build --download-taxonomy --db test_kraken_DB
k2 add-to-library --db test_kraken_DB --file metagenomic_references/* 
kraken2-build --build --db test_kraken_DB
```
⚠️ Note that the reference genome of any pathogen being analyzed with WEPP must be included in the Kraken database.

**Step 4:**  Run the pipeline.
```
snakemake --config DIR=real_metagenomic_sample KRAKEN_DB=test_kraken_DB --resources mess_slots=1 --cores 32
```

**Step 5:**  Analyze Results.

Taxonomic classification results can be viewed at `results/real_metagenomic_sample/classification_proportions.png`, while all results generated by WEPP for RSV-A variant analysis are present in `WEPP/results/rsva_a_real_metagenomic_sample`. 

Expected metagenomic classification proportions:

- 57.01% Human respiratory syncytial virus
- 42.99% Unclassified

Expected RSV-A lineage abundance results, which can be viewed at `WEPP/results/rsv_a_real_metagenomic_sample/NC_038235_lineage_abundance.csv`:
- A.D.3,0.450398
- A.D.1,0.419891
- A.D,0.094283
- A.D.5,0.035429

## <a name="guide"></a> User Guide

### <a name="data"> Data Organization
We assume that all wastewater samples are stored in the data directory, each within its own subdirectory given by DIR argument (see Run Command). For each pathogen to be analyzed for variants with WEPP, place its reference genome, corresponding MAT file, and the config.yaml (optional) in the folder: `data/pathogens_for_wepp/<pathogen_name>`. Moreover, each created `DIR` inside `data` is expected to contain either of these files.
1. Sequencing Reads: Ending with `*_R{1/2}.fastq.gz` for paired-ended reads and `*.fastq.gz` for single-ended.
OR
2. Genomes: Single `fasta` file containing all the genomes to be simulated in the sample. MeSS will generate reads and place it in the same folder.

For each sample, metaWEPP stores metagenomic analysis results in corresponding subdirectories under `results`. Variant-specific analysis outputs are located within the respective pathogen directories under `WEPP/results`. 

Visualization of metaWEPP's workflow directories
```
📁 metaWEPP             
└───📁data                                       # [User Created] Contains data to analyze 
     ├───📁pathogens_for_wepp                    # [User Created] Pathogens for Variant Analysis with WEPP
          ├───📁SARS_COV_2                
               ├───SARS_COV_2_ref_genome.fa     
               ├───SARS_COV_2_mat.pb.gz 
               ├───config.yaml                   # Customized WEPP config for SARS-COV-2

          ├───📁RSV_A                   
               ├───RSV_A_ref_genome.fa     
               ├───RSV_A_mat.pb.gz 
               ├───config.yaml                   #Customized WEPP config for RSV-A

     ├───📁real_metagenomic_sample               # [User Created] Contains metagenomic sample for analysis
          ├───metagenomic_reads_R1.fastq.gz      
          ├───metagenomic_reads_R2.fastq.gz

     ├───📁simulated_metagenomic_sample          # [User Created] Contains metagenomic fasta for simulating reads with MeSS               
          ├───metagenomic_reference.fa           
          ├───MeSS_R1.fastq.gz                   # [MeSS Generated] Simulated reads by MeSS
          ├───MeSS_R2.fastq.gz

└───📁results                                    # [metaWEPP Generated] 
      ├───📁real_metagenomic_sample                 
           ├───📁SARS_COV_2
                ├───SARS_COV_2_R1.fastq.gz    
                ├───SARS_COV_2_R2.fastq.gz

           ├───📁RSV_A
                ├───RSV_A_R1.fastq.gz         
                ├───RSV_A_R2.fastq.gz
                
          ├───📁Other_Pathogens
                ├───📁Pathogen_1
                ├───Pathogen_1_R1.fastq.gz         
                ├───Pathogen_1_R2.fastq.gz

                ├───📁Pathogen_2
                ├───Pathogen_2_R1.fastq.gz         
                ├───Pathogen_2_R2.fastq.gz

      ├───📁simulated_metagenomic_sample                        
           ├───📁SARS_COV_2
                ├───SARS_COV_2_R1.fastq.gz    
                ├───SARS_COV_2_R2.fastq.gz

           ├───📁RSV_A
                ├───RSV_A_R1.fastq.gz         
                ├───RSV_A_R2.fastq.gz
     
```

### <a name="arg"> Arguments

metaWEPP requries the following arguments, either through the 'config/config.yaml' or as command-line arguments.:

1. `KRAKEN_DB` - Name of the Kraken database.
2. `DIR` - Folder containing either metagenomic reads for analysis or reference FASTA files for simulating reads with MeSS. 
3. `SIMULATION_TOOL` - Use `MESS` to simulate reads, or leave it blank to analyze given reads.
4. `COVERAGE` - MESS's genomic coverage - Learn more about MESS's coverage calculation [here](https://metagenlab.github.io/MeSS/guide/simulate/coverage/).
5. `SEQUENCING_TYPE` - Sequencing read type (s:Illumina single-ended, d:Illumina double-ended, or n:ONT long reads)

Arguments required to run WEPP. These may vary by pathogens and can be placed within each pathogen-specific directory inside `pathogens_for_wepp`, e.g. `data/pathogens_for_wepp/<pathogen_name>/config.yaml`. See [Data Organization](#data).

6. `PRIMER_BED` - BED file for primers, which should be present in the `WEPP/primers` directory.
7. `MIN_AF` - Alleles with an allele frequency below this threshold in the reads will be masked.
8. `MIN_Q` - Alleles with a Phred score below this threshold in the reads will be masked.
9. `MAX_READS` - Maximum number of reads considered by WEPP from the sample. Helpful in the reducing runtime.
10. `CLADE_LIST` - List the clade annotation schemes stored in the MAT. SARS-CoV-2 MAT uses both nextstrain and pango lineage naming systems, so use "nextstrain,pango" for it.
11. `CLADE_IDX` - Index used for assigning clades to selected haplotypes from MAT. Use '1' for Pango naming and '0' for Nextstrain naming for SARS-CoV-2. Other pathogens usually follow a single lineage annotation system, so work with '0'. In case of NO lineage annotations, use '-1'. Lineage Annotations could be checked by running: "matUtils summary -i {TREE} -C {FILENAME}" -> Use '0' for annotation_1 and '1' for annotation_2.


⚠️ Example of pathogen specific `config.yaml` can be found in the [Quick Start](#example).


### <a name="run"> Run Command

metaWEPP requires `KRAKEN_DB` and `DIR` to be specified as command-line config arguments. Other parameters can also be provided via the config file. It also accepts `--cores`  to control the number of threads used during execution, and uses `--resources mess_slots=1` to  ensure the MeSS pipeline runs serially, preventing concurrency-related issues.

Examples:
1. Using all parameters from the config file:
```
snakemake --config KRAKEN_DB=test_kraken_DB DIR=simulated_metagenomic_sample --resources mess_slots=1 --cores 32
```

2. Overriding `CLADE_IDX` and `SIMULATE_TOOL` through the command line:
```
snakemake --config KRAKEN_DB=test_kraken_DB DIR=simulated_metagenomic_sample CLADE_IDX=1 SIMULATE_TOOL=MESS --resources mess_slots=1 --cores 32
```

### <a name="mat"> MAT Download
Mutation-annotated trees (MAT) for different pathogens are maintained by the UShER team, which can be found [here](https://dev.usher.bio). You can also create your own MAT for any pathogen from the consensus genome assemblies using [viral_usher](https://github.com/AngieHinrichs/viral_usher).

##  <a name="build-database"></a> Building Kraken Database

**Step 1:** Install the taxonomy. This is necessary for building a Kraken database. Replace "$DBNAME" above with your preferred database name.
```
kraken2-build --download-taxonomy --db $DBNAME
```

**Step 2:** Add sequences to the database's genome library.
```
kraken2-build --add-to-library <reference_genomes.fa> --db $DBNAME
```

**Step 3:** Build the database. 
```
kraken2-build --build --db $DBNAME
```

You can also customize kmer lengths with `--kmer-len` and `--minimizer-len`. For example,  
```
kraken2-build --build --db $DBNAME --kmer-len 21 --minimizer-len 15
```

⚠️ If you would like to save disk memory, perform the following commands:
```
rm -rf test_kraken_DB/taxonomy 
rm -rf test_kraken_DB/library  
```
More information about creating custom databases can be found [here](https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases).
