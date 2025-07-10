# metaWEPP: Metagenomic Wastewater-Based Epidemiology using Phylogenetic Placements

## Table of Contents
- [Introduction](#intro)
- [Installation](#install)
  - [Option-1: Install via Shell Commands](#shell)
- [Quick Start](#example)
  - [Example-1: Simulated Data](#mess)
  - [Example-2: Real World Data](#real-world) 
- [User Guide](#guide)
  - [Data Organization](#data)
  - [Arguments](#arg)
  - [Run Command](#run) 
- [Building Kraken Database](#build-database)

<br>


## <a name="intro"></a> Introduction

metaWEPP is a Snakemake-based bioinformatics pipeline designed to enable rapid classification and haplotype-level analysis of mixed-pathogen metagenomic samples. Developed for flexible, high-throughput use in public health surveillance, metaWEPP integrates [Kraken2](https://github.com/DerrickWood/kraken2) for taxonomic classification and routes identified pathogen reads into [WEPP](https://github.com/TurakhiaLab/WEPP) for phylogenetic placement and haplotype inference. The pipeline automates the end-to-end workflow—from raw mixed reads to lineage-level analysis—with optional support for simulated read generation using [MeSS](https://github.com/metagenlab/MeSS). 

<div align="center">
    <img src="docs/images/metawepp-figure.png" width="600">
    <div><b>Figure 1: metaWEPP Pipeline Visual</b></div>
</div>


## <a name="install"></a> Installation

### <a name="shell"></a> Option-1: Install via Shell Commands.

**Step 1:** Clone the repository.
```
git clone https://github.com/TurakhiaLab/metaWEPP.git
cd metaWEPP
```

**Step 2:** Install Kraken.
The following commands install kraken and also update the `$PATH` variable for easily running the tool.
```
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh .
echo -e "\nexport PATH=\"$(pwd):\$PATH\"" >> ~/.bashrc
source ~/.bashrc
cd ..
```
**Step 4:** Install WEPP.

```
git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git
```
View the WEPP installation guid starting from option 3 on the [WEPP repo](https://github.com/TurakhiaLab/WEPP/tree/main?tab=readme-ov-file#-option-3-install-via-shell-commands-requires-sudo-access).

**Step 5:** Install MeSS.

Follow the MeSS installation guide on the [MeSS Quick Start](https://github.com/metagenlab/MeSS?tab=readme-ov-file#zap-quick-start).

---

##  <a name="example"></a> Quick Start

### <a name="mess"></a> Example - 1 SARS-CoV-2 Dataset: Run the pipeline with MeSS simulated data

**Step 1:** Download the metagenomic sample along with the MAT and reference files for SARS-CoV-2.
```
wget https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2021/12/05/public-2021-12-05.all.masked.pb.gz

mkdir -p data/pathogens_for_wepp/sars_cov_2
mkdir data/simulated_metagenomic_sample

mv public-2021-12-05.all.masked.pb.gz data/pathogens_for_wepp/sars_cov_2
cp metagenomic_references/NC_045512.2.fasta data/pathogens_for_wepp/sars_cov_2
cp metagenomic_example.fa data/simulated_metagenomic_sample
```
**Step 2:** Prepare the config.yaml for SARS-CoV-2.
```
cat <<EOF > data/pathogens_for_wepp/sars_cov_2/config.yaml
PRIMER_BED: none.bed
CLADE_IDX: 1
SEQUENCING_TYPE: d
CLADE_LIST: nextstrain,pango
EOF
```

**Step 3:** Build the Kraken database.
```
mkdir test_kraken_DB
kraken2-build --download-taxonomy --db test_kraken_DB
k2 add-to-library --db test_kraken_DB --file metagenomic_references/*
kraken2-build --build --db test_kraken_DB
```
⚠️ Note that you must add the reference genome (in this example, the SARS-COV-2 reference genome file) into the custom database for the pipeline to work. This is done for us in the third line of Step 3.


**Step 4:**  Run the pipeline.
```
snakemake --config DIR=simulated_metagenomic_sample SIMULATION_TOOL=MESS KRAKEN_DB=test_kraken_DB --resources mess_slots=1 --cores 32
```

**Step 5:**  Analyze Results.

The classification distribuction can be found in `results/simulated_metagenomic_sample/classification_proportions.png`. All WEPP results can be found in the `WEPP/results/sars_cov_2_simulated_metagenomic_sample` directory. 

### <a name="real-world"></a> Example - 2: Real World Data

This example will take our own metagenomic wastewater reads and use them as input FQ file for the pipeline running on a RSV-A dataset.

**Step 1:** Download the metagenomic sample along with the MAT and reference files for RSV-A.
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

**Step 2:** Prepare the config.yaml for RSV-A.
```
cat <<EOF > data/pathogens_for_wepp/rsv_a/config.yaml
PRIMER_BED: RSVA_all_primers_best_hits.bed
CLADE_IDX: 0
SEQUENCING_TYPE: d
CLADE_LIST: annotation_1
EOF
```

**Step 3:** Build the Kraken database.
```
mkdir test_kraken_DB
kraken2-build --download-taxonomy --db test_kraken_DB
k2 add-to-library --db test_kraken_DB --file metagenomic_references/* 
kraken2-build --build --db test_kraken_DB
```
⚠️ Note that you must add the reference genome (in this example, the RSV-A reference genome file) into the custom database for the pipeline to work. This is done for us in the third line of Step 3.

**Step 4:**  Run the pipeline.
```
snakemake --config DIR=real_metagenomic_sample KRAKEN_DB=test_kraken_DB --resources mess_slots=1 --cores 32
```

**Step 5:**  Analyze Results.

The classification distribuction can be found in `results/real_metagenomic_sample/classification_proportions.png`. WEPP results can be found in the `WEPP/results/rsva_a_real_metagenomic_sample` directory. 


## <a name="guide"></a> User Guide

### <a name="data"> Data Organization
We assume that all wastewater samples are organized in the data directory, each within its own subdirectory given by DIR argument (see Run Command). For each pathogen to be analyzed for variants with WEPP, place its reference genome, corresponding MAT, and the config.yaml (optional) in the folder: `data/pathogens_for_wepp/<pathogen_name>`. Moreover, each created `DIR` inside `data` is expected to contain either of these files.
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
               ├───config.yaml                   # customized WEPP config for SARS-COV-2

          ├───📁RSV_A                   
               ├───RSV_A_ref_genome.fa     
               ├───RSV_A_mat.pb.gz 
               ├───config.yaml                   # customized WEPP config for RSV-A

     ├───📁real_metagenomic_sample               # [User Created] Contains metagenomic sample for analysis
          ├───metagenomic_reads_R1.fastq.gz      
          ├───metagenomic_reads_R2.fastq.gz

     ├───📁simulated_metagenomic_sample          # [User Created] Contains metagenomic fasta for simulating reads with MeSS               
          ├───metagenomic_reference.fa           
          ├───MeSS_R1.fastq.gz                   # [MeSS Generated] Simulated reads
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

metaWEPP requries the following arguments for config/config.yaml:

1. `KRAKEN_DB` - Name of the Kraken database.
2. `DIR` - Folder containing either metagenomic reads for analysis or reference FASTA files for read simulation using MeSS. 
3. `SIMULATION_TOOL` - Input `MESS` to simulate reads with MeSS, or leave it blank to analyze your own reads.
4. `COVERAGE` - MESS's genomic coverage - Learn more about MESS's coverage calculation [here](https://metagenlab.github.io/MeSS/guide/simulate/coverage/).
5. `SEQUENCING_TYPE` - Sequencing read type (s:Illumina single-ended, d:Illumina double-ended, or n:ONT long reads)

Arguments required to run WEPP. These may vary by pathogen and can be placed within each pathogen-specific directory inside `pathogens_for_wepp`, e.g. `data/pathogens_for_wepp/<pathogen_name>/config.yaml`. See [Data Organization](#data).

6. `PRIMER_BED` - BED file for primers. These should be present in the `WEPP/primers` directory.
7. `MIN_AF` - Alleles with an allele frequency below this threshold in the reads will be masked.
8. `MIN_Q` - Alleles with a Phred score below this threshold in the reads will be masked.
9. `MAX_READS` - Maximum number of reads considered by WEPP from the sample. Helpful for reducing runtime
10. `CLADE_LIST` - List the clade annotation schemes used in the MAT. SARS-CoV-2 MAT uses both nextstrain and pango lineage naming systems, so use "nextstrain,pango" for it.
11. `CLADE_IDX` - Index used for assigning clades to selected haplotypes from MAT. Use '1' for Pango naming and '0' for Nextstrain naming for SARS-CoV-2. Other pathogens usually follow a single lineage annotation system, so work with '0'. In case of NO lineage annotations, use '-1'. Lineage Annotations could be checked by running: "matUtils summary -i {TREE} -C {FILENAME}" -> Use '0' for annotation_1 and '1' for annotation_2.


⚠️ Example of pathogen specific `config.yaml` can be found in [Quick Start](#example).


### <a name="run"> Run Command

metaWEPP requires `KRAKEN_DB` and `DIR` to be specified as command-line config arguments. Other parameters can be provided via the config file. It also accepts `--cores`  to control the number of threads used during execution, and uses `--resources mess_slots=1` to  ensure the MeSS pipeline runs serially, preventing concurrency-related issues.

Examples:
1. Using all parameters from the config file:
```
snakemake --config KRAKEN_DB=test_kraken_DB DIR=simulated_metagenomic_sample --resources mess_slots=1 --cores 32
```

2. Overriding `CLADE_IDX` and `SIMULATE_TOOL` through the command line:
```
snakemake --config KRAKEN_DB=test_kraken_DB DIR=simulated_metagenomic_sample CLADE_IDX=1 SIMULATE_TOOL=MESS --resources mess_slots=1 --cores 32
```


##  <a name="build-database"></a> Building Kraken Database
If you would like more information on building a Kraken database, see below:

**Step 1:** Install the taxonomy. This is necessary for building a Kraken database. Replace "$DBNAME" above with your preferred database name.
```
kraken2-build --download-taxonomy --db $DBNAME
```

**Step 2:** Add sequences to the database's genomic library.
```
kraken2-build --add-to-library <reference_genomes.fa> --db $DBNAME
```

**Step 3:** Build the database 
```
kraken2-build --build --db $DBNAME
```

You can also customize kmer with `--kmer-len` and `--minimizer-len` option if needed. For example,  
```
kraken2-build --build --db $DBNAME --kmer-len 21 --minimizer-len 15
```

⚠️ If you would like to save disk memory, perform the following commands:
```
rm -rf test_kraken_DB/taxonomy 
rm -rf test_kraken_DB/library  
```
More information about creating custom databases can be found [here](https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases).
