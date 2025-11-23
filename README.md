<div align="center">

# metaWEPP: Improving the resolution of metagenomic analysis using WEPP

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
  - [Option-2: Install via Shell Commands (requires sudo access)](#shell)
- [Quick Start](#example)
  - [Example-1: Simulated metagenomic sample](#simulate)
  - [Example-2: Real world metagenomic sample](#real-world) 
- [User Guide](#guide)
  - [Data Organization](#data)
  - [metaWEPP Arguments](#arg)
    - [Detailed Example](#argexample)
  - [Run Command](#run) 
  - [MAT Download](#mat) 
- [Building Kraken Database](#build-database)
  - [Downloading prebuilt database](#prebuilt)
  - [Creating custom database](#custom)

<br>


## <a name="intro"></a> Introduction

metaWEPP is a Snakemake-based bioinformatics pipeline that achieves near-haplotype resolution in metagenomic analysis. As illustrated in the figure, metaWEPP can analyze metagenomic or mixed-genome samples from environmental sources and clinical specimens. metaWEPP first uses standard taxonomic classifiers to assign sequencing reads to known species, then applies [WEPP](https://github.com/TurakhiaLab/WEPP) to phylogenetically place these reads onto updated, species-specific mutation-annotated trees built from all publicly available clinical sequences,  and finally selects the subset of haplotypes that best explains the sample. It also reports unaccounted alleles that are indicative of novel variants and includes an interactive dashboard to provide a detailed read-level visualization for each species. 

<div align="center">
    <img src="metaWEPP_Overview.png" width="1000">
    <div><b>Figure 1: metaWEPP Overview</b></div>
</div>


## <a name="install"></a> Installation
Currently, adding new pathogens for haplotype-level analysis in WEPP depends on `viral_usher`, which requires Docker. Therefore, we recommend installing WEPP using the provided Dockerfile.

### <a name="docker"></a> Option-1: Install via Dockerfile.

**Step 1:** Clone the metaWEPP repository.
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
# -p <host_port>:<container_port> → Maps container port to a port on your host (Accessing Dashboard, NOT needed otherwise)
# Replace <host_port> with your desired local port (e.g., 100 or 8080)
docker run -it -p 80:80 metawepp
```

All set to try the [examples](#example).


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

**Step 4:** Install `Minimap2`, `viral_usher`, `matplotlib`, and `snakemake`.

```
sudo apt-get install minimap2
pip install viral_usher matplotlib snakemake
```
**Step 4:** Install `WEPP`.

```
git clone --recurse-submodules https://github.com/TurakhiaLab/WEPP.git
```
View the WEPP installation guide starting from Option 3 in the [WEPP repository](https://github.com/TurakhiaLab/WEPP/tree/main?tab=readme-ov-file#-option-3-install-via-shell-commands-requires-sudo-access).

All set to try the [examples](#example).


##  <a name="example"></a> Quick Start

### <a name="simulate"></a> Example - 1: Simulated metagenomic sample 
**Step 1:** Download the MAT file for SARS-CoV-2 and RSV-A.
```
wget https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2021/12/05/public-2021-12-05.all.masked.pb.gz
wget https://hgdownload.gi.ucsc.edu/hubs/GCF/002/815/475/GCF_002815475.1/UShER_RSV-A/2025/04/25/rsvA.2025-04-25.pb.gz
```

**Step 2:** Build a custom Kraken database for taxonomic classification.
```
mkdir test_kraken_DB
kraken2-build --download-taxonomy --db test_kraken_DB
kraken2-build --download-library viral --db test_kraken_DB
kraken2-build --build --db test_kraken_DB
```

**Step 3:** Confirm the `config.yaml` for pathogens analysis.
```
PATHOGENS: "default,sars_cov_2"

CLADE_LIST: "pango,nextstrain:pango"

CLADE_IDX: "0,1"

MIN_DEPTH_FOR_WEPP: "10"
```

**Step 4:** Run the pipeline

Follow the command prompts to add the pathogens.

```
snakemake --config KRAKEN_DB=test_kraken_DB DIR=simulated_metagenomic_sample --cores 8
```

For SARS-CoV-2:
```txt
a) Type "sars cov 2" and press [enter] to search for SARS-CoV-2.
b) Select Severe acute respiratory syndrome coronavirus 2 by pressing the corresponding number, then [enter].
c) Select NC_045512.2 by pressing the corresponding number, then [enter].
d) Provide the path of the MAT "./public-2021-12-05.all.masked.pb.gz".
```
For RSV-A:
```txt
a) Type "rsv a" and press [enter] to search for RSV.
b) Select human respiratory syncytial virus by pressing the corresponding number, then [enter].
c) Select NC_038235.1
d) Provide the path of the MAT "./rsvA.2025-04-25.pb.gz".
e) 'N' to end adding pathogen.
```

**Step 5:**  Analyze Results.

Taxonomic classification results can be viewed at `results/simulated_metagenomic_sample/classification_proportions.png`. All results generated by WEPP for SARS-CoV-2 variant analysis are present in `WEPP/results/sars_cov_2_simulated_metagenomic_sample`, and it's the same for RSV-A.

Expected metagenomic classification proportions:

- 66.33%  Severe acute respiratory syndrome coronavirus 2
- 31.38%  Human respiratory syncytial virus A
- 1.94%   human respiratory syncytial virus
- 0.35%   Unclassified


### <a name="real-world"></a> Example - 2: Real world metagenomic sample
**Step 1:** Download the real world metagenomic sample.
```
prefetch -v SRR14530845
fastq-dump --split-e SRR14530845
```

**Step 2:** Download the MAT file for SARS-CoV-2.
```
wget https://hgdownload.gi.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2021/12/05/public-2021-12-05.all.masked.pb.gz
```

**Step 3:** Build a custom Kraken database for taxonomic classification.
```
mkdir test_kraken_DB
kraken2-build --download-taxonomy --db test_kraken_DB
kraken2-build --download-library viral --db test_kraken_DB
kraken2-build --build --db test_kraken_DB
```

**Step 4:** Confirm the `config.yaml` for pathogens analysis.
```
PATHOGENS: "default,sars_cov_2"

CLADE_LIST: "pango,nextstrain:pango"

CLADE_IDX: "0,1"

MIN_DEPTH_FOR_WEPP: "10"
```

**Step 5:** Run the pipeline 

```
snakemake --config KRAKEN_DB=test_kraken_DB DIR=simulated_metagenomic_sample --cores 8
```

Follow the command prompts to add the pathogens.

For SARS-CoV-2:
```txt
a) Type "sars cov 2" and press [enter] to search for SARS-CoV-2.
b) Select Severe acute respiratory syndrome coronavirus 2 by pressing the corresponding number, then [enter].
c) Select  by pressing the corresponding number, then [enter].
d) Provide the path of the MAT "./public-2021-12-05.all.masked.pb.gz".
```

For Tomato brown rugose fruit virus:
```txt
a) Type "tomato brown rugose fruit virus" and press [enter] to search for RSV.
b) Select Tomato brown rugose fruit virus by pressing the corresponding number, then [enter].
c) Select NC_028478.1
d) [enter] to build the MAT in viral usher. 
e) Follow the viral usher instruction, [enter], [enter], [enter], [enter].
f) Input "./data/pathogens_for_wepp/tomato_brown_rugose/NC_028478.1.fna" for fasta file.
g) Input "./data/pathogens_for_wepp/tomato_brown_rugose/viral_usher_build" for MAT target directory.
h) 'N' to end adding pathogen.
```

**Step 6:**  Analyze Results.

Taxonomic classification results can be viewed at `results/real_metagenomic_sample/classification_proportions.png`, while all results generated by WEPP for SARS-CoV-2 variant analysis are present in `WEPP/results/sars_cov_2_simulated_metagenomic_sample`, and it's the same for tomato rugose virus.

Expected metagenomic classification proportions:

- 9.76%   Tomato brown rugose fruit virus
- 21.71%  Severe acute respiratory syndrome coronavirus 2
- 41.46%  Unclassified


## <a name="guide"></a> User Guide

### <a name="data"> Data Organization
We assume that all metagenomic samples are stored in the `data` directory, each within its own subdirectory given by DIR argument (see [Run Command](#run)). Pathogen species to be analyzed with WEPP at the haplotype level can be added either manually or automatically fetched by metaWEPP based on user inputs. For manual addition, create a directory for each species at: 
```
data/pathogens_for_wepp/<pathogen_species_name>
```
and place both its reference genome and the MAT file inside. Additionally, add the `Taxonomy ID` fof each added species to `data/pathogens_for_wepp/added_taxons.csv` on a new line using the format: 
```
<Taxonomy ID>,<pathogen_species_name>
``` 

Each created `DIR` inside `data` is expected to contain only the metagenomic sequencing reads, with filenames ending in `*_R{1/2}.fastq.gz` for paired-end reads, and `*.fastq.gz` for single-end reads. For each metagenomic sample, metaWEPP generates species-level results in the corresponding sample subdirectories under `results`. Haplotype-level results for each pathogen in a sample are located within the respective pathogen directories under `WEPP/results` as shown below. 

Visualization of metaWEPP's workflow directories
```
📁 metaWEPP
└───📁config
     ├───📁config.yaml                           # configuration for metaWEPP
└───📁data                                    
     ├───📁pathogens_for_wepp                    # [metaWEPP Generated/User Provided] Pathogens for haplotype analysis
          ├───📁RSV_A                   
               ├───RSV_A_ref_genome.fa     
               ├───RSV_A_mat.pb.gz
          ├───📁Tomato_Brown_Rugose_Fruit_Virus
               ├───Rugose_Virus_ref_genome.fa     
               ├───Rugose_Virus_mat.pb.gz
          ├───📁SARS_COV_2                
               ├───SARS_COV_2_ref_genome.fa     
               ├───SARS_COV_2_mat.pb.gz 

     ├───📁simulated_metagenomic_sample          # [User Provided] Contains metagenomic reads
          ├───sample_R1.fastq.gz      
          ├───sample_R2.fastq.gz
     ├───📁real_metagenomic_sample               # [User Provided] Contains metagenomic reads
          ├───sample_R1.fastq.gz      
          ├───sample_R2.fastq.gz


└───📁results                                    # [metaWEPP Generated] 
      ├───📁real_metagenomic_sample
           ├───📁Tomato_Brown_Rugose_Fruit_Virus
                ├───Tomato_Brown_Rugose_Fruit_Virus_R1.fastq.gz         
                ├───Tomato_Brown_Rugose_Fruit_Virus_R2.fastq.gz             
           ├───📁SARS_COV_2
                ├───SARS_COV_2_R1.fastq.gz    
                ├───SARS_COV_2_R2.fastq.gz
                
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
     
└───📁WEPP                                    
      ├───📁results                              # [metaWEPP Generated]
           ├───📁Tomato_Brown_Rugose_Fruit_Virus_real_metagenomic_sample
                ├───file_1  
                ├───file_2
           ├───📁SARS_COV_2_real_metagenomic_sample
                ├───file_1  
                ├───file_2
           ├───📁SARS_COV_2_simulated_metagenomic_sample
                ├───file_1  
                ├───file_2
           ├───📁RSV_A_simulated_metagenomic_sample
                ├───file_1  
                ├───file_2
```

### <a name="arg"> metaWEPP Arguments

The metaWEPP Snakemake pipeline requires the following arguments, which can be provided either via the configuration file (`config/config.yaml`) or passed directly on the command line using the `--config` argument. The command line arguments take precedence over the config file.

1. `DIR` - Folder containing the metagenomic reads.
2. `KRAKEN_DB` - Folder containing the Kraken2 database. 
3. `SEQUENCING_TYPE` - Sequencing read type (s:Illumina single-ended, d:Illumina double-ended, or n:ONT long reads)
4. `PRIMER_BED` - BED file for primers, which should be present in the `WEPP/primers` directory.
5. `MIN_AF` - Alleles with an allele frequency below this threshold in the reads will be masked (Illumina: 0.5%, Ion Torrent: 1.5%, ONT: 2%) by WEPP.
6. `MIN_DEPTH` - Sites with read depth below this threshold will be masked by WEPP.
7. `MIN_Q` - Alleles with a Phred score below this threshold in the reads will be masked by WEPP.
8. `MIN_PROP` -  Minimum Proportion of haplotypes detected by WEPP (Wastewater Samples: 0.5%, Clinical Samples: 5%).
9. `MIN_LEN` -  Minimum read length to be considered after ivar trim (Default: 80).
10. `MAX_READS` - Maximum number of reads considered by WEPP from the sample. Helpful for reducing runtime.
11. `DASHBOARD_ENABLED` - Enables WEPP dashboard for visualization of haplotype results
12. `PATHOGENS` - List of pathogens with custom WEPP settings. Any species not listed here will use the `default` settings.
12. `CLADE_LIST` - Comma-separated list of clade annotation schemes present in the MAT file. Each element corresponds to the pathogen species in the order specified in `PATHOGENS`. If there is no clade annotation for a pathogen species, do not provide any value for that species.
13. `CLADE_IDX` - Comma-separated list of clade indices for each pathogen. If a pathogen has no lineage annotations, use `-1`. Each element corresponds to the pathogen species in the order specified by `PATHOGENS`.
14. `MIN_DEPTH_FOR_WEPP` - Minimum read coverage required to run WEPP for any pathogen species.
15. `MIN_PROP_FOR_WEPP` - Minimum relative abundance of a species before metaWEPP prompts to add it for haplotype-level analysis.


#### <a name="argexample"> Detailed Example:
```
PATHOGENS: "default,sars_cov_2"
CLADE_LIST: "pango,nextstrain:pango"
CLADE_IDX: "0,1"
MIN_DEPTH_FOR_WEPP: "10"
```
This means:

**default pathogens** (all pathogens not named explicitly):

Use the first element of the list in CLADE_IDX "0" and CLADE_LIST "pango"

Use depth threshold 10 (only one value provided -> applies to all pathogens)

**sars_cov_2**:

Use the second element of the list in CLADE_IDX "1" and CLADE_LIST "nextstrain:pango"(will be convert to "nextstrain,pango" when running WEPP)

Use depth threshold 10 (only one value provided -> applies to all pathogens)

| Pathogen   | Clade Index | Clade List        | Minimum depth |
| ---------- | ---------   | ----------------- | --------------|
| default    | 0           | pango             | 10            |
| sars_cov_2 | 1           | nextstrain:pango  | 10            |


### <a name="run"> Run Command

metaWEPP requires `KRAKEN_DB` and `DIR` to be specified as command-line config arguments. Other parameters can also be provided via the config file. It also accepts `--cores`  to control the number of threads used during execution.

Examples:
1. Using all parameters from the config file:
```
snakemake --config KRAKEN_DB=<path to the Kraken database> DIR=<name of the sample folder for analysis> --cores 32
```


### <a name="mat"> MAT Download
Mutation-annotated trees (MAT) for different pathogens are maintained by the UShER team, which can be found [here](https://dev.usher.bio). You can also create your own MAT for any pathogen from the consensus genome assemblies using [viral_usher](https://github.com/AngieHinrichs/viral_usher).

##  <a name="build-database"></a> Building Kraken Database
Users can either download prebuilt databases available online or create custom ones using their genomes.

### <a name="prebuilt"> Downloading prebuilt database
**Step 1:** Get the link of `.tar.gz` file of the genome collection you want from [here](https://benlangmead.github.io/aws-indexes/k2). Download the database using either `wget` or `curl` command.


**Step 2:** Unzip the downloaded database in a new directory using `tar`.
```
tar -xzf file.tar.gz
```


### <a name="custom"> Creating custom database
**Step 1:** Install the taxonomy. This is necessary for building a Kraken database. Replace "$DBNAME" above with your preferred database name.
```
kraken2-build --download-taxonomy --db $DBNAME
```

**Step 2:** Add sequences to the database's genome library.
```
k2 add-to-library --db $DBNAME --file <reference_genomes.fa> 
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
