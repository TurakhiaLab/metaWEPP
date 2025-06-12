# Wastewater-Based Epidemiology using Phylogenetic Placements: Classification Pipeline

## Table of Contents
- [Introduction](#intro)
- [Installation](#install)
  - [Option-1: Install via Conda](#dockerhub)
- [Quick Start](#example)

<br>


## <a name="intro"></a> Introduction

From [WEPP Main Repo](https://github.com/TurakhiaLab/WEPP):
> WEPP (**W**astewater-Based **E**pidemiology using **P**hylogenetic **P**lacements) is a phylogeny-based pipeline that estimates haplotype proportions from wastewater sequencing reads using a mutation-annotated tree (MAT) (Figure 1A). By improving the resolution of pathogen variant detection, WEPP enables critical epidemiological applications previously feasible only through clinical sequencing. It also flags potential novel variants via unaccounted mutations, which can be examined at the read level using the interactive dashboard (Figure 1C).
>
> WEPP’s algorithm begins with parsimonious placement of all reads onto the MAT, followed by identifying candidate haplotype nodes, or “Peaks” (Figure 1B). This set is expanded with neighboring haplotypes of selected Peaks to form a candidate pool, which is passed to a deconvolution algorithm to estimate haplotype abundances. This pool is iteratively refined by retaining haplotypes above a threshold and adding their neighbors until convergence.


<div align="center">
    <img src="images/WEPP_Overview.png" width="600">
    <div><b>Figure 1: Overview of WEPP</b></div>
</div>

This repo serves as an extension to WEPP by utilizing [Kraken2](https://github.com/DerrickWood/kraken2) to classify genomes that could then be inputted into WEPP. The primary addition is a snakemake workflow that takes in mixed unclassified reads, inputs them into Kraken2, separates the classified reads into its respected pathogens, then inputting the desired reads into WEPP. Some additional features include read simulation through [art_illumina](https://manpages.debian.org/testing/art-nextgen-simulation-tools/art_illumina.1.en.html) and [MeSS](https://github.com/metagenlab/MeSS).


## <a name="install"></a> Installation
This extension to WEPP currently only supports one installation method through `conda install`.

### <a name="dockerhub"></a> Option-1: Install via Conda.


**Step 0:**  Build a Kraken database
1. Install a taxonomy. 
```
kraken2-build --download-taxonomy --db $DBNAME
```
(Replace "$DBNAME" above with your preferred database name/location. The database will use approximately 100 GB of disk space during creation. )

2.
(If needed) Install one or more reference libraries: https://github.com/DerrickWood/kraken2/wiki/Manual#standard-kraken-2-database
```
kraken2-build --download-library bacteria --db $DBNAME
```

Add sequence to the database's genomic library using the --add-to-library switch, e.g.:
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
3. Build the database 
```
kraken2-build --build --db $DBNAME
```
Customized kmer with `--kmer-len` and `--minimizer-len` option if needed.

More information in https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases

**Step 1:** Clone the extension repository.
```bash
git clone https://github.com/TurakhiaLab/metagenomic-WBE.git
```
**Step 2:** Clone the WEPP repository in the extension repository.
```bash
cd metagenomic-WBE \
git clone https://github.com/TurakhiaLab/WEPP.git
```
**Step 3:** Install the conda environment package.
```bash
conda env create -f environment.yaml
conda activate wepp_ext
```


### <a name="arguments"></a> WEPP Arguments

From the WEPP repository:
> The WEPP Snakemake pipeline requires the following arguments, which can be provided either via the configuration file (`config/config.yaml`) or passed directly on the command line using the `--config` argument. The command line arguments take precedence over the config file.
> 1. `DIR` - Folder name containing the wastewater reads
> 2. `FILE_PREFIX` - File Prefix for all intermediate files 
> 3. `REF` - Reference Genome in fasta
> 4. `TREE` - Mutation-Annotated Tree
> 5. `SEQUENCING_TYPE` - Sequencing read type (s:Illumina single-ended, d:Illumina double-ended, or n:ONT long reads)
> 6. `PRIMER_BED` - BED file for primers from the `primers` folder
> 7. `MIN_AF` - Alleles with an allele frequency below this threshold in the reads will be masked. 
> 8. `MIN_Q` - Alleles with a Phred score below this threshold in the reads will be masked.
> 9. `MAX_READS` - Maximum number of reads considered by WEPP from the sample. Helpful for reducing runtime
> 10. `CLADE_IDX` - Index used for assigning clades to selected haplotypes from MAT. Generally '1' for SARS-CoV-2 MATs and '0' for others. > Could be checked by running: "matUtils summary -i {TREE} -C {FILENAME}" -> Use '0' for annotation_1 and '1' for annotation_2 

### <a name="example"></a> Usage
The entire pipeline can be ran through the following command:
```
snakemake \
  --snakefile Snakefile \
  --cores 32 \
  --configfile config.yaml
```
This will run the full pipeline and run WEPP for each taxid in `target_taxids` in the `config.yaml` file (in the wepp extension directory).

