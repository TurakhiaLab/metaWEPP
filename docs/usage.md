# User Guide <a name="guide"></a>

## <b>Data Organization</b> <a name="data"></a>

We assume that all metagenomic samples are stored in the `data` directory, each within its own subdirectory given by `DIR` argument (see [Run Command](#run)). Pathogen species to be analyzed with WEPP at the haplotype level can be added either manually or automatically through the metaWEPP workflow based on user input during the analysis. For manual addition, create a separate directory for each species and place the corresponding reference genome FASTA file and mutation-annotated tree (MAT) file in that directory. The files should be placed under:
```text
data/pathogens_for_wepp/<pathogen_species_name>
```

After adding a new species manually, delete the file `data/pathogens_for_wepp/added_taxons.csv` to ensure that metaWEPP refreshes the list of available species on the next run.

Each created `DIR` inside `data` is expected to contain only the metagenomic sequencing reads, with filenames ending in `*_R{1/2}.fastq.gz` for paired-end reads, and `*.fastq.gz` for single-end reads. For each metagenomic sample, metaWEPP generates species-level results in the corresponding sample subdirectories under `results`. Haplotype-level results for each pathogen in a sample are located within the respective pathogen directories under `WEPP/results` as shown below.

Visualization of metaWEPP's workflow directories:
```text
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

## <b>metaWEPP Arguments</b> <a name="arguments"></a>

The metaWEPP Snakemake pipeline requires the following arguments, which can be provided either via the configuration file (`config/config.yaml`) or passed directly on the command line using the `--config` argument. The command line arguments take precedence over the config file.

1. `DIR` - Folder(s) containing metagenomic reads. Multiple samples can be analyzed in parallel by specifying a comma-separated list of folders.
2. `KRAKEN_DB` - Either (a) the path to an existing Kraken2 database folder, or (b) an HTTP(S) URL to a `.tar.gz` archive (e.g. `https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20251015.tar.gz`). When a URL is given, the archive is downloaded into a directory named after the tarball stem on first run and reused on subsequent runs.
3. `SEQUENCING_TYPE` - Sequencing read type (s:Illumina single-ended, d:Illumina double-ended, or n:ONT long reads)
4. `MIN_AF` - Alleles with an allele frequency below this threshold in the reads will be masked (Illumina: 0.5%, Ion Torrent: 1.5%, ONT: 2%) by WEPP.
5. `MIN_DEPTH` - Sites with read depth below this threshold will be masked by WEPP.
6. `MIN_Q` - Alleles with a Phred score below this threshold in the reads will be masked by WEPP.
7. `MIN_PROP` - Minimum Proportion of haplotypes detected by WEPP (Wastewater Samples: 0.5%, Clinical Samples: 5%).
8. `MIN_LEN` - Minimum read length to be considered after ivar trim (Default: 80).
9. `MAX_READS` - Maximum number of reads considered by WEPP from the sample. Helpful for reducing runtime.
10. `DASHBOARD_ENABLED` - Enables WEPP dashboard for visualization of haplotype results.
11. `ADD_SPECIES_RUNTIME` - Asks users to add pathogen species at runtime when enabled.
12. `PATHOGENS` - List of pathogens with custom WEPP settings. Any species not listed here will use the `default` settings.
13. `PRIMER_BED` - Absolute path(s) to the BED file(s) for primers. Comma-separated paths ordered to match `PATHOGENS` for per-species primer trimming. Leave a slot blank for species without primer trimming.
14. `CLADE_LIST` - Clade annotation schemes in the MAT file. Either a single value broadcast to every species, or use comma-separated values ordered to match `PATHOGENS`. Leave a slot blank for species without clade annotations.
15. `CLADE_IDX` - Clade indices for each pathogen. Either a single value broadcast to every species, or use comma-separated values ordered to match `PATHOGENS`. Use `-1` for species without lineage annotations.
16. `MIN_DEPTH_FOR_WEPP` - Minimum read coverage required to run WEPP for any pathogen species.
17. `MIN_PROP_FOR_WEPP` - Minimum relative abundance of a species before metaWEPP prompts to add it for haplotype-level analysis.
18. `CORES_PER_PATHOGEN` - Cores allocated to each concurrent WEPP species run. Either a single value is broadcast to every species, or comma-separated values ordered to match `PATHOGENS` are used for per-species allocation. Empty values in per-species allocation list, or 'auto' evenly divides whatever threads remain from the `--cores` after subtracting the explicit allocations.

### <b>Example of species-specific arguments</b> <a name="argexample"></a>

```text
PATHOGENS: "default,SARS_COV_2,RSV_A"
CLADE_LIST: ",nextstrain:pango,nextstrain"
CLADE_IDX: "-1,1,0"
CORES_PER_PATHOGEN: ",20,8"
PRIMER_BED: ",/path_to_sars/sars.bed,/path_to_rsv/rsv.bed"
```

These arguments are applied by first checking whether the species listed in `PATHOGENS` are present in `data/pathogens_for_wepp`. Then, the corresponding `CLADE_LIST` and `CLADE_IDX` values are used for WEPP analysis as follows:

**SARS_COV_2** - `CLADE_LIST` is **nextstrain,pango**, and `CLADE_IDX` is **1**.

**RSV_A** - `CLADE_LIST` is **nextstrain** and `CLADE_IDX` is **0**.

**All other species** - `CLADE_LIST` is not passed, and `CLADE_IDX` is **-1**.

!!!Note
    ⚠️ When generating MAT with viral_usher, DO NOT use pre-built trees. These trees have been re-rooted and are constructed using reference genomes that differ from those selected within the metaWEPP workflow, which can lead to inconsistencies in the downstream analyses.

## <b>Run Command</b> <a name="run"></a>

metaWEPP requires `KRAKEN_DB` and `DIR` to be specified as command-line arguments. Other parameters can be provided via the config file. Additionally, the `--cores` argument must be supplied on the command line to specify the number of threads used by the workflow.

Examples:

1. Using all the parameters from the config file:
```bash
run-metawepp --config KRAKEN_DB=viral_kraken_db DIR=simulated_metagenomic_sample --cores 32
```

2. Overriding MIN_Q and MIN_PROP_FOR_WEPP through command line:
```bash
run-metawepp --config KRAKEN_DB=viral_kraken_db DIR=simulated_metagenomic_sample MIN_Q=25 MIN_PROP_FOR_WEPP=0.05 --cores 32
```

## <b>MAT for pathogen species</b> <a name="mat"></a>

Mutation-annotated trees (MATs) for some pathogen species are maintained by the UShER team and can be found [here](https://dev.usher.bio). Alternatively, metaWEPP includes an internal pipeline to generate MATs for any virus with sequences available on GenBank.

## <b>Kraken2 Database</b> <a name="build-database"></a>

Users can either download prebuilt databases available online or create custom databases using their own genome sequences.

### <b>Passing URL of database</b> <a name="url"></a>

For prebuilt databases, the easiest option is to pass the `.tar.gz` URL directly through the `KRAKEN_DB=` parameter when running metaWEPP. For example:
```bash
KRAKEN_DB=https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20251015.tar.gz
```

On the first run, metaWEPP will automatically download and extract the database into a local directory, and will reuse the downloaded database in subsequent runs. The manual `wget` and `tar` commands provided below achieve the same result while allowing you to choose a custom download location and directory name.

### <b>Downloading prebuilt database</b> <a name="prebuilt"></a>

**Step 1:** Get the link of `.tar.gz` file of the genome collection you want from [here](https://benlangmead.github.io/aws-indexes/k2). Download the database using either `wget` or `curl` command.

**Step 2:** Unzip the downloaded database in a new directory using `tar`.
```bash
mkdir -p <prebuilt_database>
tar -xvzf <database.tar.gz> -C <prebuilt_database>
```

### <b>Creating custom database</b> <a name="custom"></a>

**Step 1:** Download the taxonomy data, which is required to build a custom Kraken2 database.
```bash
kraken2-build --download-taxonomy --db <custom_database>
```

**Step 2:** Add reference sequences to the database's genome library.
```bash
k2 add-to-library --file <reference_sequence.fa> --db <custom_database>
```

**Step 3:** Build the custom database.
```bash
kraken2-build --build --db <custom_database>
```

You can also customize kmer lengths with `--kmer-len` and `--minimizer-len`. For example,
```bash
kraken2-build --build --db <custom_database> --kmer-len <kmer_length> --minimizer-len <minimizer_length>
```

More information about creating custom databases can be found [here](https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases).
