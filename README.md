# HoMi
Host-microbiome dual transcriptome pipeline

## Installation
```
git clone git@github.com:sterrettJD/HoMi.git
cd HoMi
pip install -e .
```

## Usage
```
HoMi.py <config_file> --cores <n_cores> --profile <profile_name>
```

### Config file
An example config file is provided in `tests/example_config.yaml`. This config file should contain paths to relevant files, like the metadata and databases. It can also be used to alter rule-specific resource requirements (details in example config).

### Metadata file
An example metadata file is provided in `tests/example_metadata.csv`.
Metadata files should contain (at the minimum) a Sample column (named `Sample`), a forward reads filepath column (column name specified in the config file under `fwd_reads_path`), and a reverse reads filepath column (column name specified in the config file under `rev_reads_path`). These filepaths should be relative to the directory from which you are running `HoMi`.

### Conda environment building
If conda environments have already been built, and you'd like snakemake to not build them, pass the argument `--conda_prebuilt`.

### Running with example dataset
```
# Create the mock community (too big for github)
python tests/mock_community/create_mock_community.py

# Use the example config and example metadata provided
HoMi.py tests/example_config.yaml --cores 1
```

### Unlocking a snakemake directory
Sometimes, when snakemake unexpectedly exits (e.g., due to a server connection timeout), the directory may be locked. Pass the argument `--unlock` to unlock the directory before running `HoMi.py`.

## Pipeline steps
### Preprocessing
1. Create a symbolic link to the sequencing files.
2. Trim reads using Trimmomatic
    - Trims adapter sequences from reads
    - Trims reads with a starting PHRED quality below `readstart_qual_min`
    - Trims reads with an ending PHRED quality below `readend_qual_min`
    - Trims reads with wherever there's a 4-base sliding window average PHRED quality score below 20
    - Removes reads with a length below `min_readlen`
3. Generates quality report with FastQC + MultiQC
4. Trim read ends using Seqtk
    - Consider this a second pass, in case Trimmomatic didn't catch something
5. Generates a second pass quality report with FastQC + MultiQC after the second trimming step

### Read mapping
6. Remove host reads using Hostile
    - User should pass a database in the config file. Currently supported options are `human-t2t-hla` and `human-t2t-hla-argos985`.
    - OR users can pass a filepath to a bowtie2 index without the `.bt` extensions (e.g., `index/example_index`, where files exist named `index/example_index.bt1`, `index/example_index.bt2`, etc.).
7. Run HUMAnN pipeline on nonhost reads to profile microbial reads
    - Makes microshades taxa barplot from MetaPhlan and HUMAnN outputs
8. Align all reads against the host genome
    - HoMi will by default download the GRCh38 human reference genome, but you can provide an alternative genome (fna + gtf) if it's already downloaded
    - BBmap is used to map the reads, and featureCounts is used to generate a read count table



## Main repository contents
- `snakefile` contains the bulk of the pipeline
- `HoMi/` contains `HoMi.py`, a wrapper controlling the behavior of the snakemake pipeline
- `conda_envs/` contains the conda environments for each rule in the snakemake pipeline
- `data/` contains relevant data, such as adapter sequences to be removed during trimming.
