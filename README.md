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
An example config file is provided in `tests/example_config.yaml`. 

### Metadata file
An example metadata file is provided in `tests/example_metadata.csv`.
Metadata files should contain (at the minimum) a SampleID column, a forward reads filepath column (`forward_reads`), and a reverse reads filepath column (`reverse_reads`). These filepaths should be relative to the directory from which you are running `HoMi`.

### Conda environment building
If conda environments have already been built, and you'd like snakemake to not build them, pass the argument `--conda_prebuilt`.

### Running with example dataset
```
# Create the mock community (too big for github)
python tests/mock_community/create_mock_community.py

# Use the example config and example metadata provided
HoMi.py tests/example_config.yaml --cores 1
```

## Repository contents
- `snakefile` contains the bulk of the pipeline
- `HoMi/` contains `HoMi.py`, a wrapper controlling the behavior of the snakemake pipeline
- `conda_envs/` contains the conda environments for each rule in the snakemake pipeline
- `data/` contains relevant data, such as adapter sequences to be removed during trimming.

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