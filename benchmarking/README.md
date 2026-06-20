# Benchmarking HoMi

This directory includes helper scripts to benchmark HoMi against simulated fully synthetic and semi-synthetic communities.

## Environment
`benchmarking_env.yaml` contains the conda environment used for benchmarking. To create this environment, run the following command from this `benchmarking/` directory:
```
mamba env create -n benchmarking-homi --file benchmarking_env.yaml
pip install homi-pipeline
```

**NOTE:** Users need to install the SRA toolkit outside of conda, as SRA-tools doesn't support a conda distribution. See the [SRA toolkit official instructions](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit).

This pipeline can alternatively be run on ARM using Docker via a command such as
```
docker run --platform linux/amd64 -v "$(pwd)":/workdir:rw -w /workdir snakemake/snakemake:v7.32.3 snakemake -s benchmarking.smk --cores 4 --use-conda
```

## Pipeline Structure

The benchmarking pipeline is modularized into independent workflows in the `rules/` directory:
- `01_gen_refs.smk` - Generate reference genomes and indices
- `02a_synthetic.smk` - Benchmark synthetic communities
- `02b_pereira.smk` - Benchmark mock communities (Pereira et al.)
- `02c_semi.smk` - Benchmark semi-synthetic data
- `03_compare_all.smk` - Compare results across all datasets

This structure allows you to run specific workflows independently or together, reducing computational overhead and improving debuggability.

## Running Workflows

You can now run specific workflows rather than the full pipeline:

```bash
# Run just reference generation
snakemake -s benchmarking.smk --config run=gen_refs wanted_partition=amilan --cores 2 --use-conda

# Run semi-synthetic workflow (useful for focused testing)
snakemake -s benchmarking.smk --config run=semi wanted_partition=amilan --cores 2 --use-conda

# Run all benchmarks and compare
snakemake -s benchmarking.smk --config run=compare_all wanted_partition=amilan --profile slurm --use-conda
```

**Available workflows:** `gen_refs`, `synthetic`, `semi`, `pereira`, `compare_all`

For SLURM environments such as the Alpine HPC (Colorado), add `--default-resources qos=normal` to the command.

### Legacy Scripts

`benchmark_host_read_removal_method.smk` does some alternative benchmarking of hostile using HISAT2 via updates to hostile on a forked branch of the repo. This will hopefully become redundant if HISAT2 is incorporated into hostile, but for now it relies on running this pipeline after the main `benchmarking.smk`. It can be run using `snakemake -s benchmark_host_read_removal_method.smk --profile slurm`.


## Synthetic communities
Synthetic communities were generated in two ways.

### Synthetic transcriptome simulation with Polyester

Polyester was used to simulate transcriptomes from the human transcriptome, using the GRCh38 reference. Details on polyester can be found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4635655/).

### Simulated communities from genomes including the human pangenome

A custom script (`synthetic/create_mock_community.py`) was used to create communities with reads from the [human pangenome](https://humanpangenome.org/). Errors in the reads were simulated based on a mean PHRED score of 35 with a standard deviation of 3 at each base position. This could be improved, but the main point of this simulation was to make sure HoMi recovered proper portions of host reads when using "noisier" human reads than what are provided by the GRCh reference genome/transcriptome.

## Semisynthetic transcriptomes

Samples were simulated containing real transcriptomic data combined in known portions. Specified numbers of reads were subsampled from publicly available FASTQ files from bacterial isolate studies and human colon chip samples. The SRRs used for this project can be found in `semi/sample_data.csv`. These data are downloaded and subsampled as part of `benchmarking.smk`.

## Mock communities

Mock communities were pulled from the [Pereira-Marques et al. low biomass paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10913719/). These mock communities include 97%, 90%, 70%, 10%, and 0% host cells combined with a mock microbial community, then sequenced. 

`Pereira/Pereira_data.csv` contains the sample host percents and SRR accession IDs. Running `benchmarking.smk` will pull these SRR runs using prefetch+fasterq-dump to the `Pereira/` directory.

## Directory Structure

```
benchmarking/
├── benchmarking.smk          # Main modular Snakefile (workflow manager)
├── benchmark_host_read_removal_method.smk  # Legacy HISAT2 benchmarking
├── README.md
├── conda_envs/               # Conda environment YAML files
│   ├── hostile.yaml
│   ├── hostile_dev.yaml
│   ├── bbmap.yaml
│   ├── r_env.yaml
│   └── sra_tools.yaml
├── rules/                    # Modularized workflow rules
│   ├── 01_gen_refs.smk
│   ├── 02a_synthetic.smk
│   ├── 02b_pereira.smk
│   ├── 02c_semi.smk
│   └── 03_compare_all.smk
├── synthetic/                # Synthetic community data and configs
│   ├── sample_data.csv
│   ├── create_mock_community.py
│   └── *_HoMi_config.yaml    # HoMi configs for different indices
├── semi/                     # Semi-synthetic community data and configs
│   ├── sample_data.csv
│   └── *_HoMi_config.yaml
└── Pereira/                  # Mock community data (Pereira et al.)
    ├── Pereira_data.csv
    └── *_HoMi_config.yaml
```
