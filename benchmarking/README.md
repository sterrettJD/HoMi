# Benchmarking HoMi

This directory includes helper scripts to benchmark HoMi against simulated fully synthetic and semi-synthetic communities.

## Environment
`benchmarking_env.yaml` contains the conda environment used for benchmarking. To create this environment, run the following command from this `benchmarking/` directory:
```
mamba env create -n benchmarking-homi --file benchmarking_env.yaml
pip install -e ..
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

## Run the Pipeline

The workflows are `gen_refs`, `semi`, `pereira`, `synthetic` and `compare_all`.
Run the commands from this `benchmarking/` directory.

### Set Up

1. Create the benchmarking environment:
   ```bash
   mamba env create -n benchmarking-homi --file benchmarking_env.yaml
   ```
2. Activate the environment:
   ```bash
   conda activate benchmarking-homi
   ```
3. Install HoMi from this repository:
   ```bash
   pip install -e ..
   ```
4. Create the `hostile` environment. HoMi uses this environment name for host read removal:
   ```bash
   mamba env create -n hostile --file ../src/homi_pipeline/conda_envs/hostile.yaml
   ```

### Run on a Local Computer

The local mode runs each step on this computer. HoMi makes the MultiQC reports, the Kraken2 results and the reads breakdown report. By default, it does not run MetaPhlAn or HUMAnN. It uses the Kraken2 Standard-16 database, which needs approximately 16 GB of memory. The run script sets a hard memory limit. If the pipeline uses more memory than the limit, the system stops the pipeline.

1. Turn on linger. Linger keeps the memory limit and the pipeline active after you disconnect from SSH. This is primarily useful if you SSH into this computer.
   ```bash
   loginctl enable-linger $USER
   ```
2. Start a tmux session:
   ```bash
   tmux new -s bench
   ```
3. Download the Kraken2 database:
   ```bash
   bash download_kraken_standard16.sh
   ```
4. Activate the environment:
   ```bash
   conda activate benchmarking-homi
   ```
5. Run the workflows in this sequence:
   ```bash
   bash run_benchmarking_local.sh gen_refs
   bash run_benchmarking_local.sh semi
   bash run_benchmarking_local.sh pereira
   bash run_benchmarking_local.sh synthetic
   bash run_benchmarking_local.sh compare_all
   ```
   To give more options to Snakemake, write them after the workflow name. Example: `bash run_benchmarking_local.sh semi --rerun-triggers mtime`
6. To disconnect from the tmux session, push `Ctrl+b`, then `d`. You can then log off.

To examine a run after you log in again:
* To see the output, connect to the tmux session:
  ```bash
  tmux attach -t bench
  ```
* To read the log, type `tail -f local_<workflow>.log`.
* To see the memory use, type:
  ```bash
  systemctl --user status 'homi-bench-*'
  ```

#### Memory and CPU Limits

The default limits are for a computer with 32 GB of memory and 4 CPU cores. Memory values are in MB, except `CGROUP_MAX`.

| Variable | Default | Description |
| --- | --- | --- |
| `CORES` | 4 | CPU cores for Snakemake and for each HoMi run |
| `HOMI_MEM_MB` | 24000 | Memory for each HoMi run |
| `KRAKEN_MEM_MB` | 18000 | Memory for each Kraken2 job in HoMi |
| `SKIP_HUMANN` | `true` | Set to `false` to also run MetaPhlAn and HUMAnN |
| `HUMANN_MEM_MB` | 22000 | Memory for each HUMAnN job in HoMi. Must not be more than `HOMI_MEM_MB`. |
| `SMK_MEM_MB` | 26000 | Memory for all jobs that run at the same time. Must be at least `HOMI_MEM_MB`. |
| `CGROUP_MAX` | `28G` | Hard memory limit. Keep it a few GB below the total memory. |

Example for a computer with 64 GB of memory and 8 CPU cores:
```bash
CORES=8 HOMI_MEM_MB=40000 SMK_MEM_MB=48000 CGROUP_MAX=56G bash run_benchmarking_local.sh semi
```

#### Run MetaPhlAn and HUMAnN

Set `SKIP_HUMANN=false`. You can do this after a run without HUMAnN. HoMi then runs only the MetaPhlAn and HUMAnN steps and the reads breakdown report. It does not run the read trimming, host read removal or Kraken2 steps again.
```bash
SKIP_HUMANN=false bash run_benchmarking_local.sh semi
```
The first run downloads the MetaPhlAn, ChocoPhlAn, UniRef90 and utility mapping databases to `synthetic/data/`. These databases need approximately 50 GB of disk space. HUMAnN also keeps large temporary files for each sample. To delete them after a run, use `HoMi_cleanup.py <HoMi config> <HoMi metadata> --humann`.

### Run on an HPC with SLURM

1. Install a Snakemake SLURM profile named `slurm`.
2. Download the Kraken2 database with `manual_kraken_db_download.sh`.
3. Run each workflow. Replace `<workflow>`, `<partition>` and `<qos>`:
   ```bash
   snakemake -s benchmarking.smk \
       --config run=<workflow> wanted_partition=<partition> \
       --profile slurm \
       --default-resources qos=<qos> \
       --use-conda
   ```

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
