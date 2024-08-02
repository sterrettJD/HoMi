# Benchmarking HoMi

This directory includes helper scripts to benchmark HoMi against simulated fully synthetic and semi-synthetic communities.

## Environment
`benchmarking_env.yaml` contains the conda environment used for benchmarking. To create this environment, run the following command from this `benchmarking/` directory:
```
mamba env create -n benchmarking-homi --file benchmarking_env.yaml
pip install -e ../
```

**NOTE:** Users need to install the SRA toolkit outside of conda, as SRA-tools doesn't support a conda distribution. See the [SRA toolkit official instructions](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit).

NOTE: Because `ncbi-datasets-cli` doesn't have a MacOS ARM-compatible conda-forge distribution, MacOS ARM users need to create the environment and configure the subdir with the following commands:
```
conda create -n benchmarking-homi
conda activate benchmarking-homi 
conda config --env --set subdir osx-64
conda env update --file benchmarking_env.yaml
```

Alternatively, this can be run on ARM using Docker via a command such as
```
docker run --platform linux/amd64 -v "$(pwd)":/workdir:rw -w /workdir snakemake/snakemake:v7.32.3 snakemake -s benchmarking.smk --cores 4 --use-conda
```

## Running it all

`benchmarking.smk` is the snakemake file that will run all of the benchmarking. Run it at your leisure via snakemake (e.g., `snakemake -s benchmarking.smk -c8` for 8 cores or ` snakemake -s benchmarking.smk --profile slurm` if you have a slurm profile set up and want to run it with that)

## Synthetic communities

## Semi-synthetic communities

## Mock communities

