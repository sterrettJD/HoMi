# Benchmarking HoMi

This directory includes helper scripts to benchmark HoMi against simulated fully synthetic and semi-synthetic communities.

## Environment
`benchmarking_env.yaml` contains the conda environment used for benchmarking. To create this environment, run the following command from this `benchmarking/` directory:
```
mamba env create -n benchmarking-homi --file benchmarking_env.yaml
pip install homi-pipeline
```

**NOTE:** Users need to install the SRA toolkit outside of conda, as SRA-tools doesn't support a conda distribution. See the [SRA toolkit official instructions](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit).

This can be run on ARM using Docker via a command such as
```
docker run --platform linux/amd64 -v "$(pwd)":/workdir:rw -w /workdir snakemake/snakemake:v7.32.3 snakemake -s benchmarking.smk --cores 4 --use-conda
```

## Running it all

`benchmarking.smk` is the snakemake file that will run all of the benchmarking. Run it at your leisure via snakemake (e.g., `snakemake -s benchmarking.smk -c8` for 8 cores or ` snakemake -s benchmarking.smk --profile slurm` if you have a slurm profile set up and want to run it with that)

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
