#!/bin/bash

## to capture output and place in log file run: 
## bash run_benchmarking.sh > compareAll_dryRun.log 2>&1


#First unlocl :

snakemake \
    -s benchmarking.smk \
    --config run=synthetic wanted_partition=acpu \
    --cores 2 \
    --use-conda \
    --conda-prefix /gpfs/alpine1/scratch/.xsede.org/joconnor/HoMi/.snakemake/conda/ \
    --default-resources qos=cpu-normal \
    --profile /home/.xsede.org/joconnor/.config/snakemake/slurm \
    --unlock
    
## config RUN options: "gen_refs", "synthetic", "semi", "pereira", "compare_all"
## we mainly just want to run semi since that's what's failing on john's end - need to make sure
## that splice-aware aligner (hisat2) output files aren't blank after homi (i think)
snakemake \
    -s benchmarking.smk \
    --config run=synthetic wanted_partition=acpu \
    --cores 2 \
    --use-conda \
    --conda-prefix /gpfs/alpine1/scratch/.xsede.org/joconnor/HoMi/.snakemake/conda/ \
    --default-resources qos=cpu-normal \
    --profile /home/.xsede.org/joconnor/.config/snakemake/slurm

## attempt to set default qos for benchmarking to "normal"
## will still have to fix default resources for homi to include qos as well 
    ##--default-resources qos=normal \