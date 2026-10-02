#!/bin/bash

#SBATCH --partition=acpu
#SBATCH --account=amc-general
#SBATCH --output=/scratch/alpine/joconnor@xsede.org/HoMi/benchmarking/slurm_outputs/database_db/%j.out
#SBATCH --job-name=Metaphlan
#SBATCH --nodes=1 # use 1 node 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00 # Time limit days-hrs:min:sec
#SBATCH --qos=cpu-normal
#SBATCH --mem=150G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=john.2.oconnor@cuanschutz.edu


module load miniforge

conda activate /gpfs/alpine1/scratch/.xsede.org/joconnor/HoMi/benchmarking/.snakemake/conda/20743d933e9ce901648d0deb9f9391dd_

mkdir -p synthetic/data/metaphlan_db

metaphlan --install --nproc 8 --bowtie2db synthetic/data/metaphlan_db