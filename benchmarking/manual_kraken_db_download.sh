#!/bin/bash

#SBATCH --partition=acpu
#SBATCH --account=amc-general
#SBATCH --output=/scratch/alpine/joconnor@xsede.org/HoMi/benchmarking/slurm_outputs/semi/%j.out
#SBATCH --job-name=SemiBM
#SBATCH --nodes=1 # use 1 node 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00 # Time limit days-hrs:min:sec
#SBATCH --qos=cpu-normal
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=john.2.oconnor@cuanschutz.edu



cd synthetic/data/kraken2_db_2
    

wget https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20240904.tar.gz
tar -xzvf k2_core_nt_20240904.tar.gz
rm k2_core_nt_20240904.tar.gz