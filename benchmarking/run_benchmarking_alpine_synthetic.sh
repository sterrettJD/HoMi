#!/bin/bash

#SBATCH --partition=acpu
#SBATCH --account=amc-general
#SBATCH --output=/scratch/alpine/joconnor@xsede.org/HoMi/benchmarking/slurm_outputs/semi/%j.out
#SBATCH --job-name=SynBM
#SBATCH --nodes=1 # use 1 node 
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00 # Time limit days-hrs:min:sec
#SBATCH --qos=cpu-normal
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=john.2.oconnor@cuanschutz.edu

module purge 
module load miniforge

conda activate homi_benchmarking


## to capture output and place in log file run: 
bash run_benchmarking_synthetic.sh > compareAll_ActualRun_synthetic.log 2>&1
