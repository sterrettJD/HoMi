# HoMi tutorial

## Tutorial data

Data for this tutorial are found in the `docs/tutorial_data/` directory of the HoMi repo (https://github.com/sterrettJD/HoMi/tree/main/docs/tutorial_data). 

In this directory, we have a `raw_data` folder with the fastq files for 2 samples (paired end reads).

### Get data from the cloned HoMi repo (preferred)

This assumes you have GitHub ssh keys set up. If not, see this [instructions guide](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/checking-for-existing-ssh-keys).
```
# Clone the repo
git clone git@github.com:sterrettJD/HoMi.git

# Copy the tutorial data over
cp -r HoMi/docs/tutorial_data/ ./HoMi_tutorial

# cleanup - delete the HoMi repo
rm -r HoMi
```

### ALTERNATIVE: get data using wget

If you don't have ssh keys set up for GitHub, you can use the following commands:
```
# Download the repo and unzip it
wget https://github.com/sterrettJD/HoMi/archive/refs/heads/main.zip
unzip main.zip

# Copy the tutorial data over
cp -r HoMi-main/docs/tutorial_data/ ./HoMi_tutorial

# cleanup - delete the HoMi repo
rm main.zip
rm -r HoMi-main
```

## Create your environment

```
# create and activate the conda environment
mamba create -n HoMi_tutorial python=3.11
conda activate HoMi_tutorial
# install HoMi
pip install homi-pipeline
```

## Move into the tutorial directory

```
cd HoMi_tutorial
```

## Run HoMi on a Slurm-managed cluster

### Setup your profile

You can use the provided profile setup script to setup a cookiecutter Slurm profile for snakemake, using the following command.
```
profile_setup.py --cluster_type slurm-smt --output_dir ./tutorial_slurm_profile
```

Now, in the `tutorial_slurm_profile/` directory, you have a `slurm` directory that is your profile for integrating snakemake with Slurm.

### Running HoMi
You can run HoMi with the following command. This should do everything for you.
```
HoMi.py tutorial_HoMi_config.yaml --profile tutorial_slurm_profile/slurm
```