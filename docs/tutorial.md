# HoMi tutorial

## Tutorial data

Data for this tutorial are found in the `docs/tutorial_data/` directory of the HoMi repo ([here](https://github.com/sterrettJD/HoMi/tree/main/docs/tutorial_data)). 

In this directory, we have a `raw_data` folder with the fastq files for 2 samples (paired end reads). There is also an example config file and an example metadata file. Below are steps to get these files.

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

### Set up your profile

You can use the provided profile setup script to setup a cookiecutter Slurm profile for Snakemake, using the following command. This profile interfaces with Snakemake within HoMi to schedule the jobs for each step of the pipeline with a job scheduler system on a compute cluster.

```
profile_setup.py --cluster_type slurm-smt --output_dir ./tutorial_slurm_profile
```

Now, in the `tutorial_slurm_profile/` directory, you have a `slurm` directory that is your profile for integrating snakemake with Slurm.

### Set up the metadata

At a baseline, your metadata for HoMi should be a comma separated value (.csv) file with three columns:

1. Sample - This should be named exactly `Sample` and denotes the sample ID
2. forward_reads - The exact name of this should be specified in the config file, but we recommend `forward_reads`. This column contains the filepaths to the forward reads files, which can be absolute paths or paths relative to the directory from which you are running HoMi.
3. reverse_reads - The exact name of this should be specified in the config file, but we recommend `reverse_reads`. This column contains the filepaths to the reverse reads files, which can be absolute paths or paths relative to the directory from which you are running HoMi.

Here is an example:

| Sample | forward_reads | reverse_reads |
| --- | --- | --- |
| sample_1 | raw_data/tutorial_sample_1.R1.fq.gz | raw_data/tutorial_sample_1.R2.fq.gz |
| sample_2 | raw_data/tutorial_sample_2.R1.fq.gz | raw_data/tutorial_sample_2.R2.fq.gz |


#### Host mapping

An optional column can be provided in the metadata, labeled `map_host`, which can be used to denote samples that shouldn't be mapped to a host genome/transcriptome. This is useful if you have purely metagenomic or metatranscriptomic samples that you don't care about host gene/transcript counts in.

Passing a value of `False` to `map_host` will result in those samples being run through the whole HoMi pipeline (including host decontamination), but there will be no quantification of host genes/transcripts.

For example, if we didn't care about host gene counts in our previous metadata file, we could use the following metadata file:

| Sample | forward_reads | reverse_reads | map_host |
| --- | --- | --- | --- |
| sample_1 | raw_data/tutorial_sample_1.R1.fq.gz | raw_data/tutorial_sample_1.R2.fq.gz | False |
| sample_2 | raw_data/tutorial_sample_2.R1.fq.gz | raw_data/tutorial_sample_2.R2.fq.gz | False |

#### Tutorial metadata
The tutorial metadata file should be found in your `HoMi_tutorial` directory that you've copied the tutorial data into, as `HoMi_tutorial/tutorial_metadata.csv`.

### Set up the config file

The config file for HoMi is a YAML-formatted file. A variety of parameters are expected in the config file. These are checked by `check_config.py` when you run `HoMi.py`, but you can also run `check_config.py` on its own to validate your config file.

#### Tutorial config file
The tutorial config file should be found in your `HoMi_tutorial` that you've copied the tutorial data into, as `HoMi_tutorial/tutorial_HoMi_config.yaml`. It can be helpful to use this file as a backbone for your own data, and modify the parameters as you see fit.

Parameters are explained below:

#### Project name

```
PROJ: tutorial
```

This variable (`PROJ`) is the first part of the file output path, in which all results will be stored. For example, reads that have been QCed (no trimming) and host-filtered will be put into the directory `tutorial.f0.0.r0.0.nonhost`.


#### Metadata file
```
METADATA: tutorial_metadata.csv
fwd_reads_path: forward_reads
rev_reads_path: reverse_reads
```

`METADATA` is the path to your metadata file, relative to the directory you're running HoMi from. It could also be an absolute path.

`<fwd/rev>_reads_path` are the column names in the metadata for the columns that contain paths to forward and reverse reads.


#### Quality control parameters

Parameters are explained below in the code chunk.

```
################################################
### TRIMMOMATIC ADAPTER AND QUALITY TRIMMING ###
# minimum read length after trimming
min_readlen: 50
# the minimum quality for the start of a read. If it's below this quality, trim that base pair
readstart_qual_min: 20
# the minimum quality for the end of a read. If it's below this quality, trim that base pair
readend_qual_min: 20

######################
### SEQTK TRIMMING ###
# This is a second pass hard trimming, after trimmomatic has done quality-based and adapter trimming
# Read trimming this many base pairs from the start of the read
trim_fwd: 0
trim_rev: 0

# Read truncating this many base pairs from the end of the read
trunc_fwd: 0
trunc_rev: 0
```

#### Databases

Parameters are explained below in the code chunk. If databases do not exist in the specified location, they will be downloaded.

Notably, the `host_tax_id` parameter should be a NCBI taxon ID. It is used when making a stacked barplot from the Kraken/Bracken microbial taxonomy, to filter out host reads from the visualization. It is an optional parameter.

```
#################################
### Hostile host read removal ###
# Database for host read removal
# Current supported options include:
# - human-t2t-hla
# - human-t2t-hla-argos985
# - filepath to an already downloaded and bowtie2-indexed database 
#   (with no .bt1, .bt2 etc file extensions in this argument)
hostile_db: human-t2t-hla-argos985

# Where to download a database for hostile (if applicable)
# This example will result in the following database location
# data/hostile-t2t-hla-argos985.bt1, data/hostile-t2t-hla-argos985.bt2, etc.
loc_for_hostile_db_download: data


#################################
### HUMAnN BioBakery pipeline ###

# Path to MetaPhlan bowtie database (step 1 of HUMAnN)
metaphlan_bowtie_db: data/metaphlan_db/

# Paths to Chocophlan and UniRef databases for HUMAnN
# If these aren't already downloaded, HoMi will download them
chocophlan_db: data/humann_dbs/chocophlan
uniref_db: data/humann_dbs/uniref
utility_mapping_db: data/humann_dbs/utility_mapping 


########################
### Kraken + Bracken ###
# database location
kraken_db: data/kraken2_db
# host taxon ID to not be plotted in microbial taxa barplot
host_tax_id: 9606


#################################
### Host read mapping ###
# provide the path for the host reference genome
# if these files do not exist yet, HoMi will by default download the GRCh38 human reference genome
host_ref_fna: data/GRCh38/GRCh38_full_analysis_set.fna 
host_ref_gtf: data/GRCh38/GRCh38_full_analysis_set.refseq.gtf
```

#### Host map method

A host mapping method should be specified, even if `map_host` is `False` for all samples (just pick one arbitrarily in that case). Options include BBMap and HISAT2. 

```
host_map_method: HISAT2
```

#### Resources
Rules will use the default resources outlined for them in `src/homi_pipeline/snakefile`, unless you provide other resource specifications here. These should be formatted as `<rule_name>_<resource>`. Rule descriptions can be accessed [here](https://homi-pipeline.readthedocs.io/en/latest/content/all_rules.html).


```
rule_name_partition: new_partition
rule_name_mem_mb: 10000 # (10 GB)
rule_name_runtime: 600 # (10 hours)
rule_name_threads: 8 
```

Additionally, HoMi passes default partitions to clusters, as `short` and `long`. These default partitions can be updated. For example, FIJI HPC users at the University of Colorado's BioFrontiers Institute can stick with the default partition names, but Alpine HPC users in the University of Colorado/Colorado State University system need to update this as follows:

```
default_short_partition_name: amilan
default_long_partition_name: amilan
```

For long runtime steps, or for setting default email notifications on a Slurm-managed cluster, users may want to update Slurm parameters. 

This can be done using `default_slurm_extra`, which will serve as a default for all steps, or by using `<rule_name>_slurm_extra`, which provides rule-specific slurm extra parameters. Rule-specific parameters are prioritized over default params.


These parameters should be a space-separated list, formatted as `<slurm_long_name>=<parameter> <slurm_long_name_2>=<parameter_2>`. Here is an example, for if I wanted to get an email about each job, and I wanted to run the `map_host` step with a long quality of service characteristic (relevant for Alpine HPC users).

```
default_slurm_extra: email=sample_email@colorado.edu
map_host_slurm_extra: qos=long email=sample_email@colorado.edu
```

#### Rule-specific extra arguments

Sometimes, you'll need to provide extra arguments to a specific rule. This can be done with a `<rule_name>_extra` parameter.

For example, if you are RAM-limited and don't want to load the full Kraken database into memory, you can pass the argument `--memory-mapping` to kraken2, you should include the following line in your config file:

```
run_kraken_extra: " --memory-mapping "
```


### Running HoMi
You can run HoMi with the following command. This should do everything for you. 
```
HoMi.py tutorial_HoMi_config.yaml --profile tutorial_slurm_profile/slurm
```