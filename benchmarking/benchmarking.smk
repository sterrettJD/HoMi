import os
import pandas as pd

######## MADIS NOTES ########
## 1. might need to change filepaths to python/R scripts since subworkflows are now in sep directory? 
## 2. changed all partitions for benchmarking and homi (config files) to amilan 
## 3. added qos=normal as a default resource for homi_args (i think this will work but idk) - need to check strings 
## if this doesn't work, could add `default_slurm_extra: qos=normal` to ALL homi configs 

######## CONFIG ########

## IMPORTANT!! this is what determines which parts of benchmarking are run right now ##
## options for WANT_TO_RUN: "gen_refs", "synthetic", "semi", "pereira", "compare_all" ##
## define on the command line when you call snakemake like so: --config run=gen_refs  ##
WANT_TO_RUN = config["run"]

## executor=slurm (default) submits jobs to SLURM. executor=local runs on this computer. ##
## define on the command line when you call snakemake like so: --config executor=local ##
EXECUTOR = config.get("executor", "slurm")
IS_LOCAL = EXECUTOR == "local"

## this may need to be changed based on HPC being used!!
## define on the command line when you call snakemake like so: --config wanted_partition=amilan ##
hpc_partition = config.get("wanted_partition", "local")

## local runs skip HUMAnN by default, because it needs more memory than a desktop has ##
SKIP_HUMANN = str(config.get("skip_humann", IS_LOCAL)).lower() == "true"
## cores and memory (MB) for each nested HoMi run when executor=local ##
HOMI_CORES = int(config.get("homi_cores", 4))
HOMI_MEM_MB = int(config.get("homi_mem_mb", 24000))
## kraken memory (MB) when executor=local. The Standard-16 database needs about 16 GB. ##
KRAKEN_MEM_MB = int(config.get("kraken_mem_mb", 18000))
## HUMAnN and MetaPhlAn database setup memory (MB) when executor=local and skip_humann=False ##
HUMANN_MEM_MB = int(config.get("humann_mem_mb", 22000))

## nested HoMi targets for local runs. Explicit targets keep HoMi from running its rule all, ##
## which also asks for host read mapping and nonpareil outputs.                              ##
## The MultiQC reports read the trimmed reads, so they must be made before those are deleted. ##
LOCAL_HOMI_TARGETS = "aggregate_bracken reads_breakdown multiqc_pass1 multiqc_pass2"
LOCAL_HUMANN_TARGETS = "aggregate_humann_outs_nonhost"

SLURM_HOMI_ARGS = """--profile slurm --snakemake_extra "--jobs 40 --default-resources qos=cpu-normal --until host_filter" """


def build_local_homi_args(cores, mem_mb, kraken_mem_mb, humann_mem_mb, skip_humann):
    """
    Returns the HoMi.py arguments for a local run.
    Targets must come first in --snakemake_extra, because HoMi.py adds it at the end of
    the snakemake command, and --resources and --set-resources take a list of values.
    """
    set_resources = f"run_kraken:mem_mb={kraken_mem_mb}"
    if skip_humann:
        targets_and_config = f"{LOCAL_HOMI_TARGETS} --config skip_humann=True"
    else:
        targets_and_config = f"{LOCAL_HOMI_TARGETS} {LOCAL_HUMANN_TARGETS}"
        set_resources += (f" run_humann_nonhost:mem_mb={humann_mem_mb}"
                          f" setup_metaphlan:mem_mb={humann_mem_mb}")
    # Snakemake 7 cannot create conda environments with mamba 2, so use conda
    snakemake_extra = (f"{targets_and_config} "
                       f"--conda-frontend conda "
                       f"--resources mem_mb={mem_mb} "
                       f"--set-resources {set_resources}")
    return f'--cores {cores} --snakemake_extra "{snakemake_extra}"'


if IS_LOCAL:
    homi_args = build_local_homi_args(HOMI_CORES, HOMI_MEM_MB, KRAKEN_MEM_MB,
                                      HUMANN_MEM_MB, SKIP_HUMANN)
    # the outer snakemake counts this memory for each nested HoMi run
    homi_job_mem_mb = HOMI_MEM_MB
else:
    homi_args = SLURM_HOMI_ARGS
    homi_job_mem_mb = int(8*1000)

# CONDA ENVIRONMENT PATHS!
hostile_conda = "../conda_envs/hostile.yaml"
bbmap_conda = "../conda_envs/bbmap.yaml"
r_conda = "../conda_envs/r_env.yaml"
sra_conda = "../conda_envs/sra_tools.yaml"


# some global vars here
synthetic_work_dir = "synthetic"
synthetic_communities_dir = "synthetic_communities"
semi_work_dir = "semi"

# Synthetic communities metadata
metadata_file = os.path.join(synthetic_work_dir, "sample_data.csv")
metadata = pd.read_csv(metadata_file)
samples = metadata.drop(columns=["genome", "GCF_id"]).columns
reads = ["R1", "R2"]
organisms = metadata["genome"].to_list()
microbial_organisms = [x for x in organisms if (x != "human")]


# Semisynthetic communities metadata
semi_metadata_file = os.path.join(semi_work_dir, "sample_data.csv")
semi_metadata = pd.read_csv(semi_metadata_file)
semi_samples = semi_metadata.drop(columns=["genome", "SRR"]).columns
semi_organisms = semi_metadata["genome"].to_list()
semi_microbial_organisms = [x for x in semi_organisms if (x != "human")]
# SRR IDs are formatted as period separated lists within the SRR column
semi_srr_ids = [srr_id
                for taxon_srr_ids_list in semi_metadata["SRR"].apply(
                    lambda x: x.split(".")
                    ).values
                for srr_id in taxon_srr_ids_list]
semi_homi_args = homi_args


# mock community data
pereira_df = pd.read_csv("Pereira/Pereira_data.csv")
pereira_srr_ids = pereira_df["SRR"]
pereira_homi_args = homi_args

# Per nucleotide quality score for Polyester-simulated reads
## error rate changes depending on the phred score
polyester_phred_scores = [30, 40]
polyester_error_rate=0.001

# hostile reference data
custom_hostile_index_specs = ["rna", "hisat2"]

# host removal index option
indexes = ["dna", "rna", "hisat2"]


# metaphlan taxonomy comes from HUMAnN, so it is only plotted if HUMAnN runs
taxonomy_methods = ["kraken"] if SKIP_HUMANN else ["kraken", "metaphlan"]

# HoMi writes outputs to the PROJ path in each HoMi config, e.g. synthetic/homi_outputs/dna_transcriptomes_p30
synthetic_homi_proj_suffixes = {
    "synthetic": "",
    "synthetic_transcriptomes_p30": "_transcriptomes_p30",
    "synthetic_transcriptomes_p40": "_transcriptomes_p40"
}


def get_synthetic_homi_proj(wildcards):
    suffix = synthetic_homi_proj_suffixes[wildcards.proj]
    return os.path.join(synthetic_work_dir, "homi_outputs", f"{wildcards.index}{suffix}")


## this is not dependent on anything
gen_refs_out = [expand("t2t_hla_{hostile_index_spec}_index",
                        hostile_index_spec=indexes),
                "reference_genomes_downloaded"]

## this is dependent on gen refs 
synthetic_outs = [expand(os.path.join(synthetic_work_dir, synthetic_communities_dir, "{sample}_R1.fastq.gz"),
                      sample=samples),
                  expand(os.path.join(synthetic_work_dir, synthetic_communities_dir, "{sample}_R2.fastq.gz"),
                    sample=samples),
                  os.path.join(synthetic_work_dir, "synthetic_homi_metadata.csv"),
                  expand("{index}_synthetic_communities_benchmark.pdf",
                      index=indexes),
                  expand("{index}_synthetic_communities_benchmark_lm_results.txt",
                      index=indexes),
                  expand(os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}/combined/{sample}_R1.fastq.gz"),
                      sample=samples,
                      score=polyester_phred_scores),
                  expand(os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}/combined/{sample}_R2.fastq.gz"),
                      sample=samples,
                      score=polyester_phred_scores),
                  expand(os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}_homi_metadata.csv"),
                      score=polyester_phred_scores),
                  expand("{index}_HoMi_is_done_synthetic_transcriptomes_p{score}",
                      index=indexes,
                      score=polyester_phred_scores),
                  expand("{index}_index_synthetic_transcriptomes_p{score}_benchmark.pdf",
                        score=polyester_phred_scores,
                        index=indexes),
                  expand("{index}_index_synthetic_transcriptomes_p{score}_benchmark_lm_results.txt",
                        score=polyester_phred_scores,
                        index=indexes),
                  expand("taxonomy_compared/{index}/{proj}/{method}_genus.pdf",
                      index=indexes,
                      proj=["synthetic_transcriptomes_p30", 
                            "synthetic_transcriptomes_p40", 
                            "synthetic"],
                      method=taxonomy_methods),
                  expand("taxonomy_compared/{index}/{proj}/{method}_species.pdf",
                      index=indexes,
                      proj=["synthetic_transcriptomes_p30", 
                            "synthetic_transcriptomes_p40", 
                            "synthetic"],
                      method=taxonomy_methods)]

synthetic_homi_args = homi_args

## this is only dependent on gen refs 
semi_outs = [expand(os.path.join("semi", "samples", "{sample}_{read}.fastq.gz"),
               sample=semi_samples, read=reads),
             expand("{index}_semi_benchmark.pdf",
                    index=indexes),
             expand("{index}_semi_benchmark_lm_results.txt",
                    index=indexes)]


## this is only dependent on gen refs
pereira_outs = [expand(os.path.join("Pereira", "{srr_id}_R1.fastq.gz"),
                srr_id=pereira_srr_ids),
                expand(os.path.join("Pereira", "{srr_id}_R2.fastq.gz"),
                        srr_id=pereira_srr_ids),
                expand("{index}_Pereira_benchmark.pdf",
                    index=indexes),
                expand("{index}_Pereira_benchmark_lm_results.txt",
                    index=indexes),
                expand("{index}_Pereira_benchmark_from_paper.pdf",
                    index=indexes),
                expand("{index}_Pereira_benchmark_from_paper_lm_results.txt",
                    index=indexes)]
                
## this requires all other subworkflows to be run!
compare_all_outs = ["all_hostile_out_benchmark.pdf"]
# the combined taxonomy plots use metaphlan outputs from HUMAnN
if not SKIP_HUMANN:
    compare_all_outs.extend(["taxonomy_compared/combined_taxa_boxplot_genus.pdf",
                             "taxonomy_compared/combined_taxa_boxplot_species.pdf"])


## dictionary of outputs per sub-workflow 
sub_snake_outs = {
    "gen_refs": gen_refs_out,
    "synthetic": synthetic_outs,
    "pereira": pereira_outs,
    "semi": semi_outs,
    "compare_all": compare_all_outs
}

## dictionary of paths to sub-workflows 
sub_smk_paths = {
    "gen_refs": "rules/01_gen_refs.smk",
    "synthetic": "rules/02a_synthetic.smk",
    "pereira": "rules/02b_pereira.smk",
    "semi": "rules/02c_semi.smk",
    "compare_all": "rules/03_compare_all.smk"
}

## put together hierarchy based on which sub-workflows are dependent on each other 
sub_snake_hierarchy = {
    "gen_refs": ["gen_refs"],
    "synthetic": ["gen_refs", "synthetic"],
    "pereira": ["gen_refs", "pereira"],
    "semi": ["gen_refs", "semi"],
    "compare_all": ["gen_refs", "synthetic", "semi", "pereira", "compare_all"] 
}

## pulling which subworkflows to run 
snakes_to_run = sub_snake_hierarchy[WANT_TO_RUN]

## empty rule all input list
rule_all_input_list = []

## add rule_all outputs to list and include the correct subworkflows based on the 
## sub-workflows specified (or lack of)
for snake in snakes_to_run:
    rule_all_input_list.extend(sub_snake_outs[snake])
    include: sub_smk_paths[snake]

## finally call rule all 
rule all:
    input:
        rule_all_input_list
