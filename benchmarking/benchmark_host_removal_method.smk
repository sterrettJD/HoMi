import os
import pandas as pd


######## CONFIG ########
# some global vars here
synthetic_communities_dir = "synthetic_communities"
synthetic_transcriptomes_dir = "synthetic_transcriptomes"
semi_work_dir = "semi"

sim_methods = ["synthetic_communities", "synthetic_transcriptomes", "semi"]


# hostile reference data
index_methods = ["bowtie2", "hisat2"]

rule all:
    input:
        expand("{method}_index",
                method=index_methods)


rule create_alt_hostile_index:
    output:
        ref_dir=directory("{method}_index")
    threads: 4
    conda: "conda_envs/hostile.yaml"
    resources:
        partition="short",
        mem_mb=int(12*1000), # MB
        runtime=int(6*60) # min
    params:
        script="create_decontam_ref_human.py"
    shell:
        """
        python {params.script} -m {wildcards.method} -o {output.ref_dir}
        """
