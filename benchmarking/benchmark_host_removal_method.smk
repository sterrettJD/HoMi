import os
import pandas as pd


######## CONFIG ########
# some global vars here
synthetic_dir = "synthetic"
synthetic_communities_dir = "synthetic"
synthetic_transcriptomes_dir = "synthetic_transcriptomes"
semi_work_dir = "semi"

sim_methods = ["synthetic", "synthetic_transcriptomes", "semi"]

semi_metadata_file = os.path.join(semi_work_dir, "semi_homi_metadata.csv")
semi_metadata = pd.read_csv(semi_metadata_file)


metadata_files = [os.path.join(synthetic_dir, 
                                f"{synthetic_communities_dir}_homi_metadata.csv"),
                  os.path.join(synthetic_dir, 
                                f"{synthetic_transcriptomes_dir}_homi_metadata.csv"),
                  os.path.join(semi_work_dir, 
                                "semi_homi_metadata.csv")
                ]

metadata_dfs = [pd.read_csv(file) for file in metadata_files]
sample_paths = [path.replace("_R1.fastq.gz", "").replace("_R2.fastq.gz", "") 
                for df in metadata_dfs 
                for path in df["forward_reads"].to_list()]
project = ["synthetic_communities"] * metadata_dfs[0].shape[0] + \
          ["synthetic_transcriptomes"] * metadata_dfs[1].shape[0] + \
          ["semi"] * metadata_dfs[2].shape[0]
forward_reads = [path for df in metadata_dfs for path in df["forward_reads"].to_list()]
reverse_reads = [path for df in metadata_dfs for path in df["reverse_reads"].to_list() ]


# hostile reference data
index_methods = ["bowtie2", "hisat2"]

rule all:
    input:
        expand("{method}_index",
                method=index_methods),
        expand(os.path.join("{method}_data", 
                "{sample}_R1.fastq.gz"),
                method=index_methods, sample=sample_paths),
        expand(os.path.join("{method}_data", 
                "{sample}_R2.fastq.gz"),
                method=index_methods, sample=sample_paths),
        expand(os.path.join("{method}_data",
                "{sample}.report"),
                method=index_methods, sample=sample_paths),
        "host_filtering_benchmark_aggregated.csv"


rule create_alt_hostile_indexes:
    output:
        ref_dir=directory("{method}_index")
    threads: 4
    conda: "conda_envs/hostile_dev.yaml"
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


rule run_hostile:
    input:
        index=expand("{method}_index",
                method=index_methods),
        fwd="{sample}_R1.fastq.gz",
        rev="{sample}_R2.fastq.gz"
    output:
        fwd=os.path.join("{method}_data", "{sample}_R1.fastq.gz"),
        rev=os.path.join("{method}_data", "{sample}_R2.fastq.gz"),
        report=os.path.join("{method}_data", "{sample}.report")
    threads: 12
    conda: "conda_envs/hostile_dev.yaml"
    resources:
        partition="short",
        mem_mb=int(12*1000), # MB
        runtime=int(8*60) # min
    shell:
        """
        mkdir -p {wildcards.method}_data

        hostile clean \
        --fastq1 {input.fwd} --fastq2 {input.rev} \
        -o {wildcards.method}_data/{wildcards.sample}/.. \
        --threads {threads} \
        --index {wildcards.method}_index/index \
        --debug \
        --aligner {wildcards.method} > {output.report}

        # cleanup filepath
        mv {wildcards.method}_data/{wildcards.sample}_R1.clean_1.fastq.gz {output.fwd}
        mv {wildcards.method}_data/{wildcards.sample}_R2.clean_2.fastq.gz {output.rev}
        """

rule aggregate_hostile_reports:
    input:
        reports=expand(os.path.join("{method}_data",
                "{sample}.report"),
                method=index_methods, sample=sample_paths)
    output:
        report="host_filtering_benchmark_aggregated.csv"
    threads: 1
    conda: "conda_envs/hostile_dev.yaml"
    resources:
        partition="short",
        mem_mb=int(1*1000), # MB
        runtime=int(1*60) # min
    params:
    shell:
        import json
        import pandas as pd

        vars_of_interest = ["fastq1_in_name", "aligner", 
                            "reads_removed_proportion"]
        data = dict()
        for report in input.reports:
            with open(report) as f:
                report_data = json.load(f)[0]
            
            data[report] = [report_data[var] for var in vars_of_interest]
        
        df = pd.DataFrame.from_dict(data, 
                                    orient="index", 
                                    columns=vars_of_interest)
        df.to_csv(output.report)

