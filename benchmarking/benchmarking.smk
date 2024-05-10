import os
from semisynthetic import simulate_semisynthetic as sss


# some global vars here for now BUT NEED TO BE MIGRATED TO CONFIG
synthetic_work_dir = "synthetic"
synthetic_communities_dir = "synthetic_communities"
homi_args = "--profile slurm"

#colon_sample_htx="SRP127360"
#colon_sample_htx="SRR6410603"


rule all:
    input:
        # From simulate_synthetic_communities
        os.path.join(synthetic_work_dir, synthetic_communities_dir),
        # From create_HoMi_metadata
        os.path.join(synthetic_work_dir, "synthetic_homi_metadata.csv"),
        # From run_HoMi_synthetic_communities
        "HoMi_is_done_synthetic"


rule simulate_synthetic_communities:
    input:
        sample_data=os.path.join(synthetic_work_dir, "sample_data.csv")
    output:
        communities=directory(os.path.join(synthetic_work_dir, synthetic_communities_dir)),
        done="communities_created"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(1*60) # min
    params:
        script=os.path.join(synthetic_work_dir, "create_mock_community.py"),
        work_dir=synthetic_work_dir,
        communities_dir=synthetic_communities_dir
    shell:
        """
        python {params.script} {input.sample_data} --work_dir {params.work_dir} --output_dir {params.communities_dir}
        touch {ouput.done}
        """


rule create_HoMi_metadata:
    input:
        sample_data=os.path.join(synthetic_work_dir, "sample_data.csv")
    output:
        homi_metadata=os.path.join(synthetic_work_dir, "synthetic_homi_metadata.csv")
    threads: 1
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB
        runtime=int(10) # min
    params:
        work_dir=synthetic_work_dir,
        communities_dir=synthetic_communities_dir
    run:
        import pandas as pd
        df = pd.read_csv(input.sample_data)
        genome_names = df["genome"].to_list()
        df = df.drop(["genome", "GCF_id"], axis=1).transpose()
        df.columns = genome_names

        df["forward_reads"] = [os.path.join(params.work_dir, params.communities_dir, f"{sample}_R1.fastq.gz") 
                                for sample in df.index]
        df["reverse_reads"] = [os.path.join(params.work_dir, params.communities_dir, f"{sample}_R1.fastq.gz") 
                                for sample in df.index]

        df.to_csv(output.homi_metadata, index_label="Sample")

rule run_HoMi_synthetic_communities:
    input:
        homi_metadata=os.path.join(synthetic_work_dir, "synthetic_homi_metadata.csv"),
        homi_config=os.path.join(synthetic_work_dir, "synthetic_HoMi_config.yaml"),
        communities_created="communities_created"  
    output:
        "HoMi_is_done_synthetic"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB
        runtime=int(20*60) # min
    params:
        homi_args=homi_args
    shell:
        """
        HoMi.py {input.homi_config} {params.homi_args}
        touch {output}
        """

