import os
from semisynthetic import simulate_semisynthetic as sss


# some global vars here for now BUT NEED TO BE MIGRATED TO CONFIG
synthetic_work_dir = "synthetic"
synthetic_communities_dir = "synthetic_communities"

#colon_sample_htx="SRP127360"
#colon_sample_htx="SRR6410603"


rule all:
    input:
        os.path.join(synthetic_work_dir, synthetic_communities_dir),
        os.path.join(synthetic_work_dir, "synthetic_homi_metadata.csv")


rule simulate_synthetic_communities:
    input:
        sample_data=os.path.join(synthetic_work_dir, "sample_data.csv")
    output:
        directory(os.path.join(synthetic_work_dir, synthetic_communities_dir))
    threads: 1
    params:
        script=os.path.join(synthetic_work_dir, "create_mock_community.py"),
        work_dir=synthetic_work_dir,
        communities_dir=synthetic_communities_dir
    shell:
        """
        python {params.script} {input.sample_data} --work_dir {params.work_dir} --output_dir {params.communities_dir}
        """


rule create_HoMi_metadata:
    input:
        sample_data=os.path.join(synthetic_work_dir, "sample_data.csv")
    output:
        homi_metadata=os.path.join(synthetic_work_dir, "synthetic_homi_metadata.csv")
    threads: 1
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

        df.to_csv(output.homi_metadata, index_label="SampleID")

# rule run_HoMi_synthetic_communities:
