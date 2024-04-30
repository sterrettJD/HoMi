import os
from semisynthetic import simulate_semisynthetic as sss


# some global vars here for now BUT NEED TO BE MIGRATED TO CONFIG
synthetic_work_dir = "synthetic"
synthetic_communities_dir = "synthetic_communities"

#colon_sample_htx="SRP127360"
#colon_sample_htx="SRR6410603"


rule all:
    input:
        os.path.join(synthetic_work_dir, synthetic_communities_dir)


rule simulate_synthetic_communities:
    input:
        sample_data=os.path.join(synthetic_work_dir, "sample_data.csv")
    output:
        directory(os.path.join(synthetic_work_dir, synthetic_communities_dir))
    threads: 1
    params:
        script=os.path.join(synthetic_work_dir, "create_mock_community.py"),
        
    shell:
        """
        python {params.script} {input.sample_data} --work_dir synthetic --output_dir synthetic_communities
        """