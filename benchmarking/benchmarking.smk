import os
import pandas as pd
from semisynthetic import simulate_semisynthetic as sss


# some global vars here for now BUT NEED TO BE MIGRATED TO CONFIG
synthetic_work_dir = "synthetic"
synthetic_communities_dir = "synthetic_communities"
synthetic_transcriptomes_dir = "synthetic_transcriptomes"
homi_args = "--profile slurm"
metadata_file = os.path.join(synthetic_work_dir, "sample_data.csv")
metadata = pd.read_csv(metadata_file)
samples = metadata.drop(columns=["genome", "GCF_id"]).columns
reads = ["R1", "R2"]
pereira_srr_ids = ["SRR25593833"]


#colon_sample_htx="SRP127360"
#colon_sample_htx="SRR6410603"


rule all:
    input:
        # From simulate_synthetic_communities
        os.path.join(synthetic_work_dir, synthetic_communities_dir),
        # From create_HoMi_metadata
        os.path.join(synthetic_work_dir, "synthetic_homi_metadata.csv"),
        # From run_HoMi_synthetic_communities
        "HoMi_is_done_synthetic",
        "synthetic_transcriptomes_created",
        # From synthetic transcriptomes
        expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_{read}.fastq"),
               sample=samples,
               read=reads),
        expand(os.path.join("Pereira", "{srr_id}_{lane}.fastq"),
                srr_id=pereira_srr_ids,
                lane=[1,2])


rule simulate_synthetic_communities:
    input:
        sample_data=metadata_file
    output:
        communities=directory(os.path.join(synthetic_work_dir, synthetic_communities_dir)),
        done="communities_created"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(2*60) # min
    params:
        script=os.path.join(synthetic_work_dir, "create_mock_community.py"),
        work_dir=synthetic_work_dir,
        communities_dir=synthetic_communities_dir
    shell:
        """
        python {params.script} {input.sample_data} --work_dir {params.work_dir} --output_dir {params.communities_dir}
        touch {output.done}
        """


rule simulate_synthetic_host_transcriptomes:
    input:
        sample_data=metadata_file
    output:
        data=expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_{read}.fasta"),
                    sample=samples,
                    read=reads),
        done="synthetic_transcriptomes_created"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(12*60) # min
    params:
        script=os.path.join(synthetic_work_dir, "run_polyester.R"),
        work_dir=synthetic_work_dir,
        communities_dir=os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir)
    shell:
        """
        Rscript {params.script} \
        -t synthetic/data/host_transcriptome.fna.gz \
        --transcriptome_url https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz \
        -g synthetic/data/host_transcriptome.gff.gz \
        --gtf_url https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.gff.gz \
        -s {input.sample_data} \
        -n human \
        -o {params.communities_dir}
        
        touch {output.done}
        """


rule transcriptome_fasta_to_fastq:
    input:
        data=os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_{read}.fasta")
    output:
        data=os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_{read}.fastq")
    threads: 1
    conda: "../conda_envs/bbmap.yaml"
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(4*60) # min
    shell:
        """
        reformat.sh in={input.data} out={output.data} qin=33 qout=33 qfake=40
        rm {input.data}
        """


rule create_HoMi_metadata_synthetic:
    input:
        sample_data=metadata_file
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


rule prefetch_Pereira:
    output:
        directory(os.path.join("Pereira", "{srr_id}"))
    conda: "sra_tools.yaml"
    threads: 8
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB
        runtime=int(2*60) # min
    shell:
        """
        mkdir -p Pereira
        cd Pereira
        prefetch {wildcards.srr_id}
        cd ..
        """

rule fasterq_dump_Pereira:
    input:
        os.path.join("Pereira", "{srr_id}")
    output:
        os.path.join("Pereira", "{srr_id}_{lane}.fastq")
    conda: "sra_tools.yaml"
    threads: 8
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB
        runtime=int(1*60) # min
    shell:
        """
        cd Pereira
        fasterq-dump -p {wildcards.srr_id} -e {threads}
        """
