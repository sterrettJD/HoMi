import os
import pandas as pd
from semisynthetic import simulate_semisynthetic as sss


# some global vars here for now BUT SHOULD BE MIGRATED TO CONFIG
synthetic_work_dir = "synthetic"
synthetic_communities_dir = "synthetic_communities"
synthetic_transcriptomes_dir = "synthetic_transcriptomes"
homi_args = "--profile slurm"
metadata_file = os.path.join(synthetic_work_dir, "sample_data.csv")
metadata = pd.read_csv(metadata_file)
samples = metadata.drop(columns=["genome", "GCF_id"]).columns
reads = ["R1", "R2"]
# mock community data
pereira_df = pd.read_csv("Pereira/Pereira_data.csv")
pereira_srr_ids = pereira_df["SRR"]


#colon_sample_htx="SRP127360"
#colon_sample_htx="SRR6410603"


rule all:
    input:
        # From simulate_synthetic_communities
        expand(os.path.join(synthetic_work_dir, synthetic_communities_dir, "{sample}_R1.fastq.gz"),
               sample=samples),
        expand(os.path.join(synthetic_work_dir, synthetic_communities_dir, "{sample}_R2.fastq.gz"),
               sample=samples),
        # From create_HoMi_metadata
        os.path.join(synthetic_work_dir, "synthetic_homi_metadata.csv"),
        # From run_HoMi_synthetic_communities
        "HoMi_is_done_synthetic",
        "synthetic_communities_benchmark.pdf",
        "synthetic_communities_benchmark_lm_results.txt",

        "synthetic_transcriptomes_created",
        
        # From synthetic transcriptomes
        expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_R1.fastq"),
               sample=samples),
        expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_R2.fastq"),
               sample=samples),

        # From mock communities
        expand(os.path.join("Pereira", "{srr_id}_R1.fastq.gz"),
                srr_id=pereira_srr_ids),
        expand(os.path.join("Pereira", "{srr_id}_R2.fastq.gz"),
                srr_id=pereira_srr_ids),
        "HoMi_is_done_Pereira",
        "Pereira_benchmark.pdf",
        "Pereira_benchmark_lm_results.txt",
        "Pereira_benchmark_from_paper.pdf",
        "Pereira_benchmark_from_paper_lm_results.txt"


rule pull_reference_genomes:
    input:
        sample_data=metadata_file
    output:
        done="reference_genomes_downloaded"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(4*1000), # MB
        runtime=int(2*60) # min
    params:
        script=os.path.join(synthetic_work_dir, "pull_reference_genomes.py"),
        work_dir=synthetic_work_dir,
        communities_dir=synthetic_communities_dir
    shell:
        """
        python {params.script} {input.sample_data} --work_dir {params.work_dir}
        touch {output.done}
        """


rule simulate_synthetic_communities:
    input:
        sample_data=metadata_file,
        references_downloaded="reference_genomes_downloaded"
    output:
        fwd=os.path.join(synthetic_work_dir, synthetic_communities_dir, "{sample}_R1.fastq.gz"),
        rev=os.path.join(synthetic_work_dir, synthetic_communities_dir, "{sample}_R2.fastq.gz")
    threads: 1
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(4*60) # min
    params:
        script=os.path.join(synthetic_work_dir, "create_mock_community.py"),
        work_dir=synthetic_work_dir,
        communities_dir=synthetic_communities_dir
    shell:
        """
        python {params.script} {input.sample_data} {wildcards.sample} --work_dir {params.work_dir} --output_dir {params.communities_dir}
        """


rule simulate_synthetic_host_transcriptomes:
    input:
        sample_data=metadata_file
    output:
        data=expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_unsampled_{read}.fasta"),
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
        data=os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_unsampled_{read}.fasta")
    output:
        data=os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_unsampled_{read}.fastq")
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


rule subsample_fastq_to_correct_depth:
    input:
        data=os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_unsampled_{read}.fastq")
    output:
        data=os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_{read}.fastq")
    threads: 1
    resources:
        partition="short",
        mem_mb=int(2*1000), # MB
        runtime=int(1*60) # min
    params:
        metadata=metadata_file
    run:
        import subprocess
        import pandas as pd
        metadata = pd.read_csv(params.metadata, index_col="genome")
        depth = metadata.loc["human", wildcards.sample]
        if depth > 0:
            cmd = f"seqtk sample -s 123 {input.data} {depth} > {output.data}"
        else:
            cmd = f"cp {input.data} {output.data}"
        ran = subprocess.run(cmd, shell=True)
        cleaned = subprocess.run(["rm", input.data])
        

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
        fwd=expand(os.path.join(synthetic_work_dir, synthetic_communities_dir, "{sample}_R1.fastq.gz"),
               sample=samples),
        rev=expand(os.path.join(synthetic_work_dir, synthetic_communities_dir, "{sample}_R2.fastq.gz"),
               sample=samples)
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
        HoMi.py {input.homi_config} {params.homi_args} --unlock
        touch {output}
        """

rule plot_expected_vs_actual_synthetic_communities:
    input:
        "HoMi_is_done_synthetic"
    output:
        plot="synthetic_communities_benchmark.pdf",
        model="synthetic_communities_benchmark_lm_results.txt"
    conda: "conda_envs/r_env.yaml"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(4*1000), # MB
        runtime=int(1*60) # min
    params:
        script="Plot_benchmarked_reads_breakdown.R",
        data="benchmarking_synthetic_reads_breakdown.csv"
    shell:
        """
        Rscript {params.script} -i {params.data} -o {output.plot} > {output.model}
        """


#############################################################
##### Mock communities #####
rule fastq_dump_Pereira:
    output:
        fwd=os.path.join("Pereira", "{srr_id}_R1.fastq.gz"),
        rev=os.path.join("Pereira", "{srr_id}_R2.fastq.gz")
    conda: "conda_envs/sra_tools.yaml"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB
        runtime=int(4*60) # min
    shell:
        """
        mkdir -p Pereira
        cd Pereira
        fastq-dump --gzip --readids --read-filter pass --dumpbase --split-3 --clip {wildcards.srr_id}
        mv {wildcards.srr_id}_pass_1.fastq.gz > {output.fwd}
        mv {wildcards.srr_id}_pass_2.fastq.gz > {output.rev}
        """

rule run_HoMi_mock_data:
    input:
        homi_metadata=os.path.join("Pereira", "Pereira_data.csv"),
        homi_config=os.path.join("Pereira", "mock_community_HoMi_config.yaml"),
        fwd=expand(os.path.join("Pereira", "{srr_id}_R1.fastq.gz"),
                srr_id=pereira_srr_ids),
        rev=expand(os.path.join("Pereira", "{srr_id}_R2.fastq.gz"),
                srr_id=pereira_srr_ids)
    output:
        "HoMi_is_done_Pereira"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB
        runtime=int(24*60) # min
    params:
        homi_args=homi_args
    shell:
        """
        HoMi.py {input.homi_config} {params.homi_args} --unlock
        touch {output}
        """


rule plot_expected_vs_actual_mock_data:
    input:
        "HoMi_is_done_Pereira"
    output:
        plot="Pereira_benchmark.pdf",
        model="Pereira_benchmark_lm_results.txt"
    conda: "conda_envs/r_env.yaml"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(4*1000), # MB
        runtime=int(1*60) # min
    params:
        script="Plot_benchmarked_reads_breakdown.R",
        data="benchmarking_Pereira_reads_breakdown.csv",
        jitter=0
    shell:
        """
        Rscript {params.script} -i {params.data} -j {params.jitter} -o {output.plot} > {output.model}
        """

rule plot_expected_from_paper_vs_actual_mock_data:
    input:
        "HoMi_is_done_Pereira"
    output:
        plot="Pereira_benchmark_from_paper.pdf",
        model="Pereira_benchmark_from_paper_lm_results.txt"
    conda: "conda_envs/r_env.yaml"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(4*1000), # MB
        runtime=int(1*60) # min
    params:
        script="Plot_benchmarked_reads_breakdown.R",
        data="benchmarking_Pereira_reads_breakdown.csv",
        metadata="Pereira/Pereira_data.csv",
        column_name="Pereira_percent_microbial",
        axis_name="\"Paper-derived percent microbial\"",
        jitter=0
    shell:
        """
        Rscript {params.script} -i {params.data} -m {params.metadata} \
        -c {params.column_name} -n {params.axis_name} \
        -j {params.jitter} \
        -o {output.plot} > {output.model}
        """
