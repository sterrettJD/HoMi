import os
import pandas as pd


######## CONFIG ########
# some global vars here
synthetic_work_dir = "synthetic"
synthetic_communities_dir = "synthetic_communities"
synthetic_transcriptomes_dir = "synthetic_transcriptomes"

# to run this on a Slurm-managed cluster
homi_args = "--profile slurm"

# Metadata
metadata_file = os.path.join(synthetic_work_dir, "sample_data.csv")
metadata = pd.read_csv(metadata_file)
samples = metadata.drop(columns=["genome", "GCF_id"]).columns
reads = ["R1", "R2"]
organisms = metadata["genome"].to_list()
microbial_organisms = [x for x in organisms if (x != "human")]


# mock community data
pereira_df = pd.read_csv("Pereira/Pereira_data.csv")
pereira_srr_ids = pereira_df["SRR"]

# Per nucleotide quality score for Polyester-simulated reads
polyester_error_rate=0.001
polyester_phred=30

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
        
        # From synthetic transcriptomes
        expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_R1.fastq.gz"),
               sample=samples),
        expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_R2.fastq.gz"),
               sample=samples),
        "HoMi_is_done_synthetic_transcriptomes",
        "synthetic_transcriptomes_benchmark.pdf",
        "synthetic_transcriptomes_benchmark_lm_results.txt",

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
        data=expand(os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_human", "{sample}_unsampled_{read}.fasta"),
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
        communities_dir=os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir),
        polyester_error_rate=polyester_error_rate
    shell:
        """
        Rscript {params.script} \
        -t synthetic/data/host_transcriptome.fna.gz \
        --transcriptome_url https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz \
        -g synthetic/data/host_transcriptome.gff.gz \
        --gtf_url https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.gff.gz \
        -e {params.polyester_error_rate} \
        -s {input.sample_data} \
        -n human \
        -o {params.communities_dir}_human
        
        touch {output.done}
        """


rule transcriptome_fasta_to_fastq:
    input:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_human", "{sample}_unsampled_{read}.fasta")
    output:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_human", "{sample}_unsampled_{read}.fastq")
    threads: 1
    conda: "../conda_envs/bbmap.yaml"
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(4*60) # min
    params:
        polyester_phred=polyester_phred
    shell:
        """
        reformat.sh in={input.data} out={output.data} qin={params.polyester_phred} qout={params.polyester_phred} qfake={params.polyester_phred}
        rm {input.data}
        """


rule subsample_fastq_to_correct_depth:
    input:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_human", "{sample}_unsampled_{read}.fastq")
    output:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_human", "{sample}_{read}.fastq")
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
        

rule simulate_synthetic_microbial_transcriptomes:
    input:
        sample_data=metadata_file
    output:
        done="synthetic_microbial_transcriptomes_created_{organism}"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(12*60) # min
    params:
        script=os.path.join(synthetic_work_dir, "run_polyester.R"),
        work_dir=synthetic_work_dir,
        communities_dir=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_{{organism}}_u"),
        polyester_error_rate=polyester_error_rate
    shell:
        """
        Rscript {params.script} \
        -t synthetic/data/{wildcards.organism}/genome/cds_from_genomic.fna \
        -g synthetic/data/{wildcards.organism}/genome/genomic.gff \
        -e {params.polyester_error_rate} \
        -s {input.sample_data} \
        -n {wildcards.organism} \
        -o {params.communities_dir}
        
        touch {output.done}
        """


rule transcriptome_fasta_to_fastq_microbial:
    input:
        data_created="synthetic_microbial_transcriptomes_created_{organism}",
    output:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_{{organism}}_u", "{sample}_unsampled_{read}.fastq")
    threads: 1
    conda: "../conda_envs/bbmap.yaml"
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(4*60) # min
    params:
        in_data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_{{organism}}_u", "{sample}_unsampled_{read}.fasta"),
        polyester_phred=polyester_phred
    shell:
        """
        reformat.sh in={params.in_data} out={output.data} qin={params.polyester_phred} qout={params.polyester_phred} qfake={params.polyester_phred}
        rm {params.in_data}
        """


rule subsample_fastq_to_correct_depth_microbial:
    input:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_{{organism}}_u", "{sample}_unsampled_{read}.fastq")
    output:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_{{organism}}_s", "{sample}_{read}.fastq")
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
        depth = metadata.loc[wildcards.organism, wildcards.sample]
        if depth > 0:
            cmd = f"seqtk sample -s 123 {input.data} {depth} > {output.data}"
        else:
            cmd = f"cp {input.data} {output.data}"
        ran = subprocess.run(cmd, shell=True)
        cleaned = subprocess.run(["rm", input.data])


rule combine_transcriptomes:
    input:
        expand(os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_{{organism}}_s", "{sample}_{read}.fastq"),
               organism=microbial_organisms, sample=samples, read=reads)
    output:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}", "{sample}_{read}.fastq.gz")
    threads: 1
    resources:
        partition="short",
        mem_mb=int(2*1000), # MB
        runtime=int(1*60) # min
    params:
        organisms=organisms,
        synthetic_work_dir=synthetic_work_dir,
        synthetic_transcriptomes_dir=synthetic_transcriptomes_dir
    run:
        import subprocess
        import os
	
        out_dir = os.path.join(params.synthetic_work_dir, f"{params.synthetic_transcriptomes_dir}")
        os.makedirs(out_dir, exist_ok=True)

        # Find the paths to each organism's transcriptome
        in_paths = [os.path.join(params.synthetic_work_dir, 
                                f"{params.synthetic_transcriptomes_dir}_{organism}_s", 
                                f"{wildcards.sample}_{wildcards.read}.fastq") 
                    for organism in params.organisms]
        
        in_paths_string = " ".join(in_paths)
        # Human reads don't need the _s/_u to make them distinct for snakemake, so fixing that here.
        # Could fix that in its rule, but I would need to rerun it which takes time
        in_paths_string = in_paths_string.replace("human_s", "human")

        cmd = f"cat {in_paths_string} | gzip > {output.data}"
        ran = subprocess.run(cmd, shell=True)


###################################
### HoMi on synthetic ###

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
        df["reverse_reads"] = [os.path.join(params.work_dir, params.communities_dir, f"{sample}_R2.fastq.gz") 
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
        data="benchmarking_synthetic_reads_breakdown.csv",
        label="True percent host reads (pangenome, jittered)" 
    shell:
        """
        Rscript {params.script} -i {params.data} -o {output.plot} -n {params.label} > {output.model}
        """


rule create_HoMi_metadata_synthetic_transcriptomes:
    input:
        sample_data=metadata_file
    output:
        homi_metadata=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_homi_metadata.csv")
    threads: 1
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB
        runtime=int(10) # min
    params:
        work_dir=synthetic_work_dir,
        communities_dir=synthetic_transcriptomes_dir
    run:
        import pandas as pd
        df = pd.read_csv(input.sample_data)
        genome_names = df["genome"].to_list()
        df = df.drop(["genome", "GCF_id"], axis=1).transpose()
        df.columns = genome_names

        df["forward_reads"] = [os.path.join(params.work_dir, params.communities_dir, f"{sample}_R1.fastq.gz") 
                                for sample in df.index]
        df["reverse_reads"] = [os.path.join(params.work_dir, params.communities_dir, f"{sample}_R2.fastq.gz") 
                                for sample in df.index]

        df.to_csv(output.homi_metadata, index_label="Sample")
    

rule run_HoMi_synthetic_transcriptomes:
    input:
        homi_metadata=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_homi_metadata.csv"),
        homi_config=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_HoMi_config.yaml"),
        fwd=expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_R1.fastq.gz"),
               sample=samples),
        rev=expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_R2.fastq.gz"),
               sample=samples)
    output:
        "HoMi_is_done_synthetic_transcriptomes"
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

rule plot_expected_vs_actual_synthetic_transcriptomes:
    input:
        "HoMi_is_done_synthetic_transcriptomes"
    output:
        plot="synthetic_transcriptomes_benchmark.pdf",
        model="synthetic_transcriptomes_benchmark_lm_results.txt"
    conda: "conda_envs/r_env.yaml"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(4*1000), # MB
        runtime=int(1*60) # min
    params:
        script="Plot_benchmarked_reads_breakdown.R",
        data="benchmarking_synthetic_transcriptomes_reads_breakdown.csv",
        label="True percent host reads (GRCh38, jittered)"
    shell:
        """
        Rscript {params.script} -i {params.data} -o {output.plot}  -n {params.label} > {output.model}
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
        jitter=0,
        label="Pereira-Marques percent microbial reads"
    shell:
        """
        Rscript {params.script} -i {params.data} -j {params.jitter} -o {output.plot} -n {params.label} > {output.model}
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
