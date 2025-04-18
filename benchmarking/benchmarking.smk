import os
import pandas as pd


######## CONFIG ########
# some global vars here
synthetic_work_dir = "synthetic"
synthetic_communities_dir = "synthetic_communities"
synthetic_transcriptomes_dir = "synthetic_transcriptomes"
synthetic_transcriptomes_dir_p40 = "synthetic_transcriptomes_p40"
semi_work_dir = "semi"

# to run this on a Slurm-managed cluster
homi_args = "--profile slurm"

# Syntheti communities metadata
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
semi_homi_args = """--profile slurm --snakemake_extra "--jobs 40" """


# mock community data
pereira_df = pd.read_csv("Pereira/Pereira_data.csv")
pereira_srr_ids = pereira_df["SRR"]

# Per nucleotide quality score for Polyester-simulated reads
polyester_error_rate=0.001
polyester_phred=30

polyester_error_rate_p40=0.0001
polyester_phred_p40=40

# hostile reference data
t2t_rna_hla_index = "t2t_rna_hla_index"

indexes = ["dna", "rna"]

rule all:
    input:
        # hostile reference
        t2t_rna_hla_index,
        # From simulate_synthetic_communities
        expand(os.path.join(synthetic_work_dir, synthetic_communities_dir, "{sample}_R1.fastq.gz"),
               sample=samples),
        expand(os.path.join(synthetic_work_dir, synthetic_communities_dir, "{sample}_R2.fastq.gz"),
               sample=samples),
        # From create_HoMi_metadata
        os.path.join(synthetic_work_dir, "synthetic_homi_metadata.csv"),
        # From run_HoMi_synthetic_communities
        expand("{index}_synthetic_communities_benchmark.pdf",
            index=indexes),
        expand("{index}_synthetic_communities_benchmark_lm_results.txt",
            index=indexes),
        
        # From synthetic transcriptomes
        expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_R1.fastq.gz"),
               sample=samples),
        expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_R2.fastq.gz"),
               sample=samples),
        
         # From synthetic transcriptomes PHRED 40
        expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir_p40, "{sample}_R1.fastq.gz"),
               sample=samples),
        expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir_p40, "{sample}_R2.fastq.gz"),
               sample=samples),

        # Comparison of transcriptomes to actual
        expand("{index}_index_{proj}_benchmark.pdf",
                proj=["synthetic_transcriptomes", "synthetic_transcriptomes_p40"],
                index=indexes),
        expand("{index}_index_{proj}_benchmark_lm_results.txt",
                proj=["synthetic_transcriptomes", "synthetic_transcriptomes_p40"],
            index=indexes),

        # From mock communities
        expand(os.path.join("Pereira", "{srr_id}_R1.fastq.gz"),
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
            index=indexes),

        # Boxplot comparisons of taxonomy to what it should be
        expand("taxonomy_compared/{index}/{proj}/{method}_genus.pdf",
               index=indexes,
               proj=["synthetic_transcriptomes", "synthetic_transcriptomes_p40", "synthetic"],
               method=["kraken", "metaphlan"]),
        expand("taxonomy_compared/{index}/{proj}/{method}_species.pdf",
               index=indexes,
               proj=["synthetic_transcriptomes", "synthetic_transcriptomes_p40", "synthetic"],
               method=["kraken", "metaphlan"]),
        "taxonomy_compared/combined_taxa_boxplot_genus.pdf",
        "taxonomy_compared/combined_taxa_boxplot_species.pdf",

        # Semisynthetic data files
        expand(os.path.join("semi", "samples", "{sample}_{read}.fastq.gz"),
               sample=semi_samples, read=reads),
        expand("{index}_semi_benchmark.pdf",
                index=indexes),
        expand("{index}_semi_benchmark_lm_results.txt",
                index=indexes)


rule create_alt_hostile_index:
    output:
        ref_dir=directory(t2t_rna_hla_index)
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
        python {params.script} -o {output.ref_dir}
        """


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

# Simulating with custom script from human pangenome project
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


####################################################
### Simulate transcriptomes with phred of 30
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
        --error_rate {params.polyester_error_rate} \
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
    conda: "conda_envs/bbmap.yaml"
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(4*60) # min
    params:
        polyester_phred=polyester_phred
    shell:
        """
        reformat.sh in={input.data} out={output.data} qin=33 qout=33 qfake={params.polyester_phred}
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
        --error_rate {params.polyester_error_rate} \
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
    conda: "conda_envs/bbmap.yaml"
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(4*60) # min
    params:
        in_data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_{{organism}}_u", "{sample}_unsampled_{read}.fasta"),
        polyester_phred=polyester_phred
    shell:
        """
        reformat.sh in={params.in_data} out={output.data} qin=33 qout=33 qfake={params.polyester_phred}
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
        microbes=expand(os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_{{organism}}_s", "{sample}_{read}.fastq"),
               organism=microbial_organisms, sample=samples, read=reads),
        human=expand(os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir}_human", "{sample}_{read}.fastq"),
               sample=samples, read=reads)
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


####################################################
### Simulate transcriptomes with phred of 40
# Could be done with wildcards on the other transcriptome simulations
# But doing it separately in the interest of time and not rerunning those simualtions

rule simulate_synthetic_host_transcriptome_p40:
    input:
        sample_data=metadata_file
    output:
        data=expand(os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir_p40}_human", "{sample}_unsampled_{read}.fasta"),
                    sample=samples,
                    read=reads),
        done="synthetic_transcriptomes_created_p40"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(12*60) # min
    params:
        script=os.path.join(synthetic_work_dir, "run_polyester.R"),
        work_dir=synthetic_work_dir,
        communities_dir=os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir_p40),
        polyester_error_rate=polyester_error_rate_p40
    shell:
        """
        Rscript {params.script} \
        -t synthetic/data/host_transcriptome.fna.gz \
        --transcriptome_url https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz \
        -g synthetic/data/host_transcriptome.gff.gz \
        --gtf_url https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.gff.gz \
        --error_rate {params.polyester_error_rate} \
        -s {input.sample_data} \
        -n human \
        -o {params.communities_dir}_human
        
        touch {output.done}
        """


rule transcriptome_fasta_to_fastq_p40:
    input:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir_p40}_human", "{sample}_unsampled_{read}.fasta")
    output:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir_p40}_human", "{sample}_unsampled_{read}.fastq")
    threads: 1
    conda: "conda_envs/bbmap.yaml"
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(4*60) # min
    params:
        polyester_phred=polyester_phred_p40
    shell:
        """
        reformat.sh in={input.data} out={output.data} qin=33 qout=33 qfake={params.polyester_phred}
        rm {input.data}
        """


rule subsample_fastq_to_correct_depth_p40:
    input:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir_p40}_human", "{sample}_unsampled_{read}.fastq")
    output:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir_p40}_human", "{sample}_{read}.fastq")
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
        

rule simulate_synthetic_microbial_transcriptomes_p40:
    input:
        sample_data=metadata_file
    output:
        done="synthetic_microbial_transcriptomes_p40_created_{organism}"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(12*60) # min
    params:
        script=os.path.join(synthetic_work_dir, "run_polyester.R"),
        work_dir=synthetic_work_dir,
        communities_dir=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir_p40}_{{organism}}_u_p40"),
        polyester_error_rate=polyester_error_rate_p40
    shell:
        """
        Rscript {params.script} \
        -t synthetic/data/{wildcards.organism}/genome/cds_from_genomic.fna \
        -g synthetic/data/{wildcards.organism}/genome/genomic.gff \
        --error_rate {params.polyester_error_rate} \
        -s {input.sample_data} \
        -n {wildcards.organism} \
        -o {params.communities_dir}
        
        touch {output.done}
        """


rule transcriptome_fasta_to_fastq_microbial_p40:
    input:
        data_created="synthetic_microbial_transcriptomes_p40_created_{organism}",
    output:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir_p40}_{{organism}}_u_p40", "{sample}_unsampled_{read}.fastq")
    threads: 1
    conda: "conda_envs/bbmap.yaml"
    resources:
        partition="short",
        mem_mb=int(16*1000), # MB
        runtime=int(4*60) # min
    params:
        in_data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir_p40}_{{organism}}_u_p40", "{sample}_unsampled_{read}.fasta"),
        polyester_phred=polyester_phred_p40
    shell:
        """
        reformat.sh in={params.in_data} out={output.data} qin=33 qout=33 qfake={params.polyester_phred}
        rm {params.in_data}
        """


rule subsample_fastq_to_correct_depth_microbial_p40:
    input:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir_p40}_{{organism}}_u_p40", "{sample}_unsampled_{read}.fastq")
    output:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir_p40}_{{organism}}_s_p40", "{sample}_{read}.fastq")
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


rule combine_transcriptomes_p40:
    input:
        microbes=expand(os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir_p40}_{{organism}}_s_p40", "{sample}_{read}.fastq"),
               organism=microbial_organisms, sample=samples, read=reads),
        host=expand(os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir_p40}_human", "{sample}_{read}.fastq"),
               sample=samples, read=reads)
    output:
        data=os.path.join(synthetic_work_dir, f"{synthetic_transcriptomes_dir_p40}", "{sample}_{read}.fastq.gz")
    threads: 1
    resources:
        partition="short",
        mem_mb=int(2*1000), # MB
        runtime=int(1*60) # min
    params:
        organisms=organisms,
        synthetic_work_dir=synthetic_work_dir,
        synthetic_transcriptomes_dir_p40=synthetic_transcriptomes_dir_p40
    run:
        import subprocess
        import os
	
        out_dir = os.path.join(params.synthetic_work_dir, f"{params.synthetic_transcriptomes_dir_p40}")
        os.makedirs(out_dir, exist_ok=True)

        # Find the paths to each organism's transcriptome
        in_paths = [os.path.join(params.synthetic_work_dir, 
                                f"{params.synthetic_transcriptomes_dir_p40}_{organism}_s_p40", 
                                f"{wildcards.sample}_{wildcards.read}.fastq") 
                    for organism in params.organisms]
        
        in_paths_string = " ".join(in_paths)
        # Human reads don't need the _s/_u to make them distinct for snakemake, so fixing that here.
        # Could fix that in its rule, but I would need to rerun it which takes time
        in_paths_string = in_paths_string.replace("human_s_p40", "human")

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
        homi_config=os.path.join(synthetic_work_dir, "{index}_synthetic_HoMi_config.yaml"),
        fwd=expand(os.path.join(synthetic_work_dir, synthetic_communities_dir, "{sample}_R1.fastq.gz"),
               sample=samples),
        rev=expand(os.path.join(synthetic_work_dir, synthetic_communities_dir, "{sample}_R2.fastq.gz"),
               sample=samples)
    output:
        "{index}_HoMi_is_done_synthetic"
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
        "{index}_HoMi_is_done_synthetic"
    output:
        plot="{index}_synthetic_communities_benchmark.pdf",
        model="{index}_synthetic_communities_benchmark_lm_results.txt"
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
        Rscript {params.script} -i {wildcards.index}_benchmarking_synthetic_reads_breakdown.csv -o {output.plot} -n "{params.label}" > {output.model}
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
        homi_config=os.path.join(synthetic_work_dir, "{index}_synthetic_transcriptomes_HoMi_config.yaml"),
        fwd=expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_R1.fastq.gz"),
               sample=samples),
        rev=expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir, "{sample}_R2.fastq.gz"),
               sample=samples)
    output:
        "{index}_HoMi_is_done_synthetic_transcriptomes"
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


rule create_HoMi_metadata_synthetic_transcriptomes_p40:
    input:
        sample_data=metadata_file
    output:
        homi_metadata=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p40_homi_metadata.csv")
    threads: 1
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB
        runtime=int(10) # min
    params:
        work_dir=synthetic_work_dir,
        communities_dir=synthetic_transcriptomes_dir_p40
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
    

rule run_HoMi_synthetic_transcriptomes_p40:
    input:
        homi_metadata=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p40_homi_metadata.csv"),
        homi_config=os.path.join(synthetic_work_dir, "{index}_synthetic_transcriptomes_p40_HoMi_config.yaml"),
        fwd=expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir_p40, "{sample}_R1.fastq.gz"),
               sample=samples),
        rev=expand(os.path.join(synthetic_work_dir, synthetic_transcriptomes_dir_p40, "{sample}_R2.fastq.gz"),
               sample=samples)
    output:
        "{index}_HoMi_is_done_synthetic_transcriptomes_p40"
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
        "{index}_HoMi_is_done_{proj}"
    output:
        plot="{index}_index_{proj}_benchmark.pdf",
        model="{index}_index_{proj}_benchmark_lm_results.txt"
    conda: "conda_envs/r_env.yaml"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(4*1000), # MB
        runtime=int(1*60) # min
    params:
        script="Plot_benchmarked_reads_breakdown.R",
        label="True percent host reads (GRCh38, jittered)"
    shell:
        """
        Rscript {params.script} -i {wildcards.index}_benchmarking_{wildcards.proj}_reads_breakdown.csv -o {output.plot}  -n "{params.label}" > {output.model}
        """


rule plot_taxonomy_boxplots:
    input:
        "{index}_HoMi_is_done_{proj}"
    output:
        genus_plot="taxonomy_compared/{index}/{proj}/{method}_genus.pdf",
        species_plot="taxonomy_compared/{index}/{proj}/{method}_species.pdf"
    conda: "conda_envs/r_env.yaml"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(4*1000), # MB
        runtime=int(1*60) # min
    params:
        script="Compare_taxonomy.R",
        outdir="taxonomy_compared/{index}/{proj}"
    shell:
        """
        mkdir -p {params.outdir}
        
        if [[ "{wildcards.method}" == "metaphlan" ]]; then
            Rscript {params.script} -i {wildcards.index}_benchmarking_{wildcards.proj}.f0.0.r0.0.nonhost.humann/all_bugs_list.tsv -t metaphlan -o {params.outdir}
        
        elif [[ "{wildcards.method}" == "kraken" ]]; then
            Rscript {params.script} -i {wildcards.index}_benchmarking_{wildcards.proj}.f0.0.r0.0.nonhost.kraken/Combined-taxonomy.tsv -t kraken -o {params.outdir}
        fi
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
    params:
        work_dir="Pereira"
    shell:
        """
        mkdir -p Pereira
        cd Pereira
        fastq-dump --gzip --readids --read-filter pass --dumpbase --split-3 --clip {wildcards.srr_id}
        
        cd ..
        mv {params.work_dir}/{wildcards.srr_id}_pass_1.fastq.gz {output.fwd}
        mv {params.work_dir}/{wildcards.srr_id}_pass_2.fastq.gz {output.rev}
        """

rule run_HoMi_mock_data:
    input:
        homi_metadata=os.path.join("Pereira", "Pereira_data.csv"),
        homi_config=os.path.join("Pereira", "{index}_mock_community_HoMi_config.yaml"),
        fwd=expand(os.path.join("Pereira", "{srr_id}_R1.fastq.gz"),
                srr_id=pereira_srr_ids),
        rev=expand(os.path.join("Pereira", "{srr_id}_R2.fastq.gz"),
                srr_id=pereira_srr_ids)
    output:
        "{index}_HoMi_is_done_Pereira"
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
        "{index}_HoMi_is_done_Pereira"
    output:
        plot="{index}_Pereira_benchmark.pdf",
        model="{index}_Pereira_benchmark_lm_results.txt"
    conda: "conda_envs/r_env.yaml"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(4*1000), # MB
        runtime=int(1*60) # min
    params:
        script="Plot_benchmarked_reads_breakdown.R",
        data="{index}_benchmarking_Pereira_reads_breakdown.csv",
        jitter=0,
        label="Pereira-Marques percent microbial reads"
    shell:
        """
        Rscript {params.script} -i {params.data} -j {params.jitter} \
        -o {output.plot} -n "{params.label}" --no_dotted_line > {output.model}
        """

rule plot_expected_from_paper_vs_actual_mock_data:
    input:
        "{index}_HoMi_is_done_Pereira"
    output:
        plot="{index}_Pereira_benchmark_from_paper.pdf",
        model="{index}_Pereira_benchmark_from_paper_lm_results.txt"
    conda: "conda_envs/r_env.yaml"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(4*1000), # MB
        runtime=int(1*60) # min
    params:
        script="Plot_benchmarked_reads_breakdown.R",
        data="{index}_benchmarking_Pereira_reads_breakdown.csv",
        metadata="Pereira/Pereira_data.csv",
        column_name="Pereira_percent_microbial",
        axis_name="\"Paper-derived percent microbial\"",
        jitter=0
    shell:
        """
        Rscript {params.script} -i {params.data} -m {params.metadata} \
        -c {params.column_name} -n {params.axis_name} \
        -j {params.jitter} \
        --no_dotted_line \
        -o {output.plot} > {output.model}
        """


#############################################################
##### Semisynthetic communities #####
rule fastq_dump_semi:
    output:
        fwd=os.path.join(semi_work_dir, "data", "raw_{srr_id}_R1.fastq.gz"),
        rev=os.path.join(semi_work_dir, "data", "raw_{srr_id}_R2.fastq.gz")
    conda: "conda_envs/sra_tools.yaml"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB
        runtime=int(4*60) # min
    params:
        semi_work_dir=semi_work_dir
    shell:
        """
        mkdir -p semi/data
        mkdir -p semi/samples
        cd semi
        fastq-dump --gzip --readids --read-filter pass --dumpbase --split-3 --clip {wildcards.srr_id}
        cd ..
        mv {params.semi_work_dir}/{wildcards.srr_id}_pass_1.fastq.gz {output.fwd}
        mv {params.semi_work_dir}/{wildcards.srr_id}_pass_2.fastq.gz {output.rev}
        """


rule combine_semi_srrs:
    input:
        data=expand(os.path.join(semi_work_dir, "data", "raw_{srr_id}_{read}.fastq.gz"),
            srr_id=semi_srr_ids, read=reads)
    output:
        data=os.path.join(semi_work_dir, "data", "{taxon}_{read}.fastq.gz")
    threads: 1
    resources:
        partition="short",
        mem_mb=int(2*1000), # MB
        runtime=int(1*60) # min
    params:
        metadata=semi_metadata_file,
        data_dir=os.path.join(semi_work_dir, "data")
    run:
        import subprocess
        import pandas as pd
        import os

        metadata = pd.read_csv(params.metadata, index_col="genome")
        
        # For this taxon, get the SRRs
        taxon_output_path = output.data
        srrs = metadata.loc[wildcards.taxon, "SRR"].split(".")

        # cat those respective files together
        for srr_id in srrs:
            srr_path = os.path.join(params.data_dir, f"raw_{srr_id}_{wildcards.read}.fastq.gz")
            if not os.path.exists(srr_path):
                raise FileNotFoundError(f"File not found: {srr_path}")

            # Adding to output file
            cmd = f"cat {srr_path} >> {output.data}"
            print(f"running command: {cmd}")
            ran = subprocess.run(cmd, shell=True)
            
            if ran.returncode != 0:
                raise RuntimeError(f"Command failed: {cmd}\nStderr: {ran.stderr.decode()}")


rule subsample_and_combine_semi_fastqs:
    input:
        data=expand(os.path.join(semi_work_dir, "data", "{taxon}_{read}.fastq.gz"),
            taxon=semi_organisms, read=reads)
    output:
        data=os.path.join(semi_work_dir, "samples", "{sample}_{read}.fastq.gz")
    threads: 1
    resources:
        partition="short",
        mem_mb=int(2*1000), # MB
        runtime=int(1*60) # min
    params:
        metadata=semi_metadata_file,
        data_dir=os.path.join(semi_work_dir, "data")
    run:
        import subprocess
        import pandas as pd
        import os

        metadata = pd.read_csv(params.metadata, index_col="genome")
        
        # get the random seed to use. Doing it based on the hash of sample ID so that
        # the fwd and rev reads per sample get the same seed
        sample_hash = hash(wildcards.sample)

        # For each taxon, subsample it and add it to the gzipped output file
        for taxon in metadata.index:    
            depth = metadata.loc[taxon, wildcards.sample]
            print(f"sampling {taxon} to {depth} reads")

            if depth > 0:
                unsampled_path = os.path.join(params.data_dir, f"{taxon}_{wildcards.read}.fastq.gz")
                if not os.path.exists(unsampled_path):
                    raise FileNotFoundError(f"File not found: {unsampled_path}")

                cmd = f"seqtk sample -s {sample_hash} {unsampled_path} {depth} | gzip >> {output.data}"
                print(f"running command: {cmd}")
                ran = subprocess.run(cmd, shell=True)
                
                if ran.returncode != 0:
                    raise RuntimeError(f"Command failed: {cmd}\nStderr: {ran.stderr.decode()}")



rule create_HoMi_metadata_semi:
    input:
        sample_data=semi_metadata_file
    output:
        homi_metadata=os.path.join(semi_work_dir, "semi_homi_metadata.csv")
    threads: 1
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB
        runtime=int(10) # min
    params:
        work_dir=semi_work_dir,
        communities_dir="samples"
    run:
        import pandas as pd 
        df = pd.read_csv(input.sample_data)
        genome_names = df["genome"].to_list()
        df = df.drop(["genome", "SRR"], axis=1).transpose()
        df.columns = genome_names

        df["forward_reads"] = [os.path.join(params.work_dir, params.communities_dir, f"{sample}_R1.fastq.gz") 
                                for sample in df.index]
        df["reverse_reads"] = [os.path.join(params.work_dir, params.communities_dir, f"{sample}_R2.fastq.gz") 
                                for sample in df.index]

        df.to_csv(output.homi_metadata, index_label="Sample")



rule run_HoMi_semi:
    input:
        homi_metadata=os.path.join(semi_work_dir, "semi_homi_metadata.csv"),
        homi_config=os.path.join(semi_work_dir, "{index}_semi_HoMi_config.yaml"),
        fwd=expand(os.path.join(semi_work_dir, "samples", "{sample}_{read}.fastq.gz"),
                sample=semi_samples, read=reads)
    output:
        "{index}_HoMi_is_done_semi"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(8*1000), # MB
        runtime=int(24*60) # min
    params:
        homi_args=semi_homi_args
    shell:
        """
        HoMi.py {input.homi_config} {params.homi_args} --unlock
        touch {output}
        """


rule plot_expected_vs_actual_semi:
    input:
        "{index}_HoMi_is_done_semi"
    output:
        plot="{index}_semi_benchmark.pdf",
        model="{index}_semi_benchmark_lm_results.txt"
    conda: "conda_envs/r_env.yaml"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(4*1000), # MB
        runtime=int(1*60) # min
    params:
        script="Plot_benchmarked_reads_breakdown.R",
        label="True percent host reads (transcriptome, jittered)"
    shell:
        """
        Rscript {params.script} -i {wildcards.index}_semi_reads_breakdown.csv -o {output.plot}  -n "{params.label}" > {output.model}
        """


rule plot_multi_taxonomy_boxplots:
    input:
        data=expand("{index}_HoMi_is_done_{proj}",
                index=indexes,
                proj=["synthetic_transcriptomes", "synthetic", "semi"]),
        read_lengths=os.path.join("semi", "read_lengths.csv")
    output:
        genus_plot="taxonomy_compared/combined_taxa_boxplot_genus.pdf",
        species_plot="taxonomy_compared/combined_taxa_boxplot_species.pdf"
    conda: "conda_envs/r_env.yaml"
    threads: 1
    resources:
        partition="short",
        mem_mb=int(4*1000), # MB
        runtime=int(1*60) # min
    params:
        script="Compare_all_taxonomy_results.R",
        outdir="taxonomy_compared/",
        # Not the best way to do this, but easier than alternatives
        input_files=",".join(["dna_benchmarking_synthetic_transcriptomes.f0.0.r0.0.nonhost.humann/all_bugs_list.tsv",
                             "dna_benchmarking_synthetic.f0.0.r0.0.nonhost.humann/all_bugs_list.tsv",
                             "dna_semi.f0.0.r0.0.nonhost.humann/all_bugs_list.tsv",
                             "dna_benchmarking_synthetic_transcriptomes.f0.0.r0.0.nonhost.kraken/Combined-taxonomy.tsv",
                             "dna_benchmarking_synthetic.f0.0.r0.0.nonhost.kraken/Combined-taxonomy.tsv",
                             "dna_semi.f0.0.r0.0.nonhost.kraken/Combined-taxonomy.tsv",
                             "rna_benchmarking_synthetic_transcriptomes.f0.0.r0.0.nonhost.humann/all_bugs_list.tsv",
                             "rna_benchmarking_synthetic.f0.0.r0.0.nonhost.humann/all_bugs_list.tsv",
                             "rna_semi.f0.0.r0.0.nonhost.humann/all_bugs_list.tsv",
                             "rna_benchmarking_synthetic_transcriptomes.f0.0.r0.0.nonhost.kraken/Combined-taxonomy.tsv",
                             "rna_benchmarking_synthetic.f0.0.r0.0.nonhost.kraken/Combined-taxonomy.tsv",
                             "rna_semi.f0.0.r0.0.nonhost.kraken/Combined-taxonomy.tsv"
                             ])
    shell:
        """
        mkdir -p {params.outdir}
        
        Rscript {params.script} -i {params.input_files} -l genus -r {input.read_lengths} -o {params.outdir}
        Rscript {params.script} -i {params.input_files} -l species -r {input.read_lengths} -o {params.outdir}
        
        """
