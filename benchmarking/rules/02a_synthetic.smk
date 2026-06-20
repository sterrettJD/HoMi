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
        partition=hpc_partition,
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
### Simulate transcriptomes with phreds of 30 and 40
rule simulate_synthetic_host_transcriptomes:
    input:
        sample_data=metadata_file
    output:
        data=expand(os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{{score}}/human/{sample}_unsampled_{read}.fasta"),
                    sample=samples,
                    read=reads),
        done="synthetic_transcriptomes_p{score}_created_human"
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(16*1000), # MB
        runtime=int(12*60) # min
    params:
        script=os.path.join(synthetic_work_dir, "run_polyester.R"),
        work_dir=synthetic_work_dir,
        communities_dir=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}/human"),
        polyester_phred=lambda wc: wc.get("score") 
        ##polyester_error_rate=polyester_error_rate
    shell:
        """
        ## pulling correct error rate depending on which phred is being run 
        if [[ {params.polyester_phred} == 30 ]];
        then
            errorRate=0.001
        else
            errorRate=0.0001
        fi


        Rscript {params.script} \
        -t synthetic/data/host_transcriptome.fna.gz \
        --transcriptome_url https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz \
        -g synthetic/data/host_transcriptome.gff.gz \
        --gtf_url https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.gff.gz \
        --error_rate ${{errorRate}} \
        -s {input.sample_data} \
        -n human \
        -o {params.communities_dir} && touch {output.done}
        """


rule transcriptome_fasta_to_fastq:
    input:
        data=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}/human/{sample}_unsampled_{read}.fasta")
    output:
        data=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}/human/{sample}_unsampled_{read}.fastq")
    threads: 1
    conda: bbmap_conda
    resources:
        partition=hpc_partition,
        mem_mb=int(16*1000), # MB
        runtime=int(4*60) # min
    params:
        polyester_phred=lambda wc: wc.get("score")
    shell:
        """
        reformat.sh in={input.data} out={output.data} qin=33 qout=33 qfake={params.polyester_phred}
        rm {input.data}
        """


rule subsample_fastq_to_correct_depth:
    input:
        data=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}/human/{sample}_unsampled_{read}.fastq")
    output:
        data=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}/human/{sample}_{read}.fastq")
    threads: 1
    resources:
        partition=hpc_partition,
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
        done="synthetic_microbial_transcriptomes_p{score}_created_{organism}"
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(16*1000), # MB
        runtime=int(12*60) # min
    params:
        script=os.path.join(synthetic_work_dir, "run_polyester.R"),
        work_dir=synthetic_work_dir,
        communities_dir=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}/{organism}_u"),
        polyester_phred=lambda wc: wc.get("score")
        ##polyester_error_rate=polyester_error_rate
    shell:
        """
        ## pulling correct error rate depending on which phred is being run 
        if [[ {params.polyester_phred} == 30 ]];
        then
            errorRate=0.001
        else
            errorRate=0.0001
        fi

        Rscript {params.script} \
        -t synthetic/data/{wildcards.organism}/genome/cds_from_genomic.fna \
        -g synthetic/data/{wildcards.organism}/genome/genomic.gff \
        --error_rate ${{errorRate}} \
        -s {input.sample_data} \
        -n {wildcards.organism} \
        -o {params.communities_dir} && touch {output.done}
        """


rule transcriptome_fasta_to_fastq_microbial:
    input:
        data_created="synthetic_microbial_transcriptomes_p{score}_created_{organism}",
    output:
        data=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}/{organism}_u/{sample}_unsampled_{read}.fastq")
    threads: 1
    conda: bbmap_conda
    resources:
        partition=hpc_partition,
        mem_mb=int(16*1000), # MB
        runtime=int(4*60) # min
    params:
        in_data=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}/{organism}_u/{sample}_unsampled_{read}.fasta"),
        polyester_phred=lambda wc: wc.get("score")
    shell:
        """
        reformat.sh in={params.in_data} out={output.data} qin=33 qout=33 qfake={params.polyester_phred}
        rm {params.in_data}
        """


rule subsample_fastq_to_correct_depth_microbial:
    input:
        data=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}/{organism}_u/{sample}_unsampled_{read}.fastq")
    output:
        data=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}/{organism}_s/{sample}_{read}.fastq")
    threads: 1
    resources:
        partition=hpc_partition,
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


## im not sure if this rule is going to freak out at me 
## update: it did, so I need to escape the score wildcard w double brackets
rule combine_transcriptomes:
    input:
        microbes=expand(os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{{score}}/{organism}_s/{sample}_{read}.fastq"),
                        organism=microbial_organisms, 
                        sample=samples, 
                        read=reads),
        human=expand(os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{{score}}/human/{sample}_{read}.fastq"),
                    sample=samples, 
                    read=reads)
    output:
        data=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}/combined/{sample}_{read}.fastq.gz")
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(2*1000), # MB
        runtime=int(1*60) # min
    params:
        organisms=organisms,
        synthetic_work_dir=synthetic_work_dir,
        synthetic_transcriptomes_dir="synthetic_transcriptomes_p{score}"
    run:
        import subprocess
        import os
	
        out_dir = os.path.join(params.synthetic_work_dir, f"{params.synthetic_transcriptomes_dir}")
        os.makedirs(out_dir, exist_ok=True)

        # Find the paths to each organism's transcriptome
        in_paths = [os.path.join(params.synthetic_work_dir, 
                                f"{params.synthetic_transcriptomes_dir}/{organism}_s", 
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
        partition=hpc_partition,
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
               sample=samples),
        alt_indexes_created=expand("t2t_hla_{hostile_index_spec}_index",
              hostile_index_spec=custom_hostile_index_specs)
    output:
        "{index}_HoMi_is_done_synthetic"
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(8*1000), # MB
        runtime=int(20*60), # min
        homi_runs=1
    params:
        homi_args=homi_args
    shell:
        """
        HoMi.py {input.homi_config} {params.homi_args} --unlock && touch {output}
        """

rule plot_expected_vs_actual_synthetic_communities:
    input:
        "{index}_HoMi_is_done_synthetic"
    output:
        plot="{index}_synthetic_communities_benchmark.pdf",
        model="{index}_synthetic_communities_benchmark_lm_results.txt"
    conda: r_conda
    threads: 1
    resources:
        partition=hpc_partition,
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
        homi_metadata=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}_homi_metadata.csv")
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(8*1000), # MB
        runtime=int(10) # min
    params:
        work_dir=synthetic_work_dir,
        communities_dir="synthetic_transcriptomes_p{score}/combined"
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


## need to make sure p30 homi config is renamed to include p30!
## this also might get mad at me since im only expanding one wildcard (sample) and not phred score 
rule run_HoMi_synthetic_transcriptomes:
    input:
        homi_metadata=os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{score}_homi_metadata.csv"),
        homi_config=os.path.join(synthetic_work_dir, "{index}_synthetic_transcriptomes_p{score}_HoMi_config.yaml"),
        fwd=expand(os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{{score}}/combined/{sample}_R1.fastq.gz"),
                   sample=samples),
        rev=expand(os.path.join(synthetic_work_dir, "synthetic_transcriptomes_p{{score}}/combined/{sample}_R2.fastq.gz"),
                   sample=samples),
        alt_indexes_created=expand("t2t_hla_{hostile_index_spec}_index",
                hostile_index_spec=custom_hostile_index_specs)
    output:
         "{index}_HoMi_is_done_synthetic_transcriptomes_p{score}"
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(8*1000), # MB
        runtime=int(20*60), # min
        homi_runs=1
    params:
        homi_args=homi_args
    shell:
        """
        HoMi.py {input.homi_config} {params.homi_args} --unlock && touch {output}
        """


## might need to also rename the reads_breakdown.csv file to include p30 and p40!
rule plot_expected_vs_actual_synthetic_transcriptomes:
    input:
        "{index}_HoMi_is_done_synthetic_transcriptomes_p{score}"
    output:
        plot="{index}_index_synthetic_transcriptomes_p{score}_benchmark.pdf",
        model="{index}_index_synthetic_transcriptomes_p{score}_benchmark_lm_results.txt"
    conda: r_conda
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(4*1000), # MB
        runtime=int(1*60) # min
    params:
        script="Plot_benchmarked_reads_breakdown.R",
        label="True percent host reads (GRCh38, jittered)"
    shell:
        """
        Rscript {params.script} -i {wildcards.index}_benchmarking_synthetic_transcriptomes_p{wildcards.score}_reads_breakdown.csv -o {output.plot}  -n "{params.label}" > {output.model}
        """

## use {proj} here bc it includes synthetic transcriptomes (p30/40) and synthetic communities 
rule plot_taxonomy_boxplots:
    input:
        "{index}_HoMi_is_done_{proj}"
    output:
        genus_plot="taxonomy_compared/{index}/{proj}/{method}_genus.pdf",
        species_plot="taxonomy_compared/{index}/{proj}/{method}_species.pdf"
    conda: r_conda
    threads: 1
    resources:
        partition=hpc_partition,
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