#############################################################
##### Mock communities #####
rule fastq_dump_Pereira:
    output:
        fwd=os.path.join("Pereira", "{srr_id}_R1.fastq.gz"),
        rev=os.path.join("Pereira", "{srr_id}_R2.fastq.gz")
    conda: sra_conda
    threads: 1
    resources:
        partition=hpc_partition,
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
                srr_id=pereira_srr_ids),
        alt_indexes_created=expand("t2t_hla_{hostile_index_spec}_index",
                hostile_index_spec=custom_hostile_index_specs)
    output:
        "{index}_HoMi_is_done_Pereira"
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(8*1000), # MB
        runtime=int(24*60), # min
        homi_runs=1
    params:
        homi_args=homi_args
    shell:
        """
        HoMi.py {input.homi_config} {params.homi_args} --unlock && touch {output}
        """


rule plot_expected_vs_actual_mock_data:
    input:
        "{index}_HoMi_is_done_Pereira"
    output:
        plot="{index}_Pereira_benchmark.pdf",
        model="{index}_Pereira_benchmark_lm_results.txt"
    conda: r_conda
    threads: 1
    resources:
        partition=hpc_partition,
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
    conda: r_conda
    threads: 1
    resources:
        partition=hpc_partition,
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