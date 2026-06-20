rule plot_expected_vs_actual_hostile_all:
    input:
        homi_dones=expand("{index}_HoMi_is_done_{proj}",
                index=indexes,
                proj=["synthetic_transcriptomes_p30",
                      "synthetic", 
                      "semi", 
                      "Pereira"]),
        metadata="Plot_reads_breakdown_metadata.csv"
    output:
        plot="all_hostile_out_benchmark.pdf"
    conda: r_conda
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(2*1000), # MB
        runtime=int(40) # min
    params:
        script="Plot_all_benchmarked_reads_breakdowns.R"
    shell:
        """
        Rscript {params.script} -m {input.metadata} -o {output.plot}
        """



rule plot_multi_taxonomy_boxplots:
    input:
        data=expand("{index}_HoMi_is_done_{proj}",
                index=indexes,
                proj=["synthetic_transcriptomes_p30",
                      "synthetic", 
                      "semi"]),
        read_lengths=os.path.join("semi", "read_lengths.csv")
    output:
        genus_plot="taxonomy_compared/combined_taxa_boxplot_genus.pdf",
        species_plot="taxonomy_compared/combined_taxa_boxplot_species.pdf"
    conda: r_conda
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(4*1000), # MB
        runtime=int(1*60) # min
    params:
        script="Compare_all_taxonomy_results.R",
        outdir="taxonomy_compared/",
        # Not the best way to do this, but easier than alternatives
        input_files=",".join(["dna_benchmarking_synthetic_transcriptomes_p30.f0.0.r0.0.nonhost.humann/all_bugs_list.tsv",
                             "dna_benchmarking_synthetic.f0.0.r0.0.nonhost.humann/all_bugs_list.tsv",
                             "dna_semi.f0.0.r0.0.nonhost.humann/all_bugs_list.tsv",
                             "dna_benchmarking_synthetic_transcriptomes_p30.f0.0.r0.0.nonhost.kraken/Combined-taxonomy.tsv",
                             "dna_benchmarking_synthetic.f0.0.r0.0.nonhost.kraken/Combined-taxonomy.tsv",
                             "dna_semi.f0.0.r0.0.nonhost.kraken/Combined-taxonomy.tsv",
                             "rna_benchmarking_synthetic_transcriptomes_p30.f0.0.r0.0.nonhost.humann/all_bugs_list.tsv",
                             "rna_benchmarking_synthetic.f0.0.r0.0.nonhost.humann/all_bugs_list.tsv",
                             "rna_semi.f0.0.r0.0.nonhost.humann/all_bugs_list.tsv",
                             "rna_benchmarking_synthetic_transcriptomes_p30.f0.0.r0.0.nonhost.kraken/Combined-taxonomy.tsv",
                             "rna_benchmarking_synthetic.f0.0.r0.0.nonhost.kraken/Combined-taxonomy.tsv",
                             "rna_semi.f0.0.r0.0.nonhost.kraken/Combined-taxonomy.tsv"
                             ])
    shell:
        """
        mkdir -p {params.outdir}
        
        Rscript {params.script} -i {params.input_files} -l genus -r {input.read_lengths} -o {params.outdir}
        Rscript {params.script} -i {params.input_files} -l species -r {input.read_lengths} -o {params.outdir}
        
        """