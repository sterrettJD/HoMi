## this subworkflow creates the hostile index and pulls required reference genomes  ##
## - i only changed partition (ability to specify in main snakefile) and added '&&' ##
## before touch in pull_reference_genomes                                           ##

rule create_alt_hostile_index:
    output:
        ref_dir=directory("t2t_hla_{hostile_index_spec}_index")
    threads: 4
    conda: hostile_conda
    resources:
        partition=hpc_partition,
        mem_mb=int(12*1000), # MB
        runtime=int(6*60) # min
    params:
        script="create_decontam_ref_human.py"
    shell:
        """
        if [[ "{wildcards.hostile_index_spec}" == "hisat2" ]]; then
            python {params.script} -m hisat2 -o {output.ref_dir}
        else
            python {params.script} -o {output.ref_dir}
        fi
        """

## added '&&' after python command so output file will only be created if the 
## python command doesnt error out
rule pull_reference_genomes:
    input:
        sample_data=metadata_file
    output:
        done="reference_genomes_downloaded"
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(4*1000), # MB
        runtime=int(2*60) # min
    params:
        script=os.path.join(synthetic_work_dir, "pull_reference_genomes.py"),
        work_dir=synthetic_work_dir,
        communities_dir=synthetic_communities_dir
    shell:
        """
        python {params.script} {input.sample_data} --work_dir {params.work_dir} && touch {output.done}
        """