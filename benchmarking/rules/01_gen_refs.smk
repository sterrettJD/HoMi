## this subworkflow creates the hostile index and pulls required reference genomes  ##
## - i only changed partition (ability to specify in main snakefile) and added '&&' ##
## before touch in pull_reference_genomes                                           ##

rule create_alt_hostile_index:
    output:
        ref_dir=directory("t2t_hla_{hostile_index_spec}_index")
    wildcard_constraints:
        hostile_index_spec="|".join(custom_hostile_index_specs)
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


## the dna HoMi configs use the standard hostile human-t2t-hla bowtie2 index.
## HoMi only downloads it if hostile_db is "human-t2t-hla", so it is downloaded here.
rule download_hostile_dna_index:
    output:
        ref_dir=directory("t2t_hla_dna_index")
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(2*1000), # MB
        runtime=int(4*60) # min
    params:
        url="https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla.tar"
    shell:
        """
        mkdir -p {output.ref_dir}
        wget -O {output.ref_dir}/human-t2t-hla.tar.gz {params.url}
        tar -xzf {output.ref_dir}/human-t2t-hla.tar.gz -C {output.ref_dir}
        rm {output.ref_dir}/human-t2t-hla.tar.gz
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