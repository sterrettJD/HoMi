#############################################################
##### Semisynthetic communities #####
rule fastq_dump_semi:
    output:
        fwd=os.path.join(semi_work_dir, "data", "raw_{srr_id}_R1.fastq.gz"),
        rev=os.path.join(semi_work_dir, "data", "raw_{srr_id}_R2.fastq.gz")
    conda: sra_conda
    threads: 1
    resources:
        partition=hpc_partition,
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
        partition=hpc_partition,
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
        fwd=expand(os.path.join(semi_work_dir, "data", "{taxon}_R1.fastq.gz"),
            taxon=semi_organisms),
        rev=expand(os.path.join(semi_work_dir, "data", "{taxon}_R2.fastq.gz"),
            taxon=semi_organisms)
    output:
        fwd=os.path.join(semi_work_dir, "samples", "{sample}_R1.fastq.gz"),
        rev=os.path.join(semi_work_dir, "samples", "{sample}_R2.fastq.gz")
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(2*1000), # MB
        runtime=int(4*60) # min
    params:
        metadata=semi_metadata_file,
        data_dir=os.path.join(semi_work_dir, "data")
    run:
        import subprocess
        import pandas as pd
        import os

        # clean up output files if they exist for safety
        for f in [output.fwd, output.rev]:
            if os.path.exists(f):
                os.remove(f)


        metadata = pd.read_csv(params.metadata, index_col="genome")
        
        # get the random seed to use. Doing it based on the hash of sample ID so that
        # the fwd and rev reads per sample get the same seed
        sample_hash = hash(wildcards.sample)

        # For each taxon, subsample it and add it to the gzipped output file
        for taxon in metadata.index:    
            depth = metadata.loc[taxon, wildcards.sample]
            print(f"sampling {taxon} to {depth} reads")

            if depth > 0:
                fwd_in = os.path.join(params.data_dir, f"{taxon}_R1.fastq.gz")
                rev_in = os.path.join(params.data_dir, f"{taxon}_R2.fastq.gz")

                for file in [fwd_in, rev_in]:
                    if not os.path.exists(file):
                        raise FileNotFoundError(f"File not found: {file}")

                tmp_fwd = f"tmp_{taxon}_{wildcards.sample}_R1.fastq.gz"
                tmp_rev = f"tmp_{taxon}_{wildcards.sample}_R2.fastq.gz"

                cmd = (
                    f"reformat.sh in1={fwd_in} in2={rev_in} "
                    f"out1={tmp_fwd} out2={tmp_rev} "
                    f"sampleseed={sample_hash} samplereadstarget={depth} "
                    f"ow=t"
                )

                print(f"running command: {cmd}")
                ran = subprocess.run(cmd, shell=True, check=True)

                print(f"Moving temp files to {output.fwd} and {output.rev}") 
                subprocess.run(f"cat {tmp_fwd} >> {output.fwd}", shell=True, check=True)
                subprocess.run(f"cat {tmp_rev} >> {output.rev}", shell=True, check=True)
                os.remove(tmp_fwd)
                os.remove(tmp_rev)




rule create_HoMi_metadata_semi:
    input:
        sample_data=semi_metadata_file
    output:
        homi_metadata=os.path.join(semi_work_dir, "semi_homi_metadata.csv")
    threads: 1
    resources:
        partition=hpc_partition,
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
                sample=semi_samples, read=reads),
        alt_indexes_created=expand("t2t_hla_{hostile_index_spec}_index",
                hostile_index_spec=custom_hostile_index_specs)
    output:
        "{index}_HoMi_is_done_semi"
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(8*1000), # MB
        runtime=int(24*60), # min
        homi_runs=1
    params:
        homi_args=semi_homi_args
    shell:
        """
        HoMi.py {input.homi_config} {params.homi_args} --unlock && touch {output}
        """


rule plot_expected_vs_actual_semi:
    input:
        "{index}_HoMi_is_done_semi"
    output:
        plot="{index}_semi_benchmark.pdf",
        model="{index}_semi_benchmark_lm_results.txt"
    conda: r_conda
    threads: 1
    resources:
        partition=hpc_partition,
        mem_mb=int(4*1000), # MB
        runtime=int(1*60) # min
    params:
        script="Plot_benchmarked_reads_breakdown.R",
        label="True percent host reads (transcriptome, jittered)"
    shell:
        """
        Rscript {params.script} -i {wildcards.index}_semi_reads_breakdown.csv -o {output.plot}  -n "{params.label}" > {output.model}
        """
