import pandas as pd
from os.path import join as pj
from src.snake_utils import hostile_db_to_path, get_adapters_path, get_nonpareil_rmd_path, get_nonpareil_html_path

METADATA = pd.read_csv(config['METADATA'])
SAMPLES = METADATA["Sample"].tolist()
RAW_FWD_READS = METADATA[config['fwd_reads_path']]
RAW_REV_READS = METADATA[config['rev_reads_path']]

READS = ["R1", "R2"]
PROJ = config['PROJ']

HOSTILE_DB_NAME = config['hostile_db']
HOSTILE_DB_DWNLD_PATH = config['loc_for_hostile_db_download']
HOSTILE_DB_PATH = hostile_db_to_path(HOSTILE_DB_NAME, 
                                                 HOSTILE_DB_DWNLD_PATH)

trim_trunc_path = f"{config['PROJ']}.f{config['trim_fwd']}.{config['trunc_fwd']}.r{config['trim_rev']}.{config['trunc_rev']}"

rule all:
  input: 
    # Symlinked files
    expand(pj(PROJ,
              "{sample}.R1.fq.gz"),
           sample=SAMPLES),
    expand(pj(PROJ,
              "{sample}.R2.fq.gz"),
           sample=SAMPLES),
    
    # Adapter and quality-based trimming (trimmomatic)
    # Also removal of short reads
    expand(pj(f"{PROJ}.noadpt",
              "{sample}",
              "{sample}.trimmed.R1.fq"), # forward paired
           sample=SAMPLES),
    expand(pj(f"{PROJ}.noadpt",
              "{sample}",
              "{sample}.trimmed_1U"),    # forward unpaired (discard)
           sample=SAMPLES),
    expand(pj(f"{PROJ}.noadpt",
              "{sample}",
              "{sample}.trimmed.R2.fq"), # reverse unpaired
          sample=SAMPLES),
    expand(pj(f"{PROJ}.noadpt",
              "{sample}",
              "{sample}.trimmed_2U"),     # reverse unpaired (discard)
          sample=SAMPLES),

    # first pass fastQC
    expand(pj(f"{PROJ}.noadpt.fastqc",
              "{sample}.trimmed.{read}_fastqc.zip"),
           sample=SAMPLES, read=READS),

    # first pass multiQC
    pj(f"{PROJ}.noadpt.fastqc",
        "multiqc_report"),

    # Made by SeqTK
    expand(pj(trim_trunc_path,
              "{sample}.R1.fq"),
           sample=SAMPLES),
    expand(pj(trim_trunc_path,
              "{sample}.R2.fq"),
           sample=SAMPLES),

    # second pass fastQC
    expand(pj(f"{trim_trunc_path}.fastqc",
              "{sample}.{read}_fastqc.zip"),
           sample=SAMPLES, read=READS),

    # second pass multiQC
    pj(f"{trim_trunc_path}.fastqc",
        "multiqc_report"),

    # hostile index
    multiext(HOSTILE_DB_PATH,
             ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
             ".rev.1.bt2", ".rev.2.bt2"),

    # nonhost reads
    expand(pj(f"{trim_trunc_path}.nonhost",
            "{sample}.R1.fq.gz"),
          sample=SAMPLES),
    expand(pj(f"{trim_trunc_path}.nonhost",
            "{sample}.R2.fq.gz"),
          sample=SAMPLES),

    # HUMAnN pipeline
    expand(pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_pathabundance.tsv"),
            sample=SAMPLES),
    expand(pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_pathcoverage.tsv"),
            sample=SAMPLES),
    expand(pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_genefamilies.tsv"),
            sample=SAMPLES),
    expand(pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_humann_temp", 
                "{sample}_metaphlan_bugs_list.tsv"),
            sample=SAMPLES),
    
    pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_pathabundance.tsv"),
    pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_pathcoverage.tsv"),
    pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_genefamilies.tsv"),
    pj(f"{trim_trunc_path}.nonhost.humann", 
                        "all_genefamilies_grouped.tsv"),
    pj(f"{trim_trunc_path}.nonhost.humann", 
                              "all_genefamilies_grouped_named.tsv"),
    pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_bugs_list.tsv"),
    expand(pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_humann_temp", 
                "{sample}_metaphlan_bugs_list_v3.tsv"), 
                sample=SAMPLES),

    # Nonhost coverage (via nonpareil)
    expand(pj(f"{trim_trunc_path}.nonhost.nonpareil", "{sample}.npl"),
            sample=SAMPLES),
    expand(pj(f"{trim_trunc_path}.nonhost.nonpareil", "{sample}.npo"),
            sample=SAMPLES),
    expand(pj(f"{trim_trunc_path}.nonhost.nonpareil", "{sample}.npa"),
            sample=SAMPLES),
    get_nonpareil_html_path()

rule symlink_fastqs:
  output:
    FWD=pj(PROJ,"{sample}.R1.fq.gz"),
    REV=pj(PROJ,"{sample}.R2.fq.gz")
  threads: 1
  resources:
    partition="short",
    mem_mb=int(2*1000), # MB, or 2 GB
    runtime=int(1*60) # min, or 1 hours
  params:
    metadata=METADATA,
    samples=SAMPLES,
    proj=PROJ
  run:

    from os import symlink, getcwd
    from os.path import join as pj
    
    cwd = getcwd()
    
    df = params.metadata
    proj = params.proj

    sample = wildcards.sample
    
    print(sample)

    fwd = df.loc[df["Sample"]==sample, config['fwd_reads_path']].values[0]
    rev = df.loc[df["Sample"]==sample, config['rev_reads_path']].values[0]

    fwd_full = pj(cwd, fwd)
    fwd_symlink = pj(cwd, proj,f"{sample}.R1.fq.gz")
    rev_full = pj(cwd, rev)
    rev_symlink = pj(cwd, proj,f"{sample}.R2.fq.gz")
    print(f"{fwd_full} --> {fwd_symlink}")
    print(f"{rev_full} --> {rev_symlink}")
    
    # do the symlink
    symlink(fwd_full, 
            fwd_symlink)
    symlink(rev_full, 
            rev_symlink)


rule remove_adapters:
    input:
        FORWARD=pj(PROJ,
                  "{sample}.R1.fq.gz"),
        REVERSE=pj(PROJ,
                  "{sample}.R2.fq.gz")
    output:
        pj(f"{PROJ}.noadpt", "{sample}", "{sample}.trimmed.R1.fq"), # forward paired
        pj(f"{PROJ}.noadpt", "{sample}", "{sample}.trimmed_1U"),    # forward unpaired (discard)
        pj(f"{PROJ}.noadpt", "{sample}", "{sample}.trimmed.R2.fq"), # reverse unpaired
        pj(f"{PROJ}.noadpt", "{sample}", "{sample}.trimmed_2U")     # reverse unpaired (discard)

    conda: "conda_envs/trimmomatic.yaml"
    resources:
        mem_mb=int(8*1000), # 8 GB
        partition="short",
        runtime=int(12*60)
    threads: 8
    params:
      proj=PROJ,
      adpt=get_adapters_path(),
      minlen=config['min_readlen'],
      leading=config['readstart_qual_min'],
      trailing=config['readend_qual_min']
    shell:
        """
        mkdir -p {params.proj}.noadpt/{wildcards.sample}
        trimmomatic PE -threads 8 \
            -trimlog {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.trimlog \
            -summary {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.summary \
            -validatePairs {input.FORWARD} {input.REVERSE} \
            -baseout {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.trimmed \
            ILLUMINACLIP:{params.adpt}:2:30:10 SLIDINGWINDOW:4:20 LEADING:{params.leading} TRAILING:{params.trailing} MINLEN:{params.minlen} \
            -phred33
        
        mv {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.trimmed_1P {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.trimmed.R1.fq
        mv {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.trimmed_2P {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.trimmed.R2.fq

        echo "completed adapter trimming of {wildcards.sample}"
        """


rule fastQC_pass1:
  input:
    pj(f"{PROJ}.noadpt", 
        "{sample}",
        "{sample}.trimmed.{read}.fq")
  output:
    pj(f"{PROJ}.noadpt.fastqc",
        "{sample}.trimmed.{read}_fastqc.zip")

  conda: "conda_envs/fastqc.yaml"
  resources:
        partition="short",
        mem_mb=int(2*1000), # MB, or 2 GB
        runtime=int(2*60) # min, or 2 hours
  threads: 1
  params:
    proj=PROJ
  shell:
    """
    mkdir -p {params.proj}.noadpt.fastqc
    fastqc {input} -o {params.proj}.noadpt.fastqc
    """


rule multiqc_pass1:
  input:
    expand(pj(f"{PROJ}.noadpt.fastqc",
              "{sample}.trimmed.{read}_fastqc.zip"),
           sample=SAMPLES, read=READS)
  output:
    directory(pj(f"{PROJ}.noadpt.fastqc",
                "multiqc_report"))
  conda: "conda_envs/fastqc.yaml"
  resources:
        partition="short",
        mem_mb=int(2*1000), # MB, or 2 GB
        runtime=int(0.5*60) # min, or 0.5 hours
  threads: 1
  params:
    proj=PROJ
  shell:
    """
    multiqc {params.proj}.noadpt.fastqc -o {params.proj}.noadpt.fastqc/multiqc_report
    """


rule trim_forward:
  input:
    pj(f"{PROJ}.noadpt","{sample}","{sample}.trimmed.R1.fq")
  output:
    pj(trim_trunc_path,
       "{sample}.R1.fq")

  conda: "conda_envs/seqtk.yaml"
  resources:
        partition="short",
        mem_mb=int(12*1000), # MB, or 20 GB
        runtime=int(1*60) # min, or 1 hours
  threads: 1
  params:
    trim_trunc_path=trim_trunc_path,
    trim=config['trim_fwd'],
    trunc={config['trunc_fwd']}
  shell:
    """
    mkdir -p {params.trim_trunc_path}
    seqtk trimfq -b {params.trim} -e {params.trunc} {input} > {output}
    """


rule trim_reverse:
  input:
    pj(f"{PROJ}.noadpt","{sample}","{sample}.trimmed.R2.fq")
  output:
    pj(trim_trunc_path,
       "{sample}.R2.fq")

  conda: "conda_envs/seqtk.yaml"
  resources:
    partition="short",
    mem_mb=int(12*1000), # MB, or 20 GB
    runtime=int(1*60) # min, or 1 hours
  threads: 1
  params:
    trim_trunc_path=trim_trunc_path,
    trim=config['trim_rev'],
    trunc={config['trunc_rev']}
  shell:
    """
    mkdir -p {params.trim_trunc_path}
    seqtk trimfq -b {params.trim} -e {params.trunc} {input} > {output}
    """


rule fastQC_pass2:
  input:
    pj(trim_trunc_path,
       "{sample}.{read}.fq")
  output:
    pj(f"{trim_trunc_path}.fastqc",
        "{sample}.{read}_fastqc.zip")

  conda: "conda_envs/fastqc.yaml"
  resources:
    partition="short",
    mem_mb=int(2*1000), # MB, or 2 GB
    runtime=int(2*60) # min, or 0.5 hours
  threads: 1
  params:
    trim_trunc_path=trim_trunc_path
  shell:
    """
    mkdir -p {params.trim_trunc_path}.fastqc
    fastqc {input} -o {params.trim_trunc_path}.fastqc
    """


rule multiqc_pass2:
  input:
    expand(pj(f"{trim_trunc_path}.fastqc",
              "{sample}.{read}_fastqc.zip"),
           sample=SAMPLES, read=READS)
  output:
    directory(pj(f"{trim_trunc_path}.fastqc",
                  "multiqc_report"))
  conda: "conda_envs/fastqc.yaml"
  resources:
        partition="short",
        mem_mb=int(2*1000), # MB, or 2 GB
        runtime=int(0.5*60) # min, or 0.5 hours
  threads: 1
  params:
    trim_trunc_path=trim_trunc_path
  shell:
    """
    multiqc {params.trim_trunc_path}.fastqc -o {params.trim_trunc_path}.fastqc/multiqc_report
    """

rule download_hostile_db:
  output:
    multiext(HOSTILE_DB_PATH,
             ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
             ".rev.1.bt2", ".rev.2.bt2")
  resources:
    partition="short",
    mem_mb=int(4*1000), # MB, or 4 GB,
    runtime=int(8*60) # min, or 8 hours (in the case of a bad connection)
  threads: 1
  params:
    hostile_db_name=HOSTILE_DB_NAME,
    db_parent_path=HOSTILE_DB_DWNLD_PATH
  shell:
    """
    if [[ "{params.hostile_db_name}" == "human-t2t-hla-argos985" ]]; then
      mkdir -p {params.db_parent_path}
      cd {params.db_parent_path}
      wget https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla-argos985.tar
      tar -xzvf {params.hostile_db_name}.tar
      rm {params.hostile_db_name}.tar
    elif [[ "{params.hostile_db_name}" == "human-t2t-hla" ]]; then
      mkdir -p {params.db_parent_path}
      cd {params.db_parent_path}
      wget https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla.tar
      tar -xzvf {params.hostile_db_name}.tar
      rm {params.hostile_db_name}.tar
    fi
    """


rule host_filter:
  input:
    INDEX=multiext(HOSTILE_DB_PATH,
             ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
             ".rev.1.bt2", ".rev.2.bt2"),
    FWD=pj(trim_trunc_path,
       "{sample}.R1.fq"),
    REV=pj(trim_trunc_path,
       "{sample}.R2.fq")
  output:
    FWD=pj(f"{trim_trunc_path}.nonhost",
            "{sample}.R1.fq.gz"),
    REV=pj(f"{trim_trunc_path}.nonhost",
            "{sample}.R2.fq.gz")
  conda: "conda_envs/hostile.yaml"
  resources:
    partition="short",
    mem_mb=int(12*1000), # MB, or 12 GB, hostile should max at 4 (under 8 thread example), but playing it safe
    runtime=int(20*60) # min, or 20 hours
  threads: 16
  params:
    trim_trunc_path=trim_trunc_path,
    hostile_db_path=HOSTILE_DB_PATH
  shell:
    """
    hostile clean \
    --fastq1 {input.FWD} --fastq2 {input.REV} \
    --out-dir {params.trim_trunc_path}.nonhost \
    --threads {threads} \
    --index {params.hostile_db_path} \
    --debug \
    --aligner bowtie2

    # cleanup filepaths
    mv {params.trim_trunc_path}.nonhost/{wildcards.sample}.R1.clean_1.fastq.gz {output.FWD}
    mv {params.trim_trunc_path}.nonhost/{wildcards.sample}.R2.clean_2.fastq.gz {output.REV}
    """



############# RUN BIOBAKERY HUMANN PIPELINE ON NONHOST READS #############

rule setup_metaphlan:
  output:
    directory(config['metaphlan_bowtie_db'])
  resources:
    partition="short",
    mem_mb= int(32*1000), # MB
    runtime=int(60*8) # min
  threads: 8
  conda: "conda_envs/humann.yaml"
  shell:
    """
    mkdir -p {output}
    # metaphlan --install --bowtie2db  --nproc {threads} {output}
    cd {output}
    wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.tar
    tar -xvf mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.tar
    rm mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.tar
    """


rule get_biobakery_chocophlan_db:
  output:
    directory(config['chocophlan_db'])
  resources:
    partition="short",
    mem_mb= int(4*1000), # MB
    runtime=int(60*4) # min
  threads: 1
  conda: "conda_envs/humann.yaml"
  shell:
    """
    mkdir -p {output}
    humann_databases --download chocophlan full {output} --update-config yes
    """

rule get_biobakery_uniref_db:
  output:
    directory(config['uniref_db'])
  resources:
    partition="short",
    mem_mb= int(4*1000), # MB
    runtime=int(60*4) # min
  threads: 1
  conda: "conda_envs/humann.yaml"
  shell:
    """
    mkdir -p {output}
    humann_databases --download uniref uniref90_diamond {output} --update-config yes
    """


rule concat_nonhost_reads:
  input:
    FWD=pj(f"{trim_trunc_path}.nonhost",
            "{sample}.R1.fq.gz"),
    REV=pj(f"{trim_trunc_path}.nonhost",
            "{sample}.R2.fq.gz")
  output:
    pj(f"{trim_trunc_path}.nonhost.concat",
            "{sample}.fq.gz")
  resources:
    partition="short",
    mem_mb= int(8*1000), # MB
    runtime=int(60*4) # min
  threads: 1
  params:
    dirpath=f"{trim_trunc_path}.nonhost.concat"
  shell:
    """
    mkdir -p {params.dirpath}
    cat {input.FWD} {input.REV} > {output}
    """


rule run_humann_nonhost:
  input:
    METAPHLAN_DB=config['metaphlan_bowtie_db'],
    CHOCO_DB=config['chocophlan_db'],
    UNIREF_DB=config['uniref_db'],
    NONHUMAN_READS=pj(f"{trim_trunc_path}.nonhost.concat",
        "{sample}.fq.gz") 
  output:
    PATHABUND=pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_pathabundance.tsv"),
    PATHCOV=pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_pathcoverage.tsv"),
    GENEFAMS=pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_genefamilies.tsv"),
    BUGSLIST=pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_humann_temp", 
                "{sample}_metaphlan_bugs_list.tsv")
  resources:
    partition="short",
    mem_mb=int(64*1000), # MB, or 64 GB
    runtime=int(23.9*60) # min, or 23 hours
  threads: 32
  conda: "conda_envs/humann.yaml"
  params:
    dirpath=f"{trim_trunc_path}.nonhost.humann",
    metaphlan_bowtie_db=pj(config['metaphlan_bowtie_db'],"mpa_vOct22_CHOCOPhlAnSGB_202212_bt2")
  shell:
    """
    mkdir -p {params.dirpath}
    humann -i {input.NONHUMAN_READS} -o {params.dirpath}/{wildcards.sample} \
    --threads {threads} --search-mode uniref90 \
#    --metaphlan-options="--bowtie2db /scratch/Users/jost9358/Aug_23_dual_seq/data/metaphlan_db/"

    """


rule aggregate_humann_outs_nonhost:
  input:
    PATHABUND=expand(pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_pathabundance.tsv"),
            sample=SAMPLES),
    PATHCOV=expand(pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_pathcoverage.tsv"),
            sample=SAMPLES),
    GENEFAMS=expand(pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_genefamilies.tsv"),
            sample=SAMPLES),
    BUGSLIST=expand(pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_humann_temp", 
                "{sample}_metaphlan_bugs_list.tsv"),
            sample=SAMPLES)
  output:
    PATHABUND=pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_pathabundance.tsv"),
    PATHCOV=pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_pathcoverage.tsv"),
    GENEFAMS=pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_genefamilies.tsv"),
    GENEFAMS_GROUPED=pj(f"{trim_trunc_path}.nonhost.humann", 
                        "all_genefamilies_grouped.tsv"),
    GENEFAMS_GROUPED_NAMED=pj(f"{trim_trunc_path}.nonhost.humann", 
                              "all_genefamilies_grouped_named.tsv"),
    BUGSLIST=pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_bugs_list.tsv"),
    V3_NOAGG_BUGS=expand(pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_humann_temp", 
                "{sample}_metaphlan_bugs_list_v3.tsv"),
                sample=SAMPLES)

  resources:
    partition="short",
    mem_mb=int(10*1000), # MB, or 10 GB
    runtime=60 # min
  threads: 1
  conda: "conda_envs/humann.yaml"
  params:
    dirpath=f"{trim_trunc_path}.nonhost.humann"
  shell:
    """
    humann_join_tables -i {params.dirpath} -o {output.PATHABUND} --file_name pathabundance.tsv --search-subdirectories

    humann_join_tables -i {params.dirpath} -o {output.PATHCOV} --file_name pathcoverage.tsv --search-subdirectories

    humann_join_tables -i {params.dirpath} -o {output.GENEFAMS} --file_name genefamilies.tsv --search-subdirectories
    humann_regroup_table -i {output.GENEFAMS} -g uniref90_rxn -o {output.GENEFAMS_GROUPED}
    humann_rename_table -i {output.GENEFAMS_GROUPED} -n metacyc-rxn -o {output.GENEFAMS_GROUPED_NAMED}

    python utils/aggregate_metaphlan_bugslists.py -i {params.dirpath} -o {output.BUGSLIST}

    python utils/convert_mphlan_v4_to_v3.py -i {params.dirpath} 
    """


#################################
### Estimate nonhost coverage ###

rule nonpareil:
  input:
    # Only with forward reads for now. 
    # Could also run separately with reverse, but not sure there's much reason to do so.
    FWD=pj(f"{trim_trunc_path}.nonhost",
        "{sample}.R1.fq.gz")
  output:
    pj(f"{trim_trunc_path}.nonhost.nonpareil", "{sample}.npl"),
    pj(f"{trim_trunc_path}.nonhost.nonpareil", "{sample}.npo"),
    pj(f"{trim_trunc_path}.nonhost.nonpareil", "{sample}.npa")
  resources:
    partition="short",
    mem_mb=int(30*1000), # MB
    runtime=int(60*3) # min
  threads: 16
  conda: "conda_envs/nonpareil.yaml"
  params:
    dirpath=f"{trim_trunc_path}.nonhost.nonpareil"
  shell:
    """
    mkdir -p {params.dirpath}

    # Unzip fastq
    pigz -dc -p {threads} {input.FWD} > {params.dirpath}/{wildcards.sample}_temp_unzipped_input.fq

    # fastq is recommended for kmer algorithm, so defaulting to those
    nonpareil -s {params.dirpath}/{wildcards.sample}_temp_unzipped_input.fq \
    -b {params.dirpath}/{wildcards.sample} \
    -T kmer -f fastq -t {threads}

    # remove the temp file
    rm {params.dirpath}/{wildcards.sample}_temp_unzipped_input.fq
    """


rule nonpareil_curves:
  input:
    expand(pj(f"{trim_trunc_path}.nonhost.nonpareil", "{sample}.npo"),
          sample=SAMPLES)
  output:
    get_nonpareil_html_path()
  resources:
    partition="short",
    mem_mb=int(8*1000), # MB
    runtime=int(60*1) # min
  conda: "conda_envs/r_env.yaml"
  params:
    rmd_path=get_nonpareil_rmd_path(),
    output_dir=f"{trim_trunc_path}.nonhost.nonpareil",
    metadata=pj(os.getcwd(), config['METADATA'])
  shell:
    """
    Rscript \
    -e "rmarkdown::render('{params.rmd_path}', output_dir='{params.output_dir}'', params=c(npo_path='{params.output_dir}', metadata='{params.metadata}''))"
    """
