import pandas as pd
from os.path import join as pj
from os.path import split
from src.snake_utils import hostile_db_to_path, get_adapters_path, get_nonpareil_rmd_path, get_nonpareil_html_path, get_agg_script_path, get_mphlan_conv_script_path, get_taxa_barplot_rmd_path, get_sam2bam_path, get_func_barplot_rmd_path, get_partition, get_mem, get_runtime, get_threads, get_host_mapping_samples, get_slurm_extra, get_gmm_rmd_path, get_kraken_db_loc, get_tpm_converter_path, get_host_map_method, get_rule_extra_args, get_metaphlan_index_name, get_R_installation_path, get_read_reports_path




METADATA = pd.read_csv(config["METADATA"].strip())
SAMPLES = METADATA["Sample"].tolist()
HOST_MAP_SAMPLES = get_host_mapping_samples(METADATA, sample_column="Sample")
RAW_FWD_READS = METADATA[config["fwd_reads_path"]]
RAW_REV_READS = METADATA[config["rev_reads_path"]]

READS = ["R1", "R2"]
PROJ = config["PROJ"].strip()

HOSTILE_DB_NAME = config["hostile_db"]
HOSTILE_DB_DWNLD_PATH = config["loc_for_hostile_db_download"]
HOSTILE_DB_PATH = hostile_db_to_path(HOSTILE_DB_NAME, 
                                     HOSTILE_DB_DWNLD_PATH)

trim_trunc_path = f"{config['PROJ']}.f{config['trim_fwd']}.{config['trunc_fwd']}.r{config['trim_rev']}.{config['trunc_rev']}"

rule all:
  """
  This contains most files needed for endpoints to be reached.
  """
  input: 
    # first pass multiQC
    pj(f"{PROJ}.noadpt.fastqc",
        "multiqc_report"),

    # second pass multiQC
    pj(f"{trim_trunc_path}.fastqc",
        "multiqc_report"),

    # HUMAnN pipeline
    
    pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_pathabundance.tsv"),
    pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_pathcoverage.tsv"),
    pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_genefamilies.tsv"),
    pj(f"{trim_trunc_path}.nonhost.humann", 
                        "all_genefamilies_rxn.tsv"),
    pj(f"{trim_trunc_path}.nonhost.humann", 
                              "all_genefamilies_rxn_named.tsv"),
    pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_bugs_list.tsv"),
    expand(pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_humann_temp", 
                "{sample}_metaphlan_bugs_list_v3.tsv"), 
                sample=SAMPLES),
    pj(f"{trim_trunc_path}.nonhost.humann", "Metaphlan_microshades.html"),
    pj(f"{trim_trunc_path}.nonhost.humann", "HUMAnN_microshades.html"),

    # Nonhost coverage (via nonpareil)
    pj(f"{trim_trunc_path}.nonhost.nonpareil", "nonpareil_curves.html"),

    # Host gene counts
    pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "counts.txt"),
    pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "counts.txt.summary"),
    pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "counts_tpm.tsv"),

    # Bracken combined out
    pj(f"{trim_trunc_path}.nonhost.kraken", "Combined-taxonomy.tsv"),
    pj(f"{trim_trunc_path}.nonhost.humann", "Gut_metabolic_modules.csv"),

    # reads breakdown report
    "reads_breakdown.csv"


rule symlink_fastqs:
  """
  Creates a symbolic link to the fastq files, so files can easily be accessed.
  """
  output:
    FWD=pj(PROJ,"{sample}.R1.fq.gz"),
    REV=pj(PROJ,"{sample}.R2.fq.gz")
  threads: get_threads(1, config, "symlink_fastqs")
  resources:
    partition=get_partition("short", config, "symlink_fastqs"),
    mem_mb=get_mem(int(2*1000), config, "symlink_fastqs"), # MB, or 2 GB
    runtime=get_runtime(int(1*60), config, "symlink_fastqs"), # min, or 1 hours
    slurm=get_slurm_extra(config, "symlink_fastqs")
  params:
    metadata=METADATA,
    samples=SAMPLES,
    proj=PROJ
  run:

    from os import symlink, getcwd
    from os.path import join as pj
    
    cwd = getcwd()
    
    df = params.metadata
    df["Sample"] = df["Sample"].astype(str)
    proj = params.proj

    sample = wildcards.sample
    
    print(sample)

    fwd = df.loc[df["Sample"]==sample, config["fwd_reads_path"]].values[0]
    rev = df.loc[df["Sample"]==sample, config["rev_reads_path"]].values[0]

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
  """
  Trims common sequencing adapters (using adapters found in `data/adapters.fa`) from the fastq reads.
  Also performs quality-based trimming at the start and end of reads, and reads shorter than `minlen`
  after quality-based trimming are removed entirely.
  """
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
    partition=get_partition("short", config, "remove_adapters"),
    mem_mb=get_mem(int(8*1000), config, "remove_adapters"), # 8 GB
    runtime=get_runtime(int(12*60), config, "remove_adapters"),
    slurm=get_slurm_extra(config, "remove_adapters")
  threads: get_threads(8, config, "remove_adapters")
  params:
    proj=PROJ,
    adpt=get_adapters_path(),
    minlen=config["min_readlen"],
    leading=config["readstart_qual_min"],
    trailing=config["readend_qual_min"],
    extra=get_rule_extra_args(config, "remove_adapters")
  shell:
    """
    mkdir -p {params.proj}.noadpt/{wildcards.sample}
    trimmomatic PE -threads {threads} \
        -trimlog {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.trimlog \
        -summary {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.summary \
        -validatePairs {input.FORWARD} {input.REVERSE} \
        -baseout {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.trimmed \
        ILLUMINACLIP:{params.adpt}:2:30:10 SLIDINGWINDOW:4:20 LEADING:{params.leading} TRAILING:{params.trailing} MINLEN:{params.minlen} \
        -phred33 {params.extra}
    
    mv {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.trimmed_1P {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.trimmed.R1.fq
    mv {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.trimmed_2P {params.proj}.noadpt/{wildcards.sample}/{wildcards.sample}.trimmed.R2.fq

    echo "completed adapter trimming of {wildcards.sample}"
    """


rule fastQC_pass1:
  """
  Runs FastQC to check the quality of reads in each sample.
  """
  input:
    pj(f"{PROJ}.noadpt", 
        "{sample}",
        "{sample}.trimmed.{read}.fq")
  output:
    pj(f"{PROJ}.noadpt.fastqc",
        "{sample}.trimmed.{read}_fastqc.zip")

  conda: "conda_envs/fastqc.yaml"
  resources:
        partition=get_partition("short", config, "fastQC_pass1"),
        mem_mb=get_mem(int(2*1000), config, "fastQC_pass1"), # MB, or 2 GB
        runtime=get_runtime(int(2*60), config, "fastQC_pass1"), # min, or 2 hours
        slurm=get_slurm_extra(config, "fastQC_pass1")
  threads: get_threads(1, config, "fastQC_pass1")
  params:
    proj=PROJ,
    extra=get_rule_extra_args(config, "fastQC_pass1")
  shell:
    """
    mkdir -p {params.proj}.noadpt.fastqc
    fastqc {input} -o {params.proj}.noadpt.fastqc {params.extra}
    """


rule multiqc_pass1:
  """
  Runs MultiQC to aggregate FastQC reports into one HTML file.
  """
  input:
    expand(pj(f"{PROJ}.noadpt.fastqc",
              "{sample}.trimmed.{read}_fastqc.zip"),
           sample=SAMPLES, read=READS)
  output:
    directory(pj(f"{PROJ}.noadpt.fastqc",
                "multiqc_report"))
  conda: "conda_envs/fastqc.yaml"
  resources:
        partition=get_partition("short", config, "multiqc_pass1"),
        mem_mb=get_mem(int(2*1000), config, "multiqc_pass1"), # MB, or 2 GB
        runtime=get_runtime(int(0.5*60), config, "multiqc_pass1"), # min, or 0.5 hours
        slurm=get_slurm_extra(config, "multiqc_pass1")
  threads: get_threads(1, config, "multiqc_pass1")
  params:
    proj=PROJ,
    extra=get_rule_extra_args(config, "multiqc_pass1")
  shell:
    """
    multiqc {params.proj}.noadpt.fastqc -o {params.proj}.noadpt.fastqc/multiqc_report {params.extra}
    """


rule trim_forward:
  """
  This is a secondary trimming step using SeqTK to trim/truncate 
  any set number of base pairs from the start/end of each forward read.
  """
  input:
    pj(f"{PROJ}.noadpt","{sample}","{sample}.trimmed.R1.fq")
  output:
    pj(trim_trunc_path,
       "{sample}.R1.fq")

  conda: "conda_envs/seqtk.yaml"
  resources:
        partition=get_partition("short", config, "trim_forward"),
        mem_mb=get_mem(int(12*1000), config, "trim_forward"), # MB, or 12 GB
        runtime=get_runtime(int(1*60), config, "trim_forward"), # min, or 1 hours
        slurm=get_slurm_extra(config, "trim_forward")
  threads: get_threads(1, config, "trim_forward")
  params:
    trim_trunc_path=trim_trunc_path,
    trim=config['trim_fwd'],
    trunc={config['trunc_fwd']},
    extra=get_rule_extra_args(config, "trim_forward")
  shell:
    """
    mkdir -p {params.trim_trunc_path}
    seqtk trimfq -b {params.trim} -e {params.trunc} {params.extra} {input} > {output} 
    """


rule trim_reverse:
  """
  This is a secondary trimming step using SeqTK to trim/truncate 
  any set number of base pairs from the start/end of each reverse read.
  """
  input:
    pj(f"{PROJ}.noadpt","{sample}","{sample}.trimmed.R2.fq")
  output:
    pj(trim_trunc_path,
       "{sample}.R2.fq")

  conda: "conda_envs/seqtk.yaml"
  resources:
    partition=get_partition("short", config, "trim_reverse"),
    mem_mb=get_mem(int(12*1000), config, "trim_reverse"), # MB, or 12 GB
    runtime=get_runtime(int(1*60), config, "trim_reverse"), # min, or 1 hours
    slurm=get_slurm_extra(config, "trim_reverse")
  threads: get_threads(1, config, "trim_reverse")
  params:
    trim_trunc_path=trim_trunc_path,
    trim=config['trim_rev'],
    trunc={config['trunc_rev']},
    extra=get_rule_extra_args(config, "trim_reverse")
  shell:
    """
    mkdir -p {params.trim_trunc_path}
    seqtk trimfq -b {params.trim} -e {params.trunc} {params.extra} {input} > {output}
    """


rule fastQC_pass2:
  """
  Runs FastQC to check the quality of reads in each sample.
  """
  input:
    pj(trim_trunc_path,
       "{sample}.{read}.fq")
  output:
    pj(f"{trim_trunc_path}.fastqc",
        "{sample}.{read}_fastqc.zip")

  conda: "conda_envs/fastqc.yaml"
  resources:
    partition=get_partition("short", config, "fastQC_pass2"),
    mem_mb=get_mem(int(2*1000), config, "fastQC_pass2"), # MB, or 2 GB
    runtime=get_runtime(int(2*60), config, "fastQC_pass2"), # min, or 0.5 hours
    slurm=get_slurm_extra(config, "fastQC_pass2")
  threads: get_threads(1, config, "fastQC_pass2")
  params:
    trim_trunc_path=trim_trunc_path,
    extra=get_rule_extra_args(config, "fastQC_pass2")
  shell:
    """
    mkdir -p {params.trim_trunc_path}.fastqc
    fastqc {input} -o {params.trim_trunc_path}.fastqc {params.extra}
    """


rule multiqc_pass2:
  """
  Runs MultiQC to aggregate FastQC reports into one HTML file.
  """
  input:
    expand(pj(f"{trim_trunc_path}.fastqc",
              "{sample}.{read}_fastqc.zip"),
           sample=SAMPLES, read=READS)
  output:
    directory(pj(f"{trim_trunc_path}.fastqc",
                  "multiqc_report"))
  conda: "conda_envs/fastqc.yaml"
  resources:
    partition=get_partition("short", config, "multiqc_pass2"),
    mem_mb=get_mem(int(2*1000), config, "multiqc_pass2"), # MB, or 2 GB
    runtime=get_runtime(int(0.5*60), config, "multiqc_pass2"), # min, or 0.5 hours
    slurm=get_slurm_extra(config, "multiqc_pass2")
  threads: get_threads(1, config, "multiqc_pass2")
  params:
    trim_trunc_path=trim_trunc_path,
    extra=get_rule_extra_args(config, "multiqc_pass2")
  shell:
    """
    multiqc {params.trim_trunc_path}.fastqc -o {params.trim_trunc_path}.fastqc/multiqc_report {params.extra}
    """


rule download_hostile_db:
  """
  This rule downloads the host index to be used for Hostile's host read filtering.
  """
  output:
    multiext(HOSTILE_DB_PATH,
             ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
             ".rev.1.bt2", ".rev.2.bt2")
  resources:
    partition=get_partition("short", config, "download_hostile_db"),
    mem_mb=get_mem(int(4*1000), config, "download_hostile_db"), # MB, or 4 GB,
    runtime=get_runtime(int(8*60), config, "download_hostile_db"), # min, or 8 hours (in the case of a bad connection)
    slurm=get_slurm_extra(config, "download_hostile_db")
  threads: get_threads(1, config, "download_hostile_db")
  params:
    hostile_db_name=HOSTILE_DB_NAME,
    db_parent_path=HOSTILE_DB_DWNLD_PATH
  shell:
    """
    if [[ "{params.hostile_db_name}" == "human-t2t-hla-argos985" ]]; then
      mkdir -p {params.db_parent_path}
      cd {params.db_parent_path}
      wget https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla-argos985.tar
      tar -xvf {params.hostile_db_name}.tar
      rm {params.hostile_db_name}.tar
    elif [[ "{params.hostile_db_name}" == "human-t2t-hla" ]]; then
      mkdir -p {params.db_parent_path}
      cd {params.db_parent_path}
      wget https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla.tar
      tar -xvf {params.hostile_db_name}.tar
      rm {params.hostile_db_name}.tar
    fi
    """


rule host_filter:
  """
  This removes host reads using the tool Hostile.
  """
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
            "{sample}.R2.fq.gz"),
    REPORT=pj(f"{trim_trunc_path}.nonhost",
            "{sample}.report")
  conda: "conda_envs/hostile.yaml"
  resources:
    partition=get_partition("short", config, "host_filter"),
    mem_mb=get_mem(int(12*1000), config, "host_filter"), # MB, or 12 GB, hostile should max at 4 (under 8 thread example), but playing it safe
    runtime=get_runtime(int(20*60), config, "host_filter"), # min, or 20 hours
    slurm=get_slurm_extra(config, "host_filter")
  threads: get_threads(16, config, "host_filter")
  params:
    trim_trunc_path=trim_trunc_path,
    hostile_db_path=HOSTILE_DB_PATH,
    extra=get_rule_extra_args(config, "host_filter")
  shell:
    """
    hostile clean \
    --fastq1 {input.FWD} --fastq2 {input.REV} \
    --out-dir {params.trim_trunc_path}.nonhost \
    --threads {threads} \
    --index {params.hostile_db_path} \
    --debug \
    --aligner bowtie2 \
    {params.extra} > {output.REPORT}

    # cleanup filepaths
    mv {params.trim_trunc_path}.nonhost/{wildcards.sample}.R1.clean_1.fastq.gz {output.FWD}
    mv {params.trim_trunc_path}.nonhost/{wildcards.sample}.R2.clean_2.fastq.gz {output.REV}
    """


############# INSTALL R PACKAGES NEEDED FOR NONHOST ANALYSIS #############

rule install_R_packages:
  """
  This installs R packages needed 
  """
  output:
    "R_packages_installed"
  conda: "conda_envs/r_env.yaml"
  resources:
    partition=get_partition("short", config, "install_R_packages"),
    mem_mb=get_mem(int(12*1000), config, "install_R_packages"), # MB, or 12 GB
    runtime=get_runtime(int(4*60), config, "install_R_packages"), # min, or 4 hours
    slurm=get_slurm_extra(config, "install_R_packages")
  threads: get_threads(1, config, "install_R_packages")
  params:
    installation_script=get_R_installation_path()
  shell:
    """
    Rscript {params.installation_script}
    """

############# RUN BIOBAKERY HUMANN PIPELINE ON NONHOST READS #############

rule setup_metaphlan:
  """
  This rule installs the metaphlan database. If the user specifies "latest" or nothing as metaphlan_index_name in the config, 
  it will just download the default latest metaphlan database. Otherwise, users can specify specific versions. 
  The database will be downloaded in the directory specified by metaphlan_bowtie_db in the config file.
  """
  output:
    loc=directory(config['metaphlan_bowtie_db'])
  resources:
    partition=get_partition("short", config, "setup_metaphlan"),
    mem_mb=get_mem(int(32*1000), config, "setup_metaphlan"), # MB
    runtime=get_runtime(int(8*60), config, "setup_metaphlan"), # min
    slurm=get_slurm_extra(config, "setup_metaphlan")
  threads: get_threads(8, config, "setup_metaphlan")
  conda: "conda_envs/humann.yaml"
  params:
    index_name=get_metaphlan_index_name(config),
    extra=get_rule_extra_args(config, "setup_metaphlan")
  shell:
    """
    mkdir -p {output.loc}

    if [ "{params.index_name}" = "latest" ]; then
      metaphlan --install --nproc {threads} --bowtie2db {output.loc} {params.extra}
    
    else
      metaphlan --install --nproc {threads} --bowtie2db {output.loc} --index {params.index_name} {params.extra}
    
    fi
    
    # Option to do it manually if --install doesn't seem to work
    # cd {output}
    # Can specify whatever version you want here
    # wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.tar
    # tar -xvf mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.tar
    # rm mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.tar
    """


rule get_biobakery_chocophlan_db:
  """
  This rule installs the ChocoPhlAn database used for HUMAnN.
  """
  output:
    directory(config['chocophlan_db'])
  resources:
    partition=get_partition("short", config, "get_biobakery_chocophlan_db"),
    mem_mb=get_mem(int(4*1000), config, "get_biobakery_chocophlan_db"), # MB
    runtime=get_runtime(int(4*60), config, "get_biobakery_chocophlan_db"), # min
    slurm=get_slurm_extra(config, "get_biobakery_chocophlan_db")
  threads: get_threads(1, config, "get_biobakery_chocophlan_db")
  conda: "conda_envs/humann.yaml"
  params:
    extra=get_rule_extra_args(config, "get_biobakery_chocophlan_db")
  shell:
    """
    mkdir -p {output}
    humann_databases --download chocophlan full {output} --update-config yes {params.extra}
    """


rule get_biobakery_uniref_db:
  """
  This rule installs the UniRef database used for HUMAnN.
  """
  output:
    directory(config['uniref_db'])
  resources:
    partition=get_partition("short", config, "get_biobakery_uniref_db"),
    mem_mb=get_mem(int(4*1000), config, "get_biobakery_uniref_db"), # MB
    runtime=get_runtime(int(4*60), config, "get_biobakery_uniref_db"), # min
    slurm=get_slurm_extra(config, "get_biobakery_uniref_db")
  threads: get_threads(1, config, "get_biobakery_uniref_db")
  conda: "conda_envs/humann.yaml"
  params:
    extra=get_rule_extra_args(config, "get_biobakery_uniref_db")
  shell:
    """
    mkdir -p {output}
    humann_databases --download uniref uniref90_diamond {output} --update-config yes {params.extra}
    """


rule get_utility_mapping_db:
  """
  This rule installs the utility mapping database used for HUMAnN.
  This database is used to group/rename UniRef protein families to other
  versions, such as Enzyme Commission numbers, KEGG Orthologies, etc.
  """
  output:
    directory(config['utility_mapping_db'])
  resources:
    partition=get_partition("short", config, "get_utility_mapping_db"),
    mem_mb=get_mem(int(4*1000), config, "get_utility_mapping_db"), # MB
    runtime=get_runtime(int(1*60), config, "get_utility_mapping_db"), # min
    slurm=get_slurm_extra(config, "get_utility_mapping_db")
  threads: get_threads(1, config, "get_utility_mapping_db")
  conda: "conda_envs/humann.yaml"
  params:
    extra=get_rule_extra_args(config, "get_utility_mapping_db")
  shell:
    """
    mkdir -p {output}
    humann_databases --download utility_mapping full {output} --update-config yes {params.extra}
    """


rule concat_nonhost_reads:
  """
  This rule installs concatenates forward and reverse reads into the same file for HUMAnN.
  """
  input:
    FWD=pj(f"{trim_trunc_path}.nonhost",
            "{sample}.R1.fq.gz"),
    REV=pj(f"{trim_trunc_path}.nonhost",
            "{sample}.R2.fq.gz")
  output:
    pj(f"{trim_trunc_path}.nonhost.concat",
            "{sample}.fq.gz")
  resources:
    partition=get_partition("short", config, "concat_nonhost_reads"),
    mem_mb=get_mem(int(8*1000), config, "concat_nonhost_reads"), # MB
    runtime=get_runtime(int(4*60), config, "concat_nonhost_reads"), # min
    slurm=get_slurm_extra(config, "concat_nonhost_reads")
  threads: get_threads(1, config, "concat_nonhost_reads")
  params:
    dirpath=f"{trim_trunc_path}.nonhost.concat"
  shell:
    """
    mkdir -p {params.dirpath}
    cat {input.FWD} {input.REV} > {output}
    """


rule run_humann_nonhost:
  """
  This rule runs the HUMAnN pipeline, including MetaPhlAn.
  """
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
    partition=get_partition("short", config, "run_humann_nonhost"),
    mem_mb=get_mem(int(64*1000), config, "run_humann_nonhost"), # MB, or 64 GB
    runtime=get_runtime(int(23.9*60), config, "run_humann_nonhost"), # min, or 23 hours
    slurm=get_slurm_extra(config, "run_humann_nonhost")
  threads: get_threads(32, config, "run_humann_nonhost")
  conda: "conda_envs/humann.yaml"
  params:
    dirpath=f"{trim_trunc_path}.nonhost.humann",
    metaphlan_bowtie_db=config["metaphlan_bowtie_db"],
    metaphlan_index=get_metaphlan_index_name(config),
    extra=get_rule_extra_args(config, "run_humann_nonhost")
  shell:
    """
    mkdir -p {params.dirpath}

    if [ "{params.metaphlan_index}" = "latest" ]; then
      # read index name form the latest file
      metaphlan_db_name="$(<{params.metaphlan_bowtie_db}/mpa_latest)"
      
      humann -i {input.NONHUMAN_READS} -o {params.dirpath}/{wildcards.sample} \
      --threads {threads} --search-mode uniref90 \
      --metaphlan-options "--bowtie2db {params.metaphlan_bowtie_db}" \
      {params.extra}

    else
      humann -i {input.NONHUMAN_READS} -o {params.dirpath}/{wildcards.sample} \
      --threads {threads} --search-mode uniref90 \
      --metaphlan-options "--bowtie2db {params.metaphlan_bowtie_db}" \
      {params.extra}

    fi
    """


rule aggregate_humann_outs_nonhost:
  """
  This rule aggregates the HUMAnN and MetaPhlAn outputs across all samples into tsv files with all samples.
  It also groups protein families and renames them based on Metacyc RXNs, KOs, KEGG Modules, KEGG Pathways, 
  GO terms, EGGNOG, and PFAMs.
  """
  input:
    MAPPING_DB=config['utility_mapping_db'],
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
    GENEFAMS_GROUPED_RXN=pj(f"{trim_trunc_path}.nonhost.humann", 
                        "all_genefamilies_rxn.tsv"),
    GENEFAMS_GROUPED_NAMED_RXN=pj(f"{trim_trunc_path}.nonhost.humann", 
                              "all_genefamilies_rxn_named.tsv"),
    GENEFAMS_GROUPED_KO=pj(f"{trim_trunc_path}.nonhost.humann", 
                        "all_genefamilies_ko.tsv"),
    GENEFAMS_GROUPED_NAMED_KO=pj(f"{trim_trunc_path}.nonhost.humann", 
                              "all_genefamilies_ko_named.tsv"),
    GENEFAMS_GROUPED_NAMED_KP=pj(f"{trim_trunc_path}.nonhost.humann", 
                              "all_genefamilies_kp_named.tsv"),
    GENEFAMS_GROUPED_NAMED_KM=pj(f"{trim_trunc_path}.nonhost.humann", 
                              "all_genefamilies_km_named.tsv"),
    GENEFAMS_GROUPED_EN=pj(f"{trim_trunc_path}.nonhost.humann", 
                        "all_genefamilies_eggnog.tsv"),
    GENEFAMS_GROUPED_GO=pj(f"{trim_trunc_path}.nonhost.humann", 
                        "all_genefamilies_go.tsv"),
    GENEFAMS_GROUPED_NAMED_GO=pj(f"{trim_trunc_path}.nonhost.humann", 
                              "all_genefamilies_go_named.tsv"),
    GENEFAMS_GROUPED_PFAM=pj(f"{trim_trunc_path}.nonhost.humann", 
                        "all_genefamilies_pfam.tsv"),
    GENEFAMS_GROUPED_NAMED_PFAM=pj(f"{trim_trunc_path}.nonhost.humann", 
                              "all_genefamilies_pfam_named.tsv"),
    BUGSLIST=pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_bugs_list.tsv"),
    V3_NOAGG_BUGS=expand(pj(f"{trim_trunc_path}.nonhost.humann",
                "{sample}", "{sample}_humann_temp", 
                "{sample}_metaphlan_bugs_list_v3.tsv"),
                sample=SAMPLES)

  resources:
    partition=get_partition("short", config, "aggregate_humann_outs_nonhost"),
    mem_mb=get_mem(int(10*1000), config, "aggregate_humann_outs_nonhost"), # MB, or 10 GB
    runtime=get_runtime(int(1*60), config, "aggregate_humann_outs_nonhost"), # min
    slurm=get_slurm_extra(config, "aggregate_humann_outs_nonhost")
  threads: get_threads(1, config, "aggregate_humann_outs_nonhost")
  conda: "conda_envs/humann.yaml"
  params:
    dirpath=f"{trim_trunc_path}.nonhost.humann",
    agg_bugslists=get_agg_script_path(),
    mphlan_conv=get_mphlan_conv_script_path()
  shell:
    """
    humann_join_tables -i {params.dirpath} -o {output.PATHABUND} --file_name pathabundance.tsv --search-subdirectories

    humann_join_tables -i {params.dirpath} -o {output.PATHCOV} --file_name pathcoverage.tsv --search-subdirectories

    humann_join_tables -i {params.dirpath} -o {output.GENEFAMS} --file_name genefamilies.tsv --search-subdirectories
    
    # Metacyc RXN renaming
    humann_regroup_table -i {output.GENEFAMS} -g uniref90_rxn -o {output.GENEFAMS_GROUPED_RXN}
    humann_rename_table -i {output.GENEFAMS_GROUPED_RXN} -n metacyc-rxn -o {output.GENEFAMS_GROUPED_NAMED_RXN}

    # KO renaming
    humann_regroup_table -i {output.GENEFAMS} -g uniref90_ko -o {output.GENEFAMS_GROUPED_KO}
    humann_rename_table -i {output.GENEFAMS_GROUPED_KO} -n kegg-orthology -o {output.GENEFAMS_GROUPED_NAMED_KO}
    humann_rename_table -i {output.GENEFAMS_GROUPED_KO} -n kegg-pathway -o {output.GENEFAMS_GROUPED_NAMED_KP}
    humann_rename_table -i {output.GENEFAMS_GROUPED_KO} -n kegg-module -o {output.GENEFAMS_GROUPED_NAMED_KM}

    # EGGNOG renaming
    humann_regroup_table -i {output.GENEFAMS} -g uniref90_eggnog -o {output.GENEFAMS_GROUPED_EN}

    # GO renaming
    humann_regroup_table -i {output.GENEFAMS} -g uniref90_go -o {output.GENEFAMS_GROUPED_GO}
    humann_rename_table -i {output.GENEFAMS_GROUPED_GO} -n go -o {output.GENEFAMS_GROUPED_NAMED_GO}

    # PFAM renaming
    humann_regroup_table -i {output.GENEFAMS} -g uniref90_pfam -o {output.GENEFAMS_GROUPED_PFAM}
    humann_rename_table -i {output.GENEFAMS_GROUPED_PFAM} -n pfam -o {output.GENEFAMS_GROUPED_NAMED_PFAM}

    python {params.agg_bugslists} -i {params.dirpath} -o {output.BUGSLIST}

    python {params.mphlan_conv} -i {params.dirpath} 
    """


rule taxa_barplot:
  """
  This rule runs a .Rmd file that creates microshaded taxa barplots using MetaPhlan taxonomic profiles. 
  """
  input:
    pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_bugs_list.tsv"),
    "R_packages_installed"
  output:
    pj(f"{trim_trunc_path}.nonhost.humann", 
                "Metaphlan_microshades.html")
  resources:
    partition=get_partition("short", config, "taxa_barplot"),
    mem_mb=get_mem(int(10*1000), config, "taxa_barplot"), # MB, or 10 GB
    runtime=get_runtime(int(2*60), config, "taxa_barplot"), # min
    slurm=get_slurm_extra(config, "taxa_barplot")
  threads: get_threads(1, config, "taxa_barplot")
  conda: "conda_envs/r_env.yaml"
  params:
    rmd_path=get_taxa_barplot_rmd_path(),
    bugslist=pj(os.getcwd(), f"{trim_trunc_path}.nonhost.humann",
                "all_bugs_list.tsv"),
    output_dir=pj(os.getcwd(), f"{trim_trunc_path}.nonhost.humann"),
    metadata=pj(os.getcwd(), config['METADATA'])
  shell:
    """
    Rscript \
    -e "rmarkdown::render('{params.rmd_path}', output_dir='{params.output_dir}', params=list(bugslist='{params.bugslist}', metadata='{params.metadata}', directory='{params.output_dir}'))"
    """


rule func_barplot_rxn:
  """
  This rule runs a .Rmd file that creates microshaded functional profile barplots 
  using HUMAnN functional profiles. Specifically, it uses Enzyme Commission numbers from 
  the Metacyc RXN number grouped protein families to show the hierarchical breakdown of
  enzyme classes in your sample.
  """
  input:
    GENEFAMS_RXN=pj(f"{trim_trunc_path}.nonhost.humann", 
                    "all_genefamilies_rxn_named.tsv"),
    INSTALLED="R_packages_installed"
  output:
    pj(f"{trim_trunc_path}.nonhost.humann", 
                "HUMAnN_microshades.html")
  resources:
    partition=get_partition("short", config, "func_barplot_rxn"),
    mem_mb=get_mem(int(10*1000), config, "func_barplot_rxn"), # MB, or 10 GB
    runtime=get_runtime(int(2*60), config, "func_barplot_rxn"), # min
    slurm=get_slurm_extra(config, "func_barplot_rxn")
  threads: get_threads(1, config, "func_barplot_rxn")
  conda: "conda_envs/r_env.yaml"
  params:
    rmd_path=get_func_barplot_rmd_path(),
    gene_table_rxn=pj(os.getcwd(), f"{trim_trunc_path}.nonhost.humann", 
                    "all_genefamilies_rxn_named.tsv"),
    output_dir=pj(os.getcwd(), f"{trim_trunc_path}.nonhost.humann"),
    metadata=pj(os.getcwd(), config['METADATA'])
  shell:
    """
    Rscript \
    -e "rmarkdown::render('{params.rmd_path}', output_dir='{params.output_dir}', params=list(genetable='{params.gene_table_rxn}', metadata='{params.metadata}', directory='{params.output_dir}'))"
    """


rule calc_gut_metabolic_modules:
  """
  This rule runs a .Rmd file that groups the KOs from HUMAnN into gut metabolic modules 
  with biologically relevant functions. Abundances are summed within modules and then
  written to a .csv file. 
  """
  input:
    GENEFAMS_KO=pj(f"{trim_trunc_path}.nonhost.humann", 
                    "all_genefamilies_ko_named.tsv"),
    INSTALLED="R_packages_installed"
  output:
    HTML=pj(f"{trim_trunc_path}.nonhost.humann", 
                "Gut_metabolic_modules.html"),
    GMM_TABLE=pj(f"{trim_trunc_path}.nonhost.humann", 
                "Gut_metabolic_modules.csv")
  resources:
    partition=get_partition("short", config, "calc_gut_metabolic_modules"),
    mem_mb=get_mem(int(8*1000), config, "calc_gut_metabolic_modules"), # MB, or 8 GB
    runtime=get_runtime(int(2*60), config, "calc_gut_metabolic_modules"), # min
    slurm=get_slurm_extra(config, "calc_gut_metabolic_modules")
  threads: get_threads(1, config, "calc_gut_metabolic_modules")
  conda: "conda_envs/r_env.yaml"
  params:
    rmd_path=get_gmm_rmd_path(),
    gene_table_ko=pj(os.getcwd(), f"{trim_trunc_path}.nonhost.humann", 
                    "all_genefamilies_ko_named.tsv"),
    gmm_output=pj(os.getcwd(), f"{trim_trunc_path}.nonhost.humann", 
                    "Gut_metabolic_modules.csv"),
    output_dir=pj(os.getcwd(), f"{trim_trunc_path}.nonhost.humann")
  shell:
    """
    Rscript \
    -e "rmarkdown::render('{params.rmd_path}', output_dir='{params.output_dir}', params=list(input_file='{params.gene_table_ko}', output_file='{params.gmm_output}'))"
    """



#####################################
### Kraken + Bracken for taxonomy ###

kraken_db_loc = get_kraken_db_loc(default=pj("data", "kraken2_db"), config=config)

rule get_kraken_db:
  """
  This rule downloads and builds the Kraken2 database for taxonomic classification.
  """
  output: 
    HASH=pj(kraken_db_loc, "hash.k2d"),
    OPTS=pj(kraken_db_loc, "opts.k2d"),
    SEQ2ID=pj(kraken_db_loc, "seqid2taxid.map"),
    TAXO=pj(kraken_db_loc, "taxo.k2d"),
    LIB=directory(pj(kraken_db_loc, "library")),
    TAX=directory(pj(kraken_db_loc, "taxonomy"))
  resources:
    partition=get_partition("short", config, "get_kraken_db"),
    mem_mb=get_mem(int(188*1000), config, "get_kraken_db"), # MB
    runtime=get_runtime(int(23.9*60), config, "get_kraken_db"), # min 
    slurm=get_slurm_extra(config, "get_kraken_db")
  threads: get_threads(64, config, "get_kraken_db")
  conda: "conda_envs/kraken.yaml"
  params:
    database_dir=kraken_db_loc,
    extra=get_rule_extra_args(config, "get_kraken_db")
  shell:
    """
    kraken2-build --standard --db {params.database_dir} --threads {threads} {params.extra}
    """


rule run_kraken:
  """
  This rule runs Kraken2 for each sample.
  """
  input:
    FWD=pj(f"{trim_trunc_path}.nonhost",
            "{sample}.R1.fq.gz"),
    REV=pj(f"{trim_trunc_path}.nonhost",
            "{sample}.R2.fq.gz"),
    HASH=pj(kraken_db_loc, "hash.k2d") # if the full db dir is passed here, kraken will be rerun after bracken building
                                              # because bracken-build modifies that directory
  output: 
    OUTFILE=pj(f"{trim_trunc_path}.nonhost.kraken", 
              "{sample}.kraken.txt"),
    REPORT=pj(f"{trim_trunc_path}.nonhost.kraken", 
              "{sample}.kreport2")
  resources:
    partition=get_partition("short", config, "run_kraken"),
    mem_mb=get_mem(int(84*1000), config, "run_kraken"), # MB 
    runtime=get_runtime(int(2*60), config, "run_kraken"), # min 
    slurm=get_slurm_extra(config, "run_kraken")
  threads: get_threads(16, config, "run_kraken")
  conda: "conda_envs/kraken.yaml"
  params:
    out_dir=f"{trim_trunc_path}.nonhost.kraken",
    database=kraken_db_loc,
    extra=get_rule_extra_args(config, "run_kraken")
  shell:
    """
    mkdir -p {params.out_dir}

    # This has to all be on one line due to the way kraken2 parses it...
    kraken2 --gzip-compressed --paired --db {params.database} --threads {threads} --output {output.OUTFILE} --report {output.REPORT} --classified-out {params.out_dir}/{wildcards.sample}_classified#.fq --unclassified-out {params.out_dir}/{wildcards.sample}_unclassified#.fq {params.extra} {input.FWD} {input.REV}

    """


rule build_bracken:
  """
  This rule downloads and builds the Bracken database for taxonomic abundance estimation.
  """
  input:
    pj(kraken_db_loc, "hash.k2d")
  output:
    pj(kraken_db_loc, "database150mers.kraken"),
    pj(kraken_db_loc, "database150mers.kmer_distrib")
  resources:
    partition=get_partition("short", config, "build_bracken"),
    mem_mb=get_mem(int(128*1000), config, "build_bracken"), # MB
    runtime=get_runtime(int(4*60), config, "build_bracken"), # min
    slurm=get_slurm_extra(config, "build_bracken")
  threads: get_threads(32, config, "build_bracken")
  conda: "conda_envs/kraken.yaml"
  params:
    database=kraken_db_loc,
    extra=get_rule_extra_args(config, "build_bracken")
  shell:
    """
    bracken-build -d {params.database} -t {threads} -l 150 {params.extra}
    """


rule run_bracken:
  """
  This rule runs Bracken to estimate taxonomic abundance for each sample.
  """
  input:
    KRAKEN_HASH=pj(kraken_db_loc, "hash.k2d"),
    LMERS=pj(kraken_db_loc, "database150mers.kraken"),
    LMERS_DIST=pj(kraken_db_loc, "database150mers.kmer_distrib"),
    REPORT=pj(f"{trim_trunc_path}.nonhost.kraken", 
              "{sample}.kreport2")
  output:
    REPORT=pj(f"{trim_trunc_path}.nonhost.kraken", 
              "{sample}.bracken")
  resources:
    partition=get_partition("short", config, "run_bracken"),
    mem_mb=get_mem(int(8*1000), config, "run_bracken"), # MB
    runtime=get_runtime(int(2*60), config, "run_bracken"), # min
    slurm=get_slurm_extra(config, "run_bracken")
  threads: get_threads(1, config, "run_bracken")
  conda: "conda_envs/kraken.yaml"
  params:
    database=kraken_db_loc,
    extra=get_rule_extra_args(config, "run_bracken")
  shell:
    """
    bracken -d {params.database} -i {input.REPORT} -o {output.REPORT} -r 150 -l S -t 10 {params.extra}
    """


rule aggregate_bracken:
  """
  This rule aggregates Bracken outputs across samples.
  """
  input:
    expand(pj(f"{trim_trunc_path}.nonhost.kraken", 
              "{sample}.bracken"),
           sample=SAMPLES)
  output:
    pj(f"{trim_trunc_path}.nonhost.kraken", "Combined-taxonomy.tsv")
  resources:
    partition=get_partition("short", config, "aggregate_bracken"),
    mem_mb=get_mem(int(4*1000), config, "aggregate_bracken"), # MB
    runtime=get_runtime(int(4*60), config, "aggregate_bracken"), # min
    slurm=get_slurm_extra(config, "aggregate_bracken")
  threads: get_threads(1, config, "aggregate_bracken")
  conda: "conda_envs/kraken.yaml"
  params:
    dirpath=f"{trim_trunc_path}.nonhost.kraken"
  shell:
    """
    cd {params.dirpath}
    
    # TODO: should probably put this in its own rule, but this works for now
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xzvf taxdump.tar.gz

    ( ls *.bracken ) > bracken-sample-name-map.tsv

    bit-combine-bracken-and-add-lineage -i bracken-sample-name-map.tsv -o Combined-taxonomy.tsv -d .

    rm *.dmp
    rm readme.txt
    rm taxdump.tar.gz
    """


#################################
### Estimate nonhost coverage ###

rule nonpareil:
  """
  This rule estimates coverage of the nonhost reads based on kmer redundancy.
  """
  input:
    # Only with forward reads for now. 
    # Could also run separately with reverse, but not sure there's much reason to do so.
    FWD=pj(f"{trim_trunc_path}.nonhost",
        "{sample}.R1.fq.gz"),
    INSTALLED="R_packages_installed"
  output:
    pj(f"{trim_trunc_path}.nonhost.nonpareil", "{sample}.npl"),
    pj(f"{trim_trunc_path}.nonhost.nonpareil", "{sample}.npo"),
    pj(f"{trim_trunc_path}.nonhost.nonpareil", "{sample}.npa")
  resources:
    partition=get_partition("short", config, "nonpareil"),
    mem_mb=get_mem(int(30*1000), config, "nonpareil"), # MB
    runtime=get_runtime(int(3*60), config, "nonpareil"), # min
    slurm=get_slurm_extra(config, "nonpareil")
  threads: get_threads(16, config, "nonpareil")
  conda: "conda_envs/nonpareil.yaml"
  params:
    dirpath=f"{trim_trunc_path}.nonhost.nonpareil",
    extra=get_rule_extra_args(config, "nonpareil")
  shell:
    """
    mkdir -p {params.dirpath}

    # Unzip fastq
    pigz -dc -p {threads} {input.FWD} > {params.dirpath}/{wildcards.sample}_temp_unzipped_input.fq

    # fastq is recommended for kmer algorithm, so defaulting to those
    nonpareil -s {params.dirpath}/{wildcards.sample}_temp_unzipped_input.fq \
    -b {params.dirpath}/{wildcards.sample} \
    -T kmer -f fastq -t {threads} \
    {params.extra}


    # remove the temp file
    rm {params.dirpath}/{wildcards.sample}_temp_unzipped_input.fq
    """


rule nonpareil_curves:
  """
  This rule runs a .Rmd file to visualize coverage of the nonhost fraction across samples.
  """
  input:
    expand(pj(f"{trim_trunc_path}.nonhost.nonpareil", "{sample}.npo"),
          sample=SAMPLES)
  output:
    pj(f"{trim_trunc_path}.nonhost.nonpareil", "nonpareil_curves.html")
  resources:
    partition=get_partition("short", config, "nonpareil_curves"),
    mem_mb=get_mem(int(8*1000), config, "nonpareil_curves"), # MB
    runtime=get_runtime(int(1*60), config, "nonpareil_curves"), # min
    slurm=get_slurm_extra(config, "nonpareil_curves")
  conda: "conda_envs/r_env.yaml"
  threads: get_threads(1, config, "nonpareil_curves")
  params:
    rmd_path=get_nonpareil_rmd_path(),
    output_dir=pj(os.getcwd(), f"{trim_trunc_path}.nonhost.nonpareil"),
    metadata=pj(os.getcwd(), config['METADATA'])
  shell:
    """
    Rscript \
    -e "rmarkdown::render('{params.rmd_path}', output_dir='{params.output_dir}', params=list(npo_path='{params.output_dir}', metadata='{params.metadata}'))"
    """



##########################
### Map to Host Genome ###


# pull human genome
rule pull_host_genome:
  """
  This rule downloads the host (default human GRCh38) reference genome and annotation.
  """
  output:
    GENOME=config['host_ref_fna'], 
    ANNOTATION=config['host_ref_gtf']
  resources:
    partition=get_partition("short", config, "pull_host_genome"),
    mem_mb=get_mem(int(10*1000), config, "pull_host_genome"), # MB, or 10 GB
    runtime=get_runtime(int(1*60), config, "pull_host_genome"), # min, or 1 hr
    slurm=get_slurm_extra(config, "pull_host_genome")
  threads: get_threads(1, config, "pull_host_genome")
  params:
    ref_dir=split(config['host_ref_fna'])[0]
  shell:
    """
    mkdir -p {params.ref_dir}
    cd {params.ref_dir}

    # genome
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
    gunzip GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
    mv GCA_000001405.15_GRCh38_full_analysis_set.fna GRCh38_full_analysis_set.fna

    # annotation
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz
    gunzip GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz
    mv GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf GRCh38_full_analysis_set.refseq.gtf
    """


# make index for host genome mapping
rule build_host_genome_index:
  """
  This rule builds an index for mapping the host genome.
  """
  input:
    config['host_ref_fna']
  output:
    directory(pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "ref/"))
  conda: f"conda_envs/{get_host_map_method(config).lower()}.yaml"
  resources:
    partition=get_partition("short", config, "build_host_genome_index"),
    mem_mb=get_mem(int(32*1000), config, "build_host_genome_index"), # MB, or 30 GB
    runtime=get_runtime(int(2*60), config, "build_host_genome_index"), # min, or 2 hrs
    slurm=get_slurm_extra(config, "build_host_genome_index")
  threads: get_threads(8, config, "build_host_genome_index")
  params:
    ref_dir=f"{trim_trunc_path}.{get_host_map_method(config)}",
    method=get_host_map_method(config),
    extra=get_rule_extra_args(config, "build_host_genome_index")
  shell:
    """
    mkdir -p {params.ref_dir}
    if [ "{params.method}" == "BBMap" ]; then
        bbmap.sh ref={input} path={params.ref_dir} threads={threads} -Xmx{resources.mem_mb}m {params.extra}
    elif [ "{params.method}" == "HISAT2" ]; then
        mkdir -p {output}/
        hisat2-build {input} {output}/index -p {threads} {params.extra}
    fi
    """


# Map to host genome
rule map_host:
  """
  This rule maps reads to the host genome and compresses SAM files to BAM files.
  """
  input:
    REF=pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "ref/"),
    FWD=pj(trim_trunc_path,
          "{sample}.R1.fq"),
    REV=pj(trim_trunc_path,
          "{sample}.R2.fq")
  output:
    SAM=pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "{sample}.sam"),
    BAM=pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "{sample}.bam")
  conda: f"conda_envs/{get_host_map_method(config).lower()}.yaml"
  resources:
    partition=get_partition("short", config, "map_host"),
    mem_mb=get_mem(int(210*1000), config, "map_host"), # MB, or 210 GB, can cut to ~100 for 16 threads
    runtime=get_runtime(int(23.9*60), config, "map_host"), # min, or almost 24 hrs
    slurm=get_slurm_extra(config, "map_host")
  threads: get_threads(32, config, "map_host")
  params:
    out_dir=f"{trim_trunc_path}.{get_host_map_method(config)}",
    sam2bam_path=get_sam2bam_path(),
    method=get_host_map_method(config),
    extra=get_rule_extra_args(config, "map_host")
  shell:
    """
    cd {params.out_dir}

    if [ "{params.method}" == "BBMap" ]; then
      bbmap.sh in=../{input.FWD} in2=../{input.REV} \
      out={wildcards.sample}.sam \
      trimreaddescriptions=t \
      threads={threads} \
      -Xmx{resources.mem_mb}m \
      {params.extra}

    elif [ "{params.method}" == "HISAT2" ]; then
        hisat2 -1 ../{input.FWD} -2 ../{input.REV} \
        -S {wildcards.sample}.sam \
        -x ref/index \
        -p {threads} \
        {params.extra}
    fi
    
    bash {params.sam2bam_path} {wildcards.sample}.sam
    """

rule validate_bams:
  """
  This rule checks if BAM files are valid. 
  It doesn't stop the pipeline if they aren't but you can check these files manually.
  """
  input:
    BAM=pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "{sample}.bam")
  output:
    BAM_VALID=pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "{sample}_bam_valid.tsv")
  conda: "conda_envs/featureCounts.yaml"
  resources:
    partition=get_partition("short", config, "validate_bams"),
    mem_mb=get_mem(int(2*1000), config, "validate_bams"), # MB, or 2 GB
    runtime=get_runtime(int(1*60), config, "validate_bams"), # min, or 1 hr
    slurm=get_slurm_extra(config, "validate_bams")
  threads: get_threads(16, config, "validate_bams")
  shell:
    """
    samtools flagstat --threads {threads} {input.BAM} > {output.BAM_VALID}
    """


# Assess classification of mapped reads
rule generate_feature_counts:
  """
  This rule uses featureCounts to count host transcripts.
  """
  input:
    ANNOTATION=config['host_ref_gtf'],
    VALID=expand(pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "{sample}_bam_valid.tsv"), 
                  sample=HOST_MAP_SAMPLES),
    BAM=expand(pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "{sample}.bam"), 
                sample=HOST_MAP_SAMPLES)
  output:
    COUNTS=pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "counts.txt"),
    SUMMARY=pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "counts.txt.summary")
  conda: "conda_envs/featureCounts.yaml"
  resources:
    partition=get_partition("short", config, "generate_feature_counts"),
    mem_mb=get_mem(int(16*1000), config, "generate_feature_counts"), # MB, or 8 GB
    runtime=get_runtime(int(12*60), config, "generate_feature_counts"), # min, or 2 hrs
    slurm=get_slurm_extra(config, "generate_feature_counts")
  threads: get_threads(32, config, "generate_feature_counts")
  params:
    extra=get_rule_extra_args(config, "generate_feature_counts")
  shell:
    """
    if [ -z {input.BAM} ]
    # if there is no input (no host transcriptome mapping), just touch the files
    then
      touch {output.COUNTS}
      touch {output.SUMMARY}
    else
      featureCounts -T {threads} -p --countReadPairs \
      -t exon -g gene_id -a {input.ANNOTATION} -o {output.COUNTS} {params.extra} {input.BAM}
    fi
    """


rule feature_counts_to_tpm:
  """
  This rule converts featureCounts output to transcripts per million (TPM).
  """
  input:
    COUNTS=pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "counts.txt")
  output:
    TPM=pj(f"{trim_trunc_path}.{get_host_map_method(config)}", "counts_tpm.tsv")
  conda: "conda_envs/humann.yaml"
  resources:
    partition=get_partition("short", config, "feature_counts_to_tpm"),
    mem_mb=get_mem(int(2*1000), config, "feature_counts_to_tpm"), # MB, or 2 GB
    runtime=get_runtime(int(0.25*60), config, "feature_counts_to_tpm"), # min, should need much less
    slurm=get_slurm_extra(config, "feature_counts_to_tpm")
  threads: get_threads(1, config, "feature_counts_to_tpm")
  params:
    converter_script=get_tpm_converter_path()
  shell:
    """
    python {params.converter_script} {input.COUNTS} --output {output.TPM}
    """


rule reads_breakdown:
  """
  This rule aggregates read counts at each step, including raw reads, post-QC reads, 
  host reads, non-host reads, and non-host reads that were unmapped. 
  It outputs a CSV report of this with data per sample.
  """
  input:
    METADATA=config["METADATA"].strip(),
    GENEFAMS=pj(f"{trim_trunc_path}.nonhost.humann", 
                "all_genefamilies.tsv")
  output:
    REPORT="reads_breakdown.csv"
  conda: "conda_envs/humann.yaml"
  resources:
    partition=get_partition("short", config, "reads_breakdown"),
    mem_mb=get_mem(int(2*1000), config, "reads_breakdown"), # MB, or 2 GB
    runtime=get_runtime(int(6*60), config, "reads_breakdown"), # min, should need much less
    slurm=get_slurm_extra(config, "reads_breakdown")
  threads: get_threads(1, config, "reads_breakdown")
  params:
    reporter_script=get_read_reports_path(),
    proj=PROJ, # directory for symlinked raw reads
    hostile=f"{trim_trunc_path}.nonhost"
    # proj and hostile are more technically "inputs",
    # but putting them as params here because it doesn't know which rules
    # generate those
  shell:
    """
    python {params.reporter_script} {input.METADATA} \
    --raw_reads_dir {params.proj} \
    --hostile_dir {params.hostile} \
    --genefams_filepath {input.GENEFAMS} \
    --output {output.REPORT}

    """
