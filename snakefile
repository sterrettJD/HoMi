import pandas as pd
from os.path import join as pj

METADATA = pd.read_csv(config['METADATA'])
SAMPLES = METADATA["Sample"].tolist()
RAW_FWD_READS = METADATA["forward_reads"]
RAW_REV_READS = METADATA["reverse_reads"]

READS = ["R1", "R2"]
PROJ = config['PROJ']

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
        "multiqc_report")



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

    fwd = df.loc[df["Sample"]==sample, "forward_reads"].values[0]
    rev = df.loc[df["Sample"]==sample, "reverse_reads"].values[0]

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

    conda: "conda_envs/trimmomatic"
    resources:
        mem_mb=int(8*1000), # 8 GB
        partition="short",
        runtime=int(12*60)
    threads: 8
    params:
      proj=PROJ,
      adpt=config['adapters_path'],
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

  conda: "conda_envs/fastqc"
  resources:
        partition="short",
        mem_mb=int(2*1000), # MB, or 2 GB
        runtime=int(0.5*60) # min, or 0.5 hours
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
  conda: "conda_envs/fastqc"
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

  conda: "conda_envs/seqtk"
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

  conda: "conda_envs/seqtk"
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

  conda: "conda_envs/fastqc"
  resources:
        partition="short",
        mem_mb=int(2*1000), # MB, or 2 GB
        runtime=int(0.5*60) # min, or 0.5 hours
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
  conda: "conda_envs/fastqc"
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