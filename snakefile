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
           sample=SAMPLES)

    # Made by SeqTK
#    expand(pj(trim_trunc_path,
#              "{sample}.{read}.fq"),
#           sample=SAMPLES, read=READS)

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
    
    fwd = df.loc[df["Sample"]==sample, "forward_reads"].values[0]
    rev = df.loc[df["Sample"]==sample, "reverse_reads"].values[0]

    fwd_full = pj(cwd, fwd)
    fwd_symlink = pj(cwd, proj,f"{sample}.R1.fq.gz")

    rev_full = pj(cwd, rev)
    rev_symlink = pj(cwd, proj,f"{sample}.R2.fq.gz")


    print(f"{fwd_full} --> {fwd_symlink}")
    print(f"{rev_full} --> {rev_symlink}")
    
    symlink(fwd_full, 
            fwd_symlink)

    symlink(rev_full, 
            rev_symlink)


