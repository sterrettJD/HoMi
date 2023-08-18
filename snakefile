import pandas as pd
from os.path import join as pj

df = pd.read_csv(config['METADATA'])
SAMPLES = df["Sample"].tolist()

READS = config['READS']
PROJ = config['PROJ']

trim_trunc_path = f"{config['PROJ']}.f{config['trim_fwd']}.{config['trunc_fwd']}.r{config['trim_rev']}.{config['trunc_rev']}"

rule all:
    input: 
        # Made by SeqTK
        expand(pj(trim_trunc_path,
                  "{sample}.{read}.fq"),
                sample=SAMPLES, read=READS)
    