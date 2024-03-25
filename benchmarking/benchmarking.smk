from semisynthetic import simulate_semisynthetic as sss

#colon_sample_htx="SRP127360"
colon_sample_htx="SRR6410603"

rule all:
    input:
        "ncbi_output"

rule fetch_data:
    output:
        directory("ncbi_output")
    threads: 8
    run:
        print(output[0])
        sss.fetch_data(colon_sample_htx, output[0], threads=threads)
