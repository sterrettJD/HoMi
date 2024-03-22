# Benchmarking HoMi

This directory includes helper scripts to benchmark HoMi against simulated fully synthetic and semi-synthetic communities.

## Environment
`benchmarking_env.yaml` contains the conda environment used for benchmarking. To create this environment, run:
```
conda create -n benchmarking-homi -f benchmarking_env.yaml
```

NOTE: Because `ncbi-datasets-cli` doesn't have a MacOS ARM-compatible conda-forge distribution, MacOS ARM users need to create the environment and configure the subdir with the following commands:
```
conda create -n benchmarking-homi
conda activate benchmarking-homi 
conda config --env --set subdir osx-64
conda env update --file benchmarking_env.yaml
```

## Semi-synthetic communities

