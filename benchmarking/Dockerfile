# base of Docker image is Ubuntu version 24.04
# mostly built to keep a stable version of bash/shell 
FROM ubuntu:24.04

# update OS libraries and install OS dependencies to get FastQC working
RUN apt-get update &&  \
    apt-get install -y unzip wget perl default-jre default-jdk rsync bc vim && \
    rm -rf /var/lib/apt/lists/*  # Clean up to reduce image size

# install miniforge3 independently 
RUN cd /opt/ && \
    wget "https://github.com/conda-forge/miniforge/releases/download/24.11.3-2/Miniforge3-24.11.3-2-$(uname)-$(uname -m).sh" && \
    chmod a+x Miniforge3-24.11.3-2-$(uname)-$(uname -m).sh && \
    # install in batch mode so not prompted for user input
    bash Miniforge3-24.11.3-2-$(uname)-$(uname -m).sh -b -p /opt/miniforge3 && \ 
    rm Miniforge3-24.11.3-2-$(uname)-$(uname -m).sh

# install sra toolkit for ubuntu and delete tar.gz file to save space
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.3.0/sratoolkit.3.3.0-ubuntu64.tar.gz && \
    tar -xvzf sratoolkit.3.3.0-ubuntu64.tar.gz && \
    rm sratoolkit.3.3.0-ubuntu64.tar.gz

# set working directory to benchmarking
WORKDIR /HoMi/benchmarking

# copy benchmarking conda env yaml file to container
COPY benchmarking_env.yaml benchmarking_env.yaml

# initialize conda and install benchmarking conda env
RUN /opt/miniforge3/bin/conda env create -f benchmarking_env.yaml && \
    rm benchmarking_env.yaml && \
    /opt/miniforge3/bin/conda clean --all --yes # clean up conda cache to reduce image size

# install the HoMi pipeline
# NOT NEEDED BECAUSE THIS IS INSTALLED IN THE CONDA ENVIRONMENT
# RUN /opt/miniforge3/envs/homi_benchmarking/bin/pip install -e /HoMi

# set up command in .bashrc to automatically activate benchmarking conda env when container is opened
RUN echo "source activate homi_benchmarking" >> ~/.bashrc

## do i actually still need the basic miniforge3 bin in the path??
ENV PATH="/opt/miniforge3/envs/homi_benchmarking/bin:/opt/miniforge3/bin:/opt/sratoolkit.3.3.0-ubuntu64/bin/:$PATH"