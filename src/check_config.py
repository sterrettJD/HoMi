import argparse
import pandas as pd
from HoMi_cleanup import read_config


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="path config file for pipeline (yaml)")
    parser.add_argument("metadata", help="the metadata file used for HoMi")
    return parser


def check_strings(config):
    required = {"PROJ": "PROJ is the project name and prefix given to all files created.", 
                "METADATA": "METADATA should be a filepath to the .csv file with sample metadata.", 
                "fwd_reads_path": "The name of the column in METADATA with the forward reads filepaths.", 
                "rev_reads_path": "The name of the column in METADATA with the reverse reads filepaths.",
                "hostile_db": "The hostile database that you would like to use.",
                "loc_for_hostile_db_download": "The location to store the hostile database.",
                "metaphlan_bowtie_db": "Filepath for the metaphlan database. Files will be downloaded here if not already present.",
                "chocophlan_db":  "Filepath for the HUMAnN ChocoPhlAn database. Files will be downloaded here if not already present.",
                "uniref_db":  "Filepath for the HUMANnN UniRef database. Files will be downloaded here if not already present.",
                "utility_mapping_db": "Filepath for the HUMAnN utility mapping database. Files will be downloaded here if not already present.",
                "kraken_db": "Filepath for the Kraken database. Files will be downloaded and built here if not already present.",
                "host_ref_fna": "Filepath for the host reference genome fna. If this is not present, the GRCh38 human genome will be downlaoded to this path.",
                "host_ref_gtf": "Filepath for the host reference genome gtf file. If this is not present, the GRCh38 human genome will be downlaoded to this path."}
    
    recommended = {"host_map_method": "The default host aligner (HISAT2) will be used if this is not provided."}
    
    for param in required.keys():
        conf_param = config.get(param)
        if conf_param is None:
            raise ValueError(f"{param} is missing from the config file. {required[param]}")
        if type(conf_param) != str:
            raise TypeError(f"{param} is the wrong type in the config file. It should be a string.")
    
    for param in recommended.keys():
        conf_param = config.get(param)
        if conf_param is None:
            print(f"{param} is missing from the config file. {recommended[param]}")
        if type(conf_param) != str:
            raise TypeError(f"{param} is the wrong type in the config file. It should be a string.")
    

def check_nums(config):
    required = {"min_readlen": "Reads shorter than this after trimming will be discarded.",
                "readstart_qual_min": "Base pairs below this quality will be trimmed from the beginning of each read.",
                "readend_qual_min": "Base pairs below this quality will be trimmed from the end of each read.",
                "trim_fwd": "Trim this many bases from the start of each forward read.",
                "trim_rev": "Trim this many bases from the start of each reverse read.",
                "trunc_fwd": "Trim this many bases from the end of each forward read.",
                "trunc_rev": "Trim this many bases from the end of each reverse read."}
    
    for param in required.keys():
        conf_param = config.get(param)
        if conf_param is None:
            raise ValueError(f"{param} is missing from the config file. {required[param]}")
        if type(conf_param) != int:
            raise TypeError(f"{param} is the wrong type in the config file. It should be a string.")


def check_metadata_cols(config, metadata_path):
    fwd = config["fwd_reads_path"]
    rev = config["rev_reads_path"]
    metadata = pd.read_csv(metadata_path)

    if fwd not in metadata.columns:
        raise ValueError(f"fwd_reads_path {fwd} can't be found in your metadata. Please check this.")
    if rev not in metadata.columns:
        raise ValueError(f"rev_reads_path {rev} can't be found in your metadata. Please check this.")


def main():
    args = get_args()
    config = read_config(args.config)
    check_strings(config)
    check_nums(config)
    check_metadata_cols(config, args.metadata)


if __name__=="main":
    main()
