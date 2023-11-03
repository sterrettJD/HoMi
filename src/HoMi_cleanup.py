import argparse
import yaml
import pandas as pd
import os
from shutil import rmtree


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="path config file for pipeline (yaml)")
    parser.add_argument("metadata", help="the metadata file used for HoMi")
    return parser.parse_args()


def read_config(config_filepath):
    with open(config_filepath, "r") as file:
        config_dict = yaml.safe_load(file)
    return config_dict


def read_sample_list(filepath):
    return pd.read_csv(filepath)["Sample"]


def main():
    args = get_args()
    config = read_config(args.config)
    samples = read_sample_list(args.metadata)





if __name__=="__main__":
    main()