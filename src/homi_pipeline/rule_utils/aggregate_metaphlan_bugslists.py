#!/usr/bin/env python3

import os
import argparse
import pandas as pd


def get_args():
    """
    handles arg parsing for this script

    returns the parsed args
    """
    parser = argparse.ArgumentParser(
        prog="Aggregate Metaphlan bugs lists",
        description="Aggregates metaphlan bugs lists from humann outputs for multiple samples into a single tsv"
    )
    # input directory
    parser.add_argument("-i", "--indir",
                        help="Input directory that contains all of the sample directories",
                        required=True)
    # a "padder" string that exists between the sample id and the metaphlan output name type
    parser.add_argument("-p", "--padder",
                        help="a 'padder' string that exists between the sample id and the humann output name type",
                        default="")
    # Output file path
    parser.add_argument("-o", "--outfile",
                        help="Output file path",
                        required=True)

    parsed_args = parser.parse_args()
    return parsed_args


def get_filepaths(directory, padder):
    # get list of only subdirectories
    print(f"Searching {directory}")
    subdirs = [d for d in os.listdir(directory) 
               if os.path.isdir(os.path.join(directory, d))]
    # and each of these contains a subdir sampleid.concat_humann_temp
    # and each of the bugs list is named sampleid.concat_metaphlan_bugs_list.tsv
    filepaths = [os.path.join(directory, sampid,
                              f"{sampid}{padder}_humann_temp",
                              f"{sampid}{padder}_metaphlan_bugs_list.tsv") for sampid in subdirs]
    return filepaths, subdirs


def concat_files(filepaths_list, sampids_list):
    for i, f in enumerate(filepaths_list):
        sample_id = sampids_list[i]

        df = pd.read_csv(f, sep="\t", header=4)
        # set index to be the clade name so we can merge the dataframes
        df.index = df["#clade_name"]
        df.rename(columns={"relative_abundance": sample_id}, inplace=True)
        # create our full df if it doesn't exist yet
        if i == 0:
            full = df[["NCBI_tax_id", "additional_species"]]
        # add the new sample to the full df
        full = pd.concat([full, df[sample_id]], axis=1)

    return full

if __name__ == "__main__":
    args = get_args()

    filepaths_list, sampids_list = get_filepaths(args.indir,
                                                 padder=args.padder)
    print("Aggregating the following files: ")
    for filepath in filepaths_list:
        print(filepath)
    
    full = concat_files(filepaths_list, sampids_list)

    full.to_csv(f"{args.outfile}", sep="\t")
