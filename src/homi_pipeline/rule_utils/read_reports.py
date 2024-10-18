#!/usr/bin/env python3

import json
import pandas as pd
import subprocess
import argparse
from os import path


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata", 
                        help="A csv file with sample metadata")
    parser.add_argument("--raw_reads_dir",
                        help="Path to the directory containing the raw reads")
    parser.add_argument("--hostile_dir",
                       help="Path to the directory containing hostile output reports")
    parser.add_argument("--genefams_filepath",
                        help="Filepath containing the gene families output from HUMAnN")
    parser.add_argument("--output", 
                        help="Where should the output be stored?")
    
    return parser.parse_args()


def read_hostile_report(filepath) -> dict:
    with open(filepath, 'r') as file:
        data = json.load(file)
    return data[0]


def get_reads_into_hostile(hostile_dict) -> int:
    return hostile_dict["reads_in"]


def get_nonhost_reads(hostile_dict) -> int:
    return hostile_dict["reads_out"]


def get_percent_host(hostile_dict) -> float:
    return hostile_dict["reads_removed_proportion"]


def create_df_from_hostile_reports(metadata, hostile_directory, sample_col="Sample") -> pd.DataFrame:
    df = pd.DataFrame(index=metadata[sample_col].to_list(), 
                      columns=["Reads passing QC", "Nonhost reads", "Percent host"])

    for sample in metadata[sample_col].to_list():
        report = read_hostile_report(path.join(hostile_directory, f"{sample}.report"))
        df.loc[sample, ["Reads passing QC", "Nonhost reads", "Percent host"]] = [get_reads_into_hostile(report),
                                                                                   get_nonhost_reads(report),
                                                                                   get_percent_host(report)]
    
    return df


def get_unmapped_nonhost_from_humann(genefams_filepath):
    df = pd.read_csv(genefams_filepath, sep="\t", index_col="# Gene Family")
    
    # trim off the unneeded end to each colname
    chars_to_rm = len("_Abundance-RPKs")
    df.columns = [col[:-chars_to_rm] for col in df.columns]

    # Scale from unmapped RPK to TPM
    sample_sums = df.sum(axis=0)
    unmapped_per_sample = df.loc["UNMAPPED", :]
    unmapped_tpm = (unmapped_per_sample/sample_sums)*1000000

    return unmapped_tpm.to_dict()
    

def add_unmapped_nonhost_to_hostile_reports(unmapped_nonhost: dict, hostile_reports_df: pd.DataFrame):
    # it's a bit inefficient to cast the unmapped nonhost to a dict then back to pd.Series,
    # but I think that get_unmapped_nonhost_from_humann() returning a dict makes more sense
    # as a standalone func that may be used for other things
    return hostile_reports_df.join(pd.Series(unmapped_nonhost, name="Unmapped nonhost TPM"))


def count_reads_fastq_gz(file_path) -> int:
    try:
        result = subprocess.run(f"gunzip -c {file_path} | wc -l", shell=True, capture_output=True, text=True, check=True)
        line_count = int(result.stdout.strip())
    
        # Divide by 4 to get the number of reads (assuming FASTQ format with 4 lines per read)
        read_count = line_count // 4
        return read_count
    
    except subprocess.CalledProcessError as e:
        raise ValueError(f"Error: {e}")
        return None


def add_original_read_counts_to_df(metadata, directory, df,
                                   sample_col="Sample", 
                                   forward_extension=".R1.fq.gz",
                                   reverse_extension=".R2.fq.gz") -> pd.DataFrame:
    reads_per_sample = dict()
    for sample in metadata[sample_col].to_list():
        forward_path = path.join(directory, f"{sample}{forward_extension}")
        reverse_path = path.join(directory, f"{sample}{reverse_extension}")
        reads = count_reads_fastq_gz(forward_path) + count_reads_fastq_gz(reverse_path)
        reads_per_sample[sample] = reads

    reads_per_sample_series = pd.Series(reads_per_sample, name="Raw reads")
    return df.join(reads_per_sample_series)


def main():
    args = get_args()
    metadata = pd.read_csv(args.metadata)
    
    df = create_df_from_hostile_reports(metadata, args.hostile_dir)
    
    unmapped_nonhost = get_unmapped_nonhost_from_humann(args.genefams_filepath)
    df = add_unmapped_nonhost_to_hostile_reports(unmapped_nonhost, df)

    df = add_original_read_counts_to_df(metadata, args.raw_reads_dir, df)

    df.to_csv(args.output)


if __name__ == "__main__":
    main()