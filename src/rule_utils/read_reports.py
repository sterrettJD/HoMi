import json
import pandas as pd
import subprocess
from os import path

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

    return df.loc["UNMAPPED", :].to_dict()
    

def add_unmapped_nonhost_to_hostile_reports(unmapped_nonhost: dict, hostile_reports_df: pd.DataFrame):
    # it's a bit inefficient to cast the unmapped nonhost to a dict then back to pd.Series,
    # but I think that get_unmapped_nonhost_from_humann() returning a dict makes more sense
    # as a standalone func that may be used for other things
    return hostile_reports_df.join(pd.Series(unmapped_nonhost, name="Unmapped nonhost RPK"))


def main():
    print("Not yet implemented")


if __name__ == "__main__":
    main()