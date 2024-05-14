import json
import pandas as pd
from os import path

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


def main():
    print("Not yet implemented")


if __name__ == "__main__":
    main()