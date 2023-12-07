import pandas as pd
import numpy as np
import argparse

def filter_readcounts_df(df, sample_name=None):
    """
    filter the dataframe with readcounts
    :param df: a dataFrame contains the result coming from featureCounts
    :param sample_name: a list, all sample names, same as the result of featureCounts
    """
    result = df
    if sample_name is not None:
        return result.loc[:, sample_name].copy()

    colnames = result.columns
    dont_include = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
    sample_names = list(set(colnames) - dont_include)
    sample_reads = result.loc[:, sample_names].copy()

    return sample_reads


def get_gene_lengths(df, colname="Length"):
    return df.loc[:, ['Length']]


def read_counts2tpm(read_counts, gene_lengths):
    """
    convert read counts to TPM (transcripts per million)
    :param read_counts: a dataFrame contains the result coming from featureCounts
    :param gene_lengths: a list containing gene lengths (in the same order as read_counts)
    :return: TPM
    modified from https://gist.github.com/slowkow/c6ab0348747f86e2748b?permalink_comment_id=3051443#gistcomment-3051443
    """

    rate = read_counts.values / gene_lengths.values
    tpm = rate / np.sum(rate, axis=0).reshape(1, -1) * 1e6
    return pd.DataFrame(data=tpm, columns=read_counts.columns)


def main():

    read_counts = filter_readcounts_df(df)
    gene_lengths = get_gene_lengths(df)

    read_counts2tpm(read_counts, gene_lengths)


if __name__=="__main__":
    main()

