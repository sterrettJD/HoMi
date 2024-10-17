import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
from homi_pipeline.rule_utils import counts_to_tpm as c2t

def test_tpm_conversion():
    counts = pd.DataFrame(data = {
            "S1": [150, 50],
            "S2": [100, 100],
            "S3": [3, 1]
        })
    lengths = pd.Series([150, 50])
    
    expected_tpm = pd.DataFrame(data={
        "S1": [500000.0, 500000.0],
        "S2": [250000.0, 750000.0],
        "S3": [500000.0, 500000.0]})


    tpm = c2t.read_counts2tpm(counts, gene_lengths=lengths)

    assert_frame_equal(tpm, expected_tpm)

# Sample data for testing
@pytest.fixture
def sample_featurecounts_df():
    data = {
        "Geneid": ["Gene1", "Gene2"],
        "Chr": ["chr1", "chr1"],
        "Start": [0, 200],
        "End": [150, 250],
        "Strand": ["+", "-"],
        "Length": [150, 50],
        "S1": [150, 50],
        "S2": [100, 100],
        "S3": [3, 1]
    }
    return pd.DataFrame(data)


# Test filter_readcounts_df function
def test_filter_readcounts_df_with_sample_name(sample_featurecounts_df):
    # Test when sample_name is specified
    sample_names = ["S1", "S2"]
    result = c2t.filter_readcounts_df(sample_featurecounts_df, sample_name=sample_names)

    # Assert that only the specified columns are present
    assert set(result.columns) == set(sample_names)


def test_filter_readcounts_df_without_sample_name(sample_featurecounts_df):
    # Test when sample_name is not specified
    result = c2t.filter_readcounts_df(sample_featurecounts_df)

    # Assert that non-sample columns are excluded
    excluded_columns = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
    assert set(result.columns) == set(sample_featurecounts_df.columns) - excluded_columns


# Test get_gene_lengths function
def test_get_gene_lengths(sample_featurecounts_df):
    result = c2t.get_gene_lengths(sample_featurecounts_df)

    # Assert that the returned column is the "Length" column
    assert result.equals(sample_featurecounts_df["Length"])


# Test get_gene_names function
def test_get_gene_names(sample_featurecounts_df):
    result = c2t.get_gene_names(sample_featurecounts_df)

    # Assert that the returned column is the "Geneid" column
    assert result.equals(sample_featurecounts_df["Geneid"])

def test_full_conv(sample_featurecounts_df):
    tpm = c2t.convert_dataframe(sample_featurecounts_df)

    expected = pd.DataFrame({
        "S1": [500000.0, 500000.0],
        "S2": [250000.0, 750000.0],
        "S3": [500000.0, 500000.0],
    })

    expected.index = sample_featurecounts_df["Geneid"]

    assert_frame_equal(tpm, expected)