import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
from src.rule_utils.counts_to_tpm import read_counts2tpm

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


    tpm = read_counts2tpm(counts, gene_lengths=lengths)

    assert_frame_equal(tpm, expected_tpm)
