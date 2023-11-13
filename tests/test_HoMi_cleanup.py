import pytest
import pandas as pd
import os

import src.HoMi_cleanup as hmc


def test_config_read_correctly(filepath="tests/test_data/test_cleanup_config.yaml"):
    config = hmc.read_config(filepath)
    assert config == {
        "PROJ": "test",
        "METADATA": "tests/example_metadata.csv",
        "fwd_reads_path": "forward_reads",
        "rev_reads_path": "reverse_reads",
        "min_readlen": 50,
        "readstart_qual_min": 20,
        "readend_qual_min": 20,
        "trim_fwd": 5,
        "trim_rev": 10,
        "trunc_fwd": 0,
        "trunc_rev": 0
        }