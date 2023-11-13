import pytest
import pandas as pd
import os
from pathlib import Path

import src.HoMi_cleanup as hmc

# lightweight helper for making files
#def touch(fname):
#    with open(fname, 'w'):
#        pass


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
    

def test_deletes_correct_file(tmpdir, filepath="tests/test_data/test_cleanup_config.yaml"):
    sample = "sample_1"
    config = hmc.read_config(filepath)
    prefix = f"{config['PROJ']}.f{config['trim_fwd']}.{config['trunc_fwd']}.r{config['trim_rev']}.{config['trunc_rev']}"
    humann_path = os.path.join(tmpdir, f"{prefix}.nonhost.humann")
    sample_path = os.path.join(humann_path, sample)
    sample_temp_path = os.path.join(sample_path, f"{sample}_humann_temp")

    diamond_aligned = os.path.join(sample_temp_path, f"{sample}_diamond_aligned.tsv")
    Path(diamond_aligned).parent.mkdir(parents=True, exist_ok=True)
    Path(diamond_aligned).touch()
    
    # Test that file is here
    assert os.listdir(sample_temp_path)==[f"{sample}_diamond_aligned.tsv"]
    # Move to temp dir
    start_dir = os.getcwd()
    os.chdir(tmpdir)
    # should delete the files
    hmc.clean_humann_temps(config, samples=[sample])
    # Test that file is not be here
    assert os.listdir(sample_temp_path)==[]
    # Go back to start
    os.chdir(start_dir)