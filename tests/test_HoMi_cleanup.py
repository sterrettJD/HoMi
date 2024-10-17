import pytest
import pandas as pd
import os
from pathlib import Path

import homi_pipeline.HoMi_cleanup as hmc


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
    

def test_deletes_correct_diamond_files(tmpdir, filepath="tests/test_data/test_cleanup_config.yaml"):
    sample = "sample_1"
    config = hmc.read_config(filepath)
    prefix = f"{config['PROJ']}.f{config['trim_fwd']}.{config['trunc_fwd']}.r{config['trim_rev']}.{config['trunc_rev']}"
    humann_path = os.path.join(tmpdir, f"{prefix}.nonhost.humann")
    sample_path = os.path.join(humann_path, sample)
    sample_temp_path = os.path.join(sample_path, f"{sample}_humann_temp")

    # Create the files to be deleted
    diamond_aligned = os.path.join(sample_temp_path, f"{sample}_diamond_aligned.tsv")
    diamond_unaligned = os.path.join(sample_temp_path, f"{sample}_diamond_unaligned.fa")
    path_to_not_delete = os.path.join(sample_temp_path, "DONTDELETEME.txt")
    Path(diamond_aligned).parent.mkdir(parents=True, exist_ok=True)
    Path(diamond_aligned).touch()
    Path(diamond_unaligned).touch()
    Path(path_to_not_delete).touch()
    
    # Test that file is here
    actual_predeletion = os.listdir(sample_temp_path)
    expected_predeletion = {f"{sample}_diamond_aligned.tsv",
                            f"{sample}_diamond_unaligned.fa",
                            "DONTDELETEME.txt"}
    # Using sets so order doesn't matter
    assert set(actual_predeletion)==expected_predeletion

    # Move to temp dir
    start_dir = os.getcwd()
    os.chdir(tmpdir)
    # should delete the files
    hmc.clean_humann_temps(config, samples=[sample])
    # Test that diamond files are not here, but other files aren't deleted
    assert os.listdir(sample_temp_path)==["DONTDELETEME.txt"]
    # Go back to start
    os.chdir(start_dir)


def test_deletes_tmp_directories(tmpdir, filepath="tests/test_data/test_cleanup_config.yaml"):
    sample = "sample_1"
    config = hmc.read_config(filepath)
    prefix = f"{config['PROJ']}.f{config['trim_fwd']}.{config['trunc_fwd']}.r{config['trim_rev']}.{config['trunc_rev']}"
    humann_path = os.path.join(tmpdir, f"{prefix}.nonhost.humann")
    sample_path = os.path.join(humann_path, sample)
    sample_temp_path = os.path.join(sample_path, f"{sample}_humann_temp")

    # Create the files to be deleted
    deepest_temp = os.path.join(sample_temp_path, "tmp1283748")
    Path(deepest_temp).parent.mkdir(parents=True, exist_ok=True)
    os.mkdir(deepest_temp)
    
    # Test that temp dir is here
    assert os.listdir(sample_temp_path)==["tmp1283748"]
    # Move to temp dir
    start_dir = os.getcwd()
    os.chdir(tmpdir)
    # should delete the files
    hmc.clean_humann_temps(config, samples=[sample])
    # Test that file is not be here
    assert os.listdir(sample_temp_path)==[]
    # Go back to start
    os.chdir(start_dir)