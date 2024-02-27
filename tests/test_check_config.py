import pytest
import src.check_config as cc

def test_config_checker_strings(filepath="tests/test_data/test_cleanup_config.yaml"):
    config = {
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
        "trunc_rev": 0,
        "hostile_db": "database_name",
        "loc_for_hostile_db_download": "database_loc",
        "metaphlan_bowtie_db": "mphlan_loc",
        "chocophlan_db": "choco_loc",
        "uniref_db": "uniref_loc",
        "utility_mapping_db": "util_db",
        "kraken_db": "k_loc",
        "host_ref_fna": "host_fna_loc",
        "host_ref_gtf": "host_gtf_loc",
        "host_map_method": "map_method"
        }

    assert cc.check_strings(config) is None
    assert cc.check_nums(config) is None
    assert cc.run_checker(config) is None

def test_config_checker_valueerror():
    config = {
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
        "trunc_rev": 0,
        "hostile_db": "database_name",
        "loc_for_hostile_db_download": "database_loc",
        "metaphlan_bowtie_db": "mphlan_loc",
        "chocophlan_db": "choco_loc",
        "uniref_db": "uniref_loc",
        "utility_mapping_db": "util_db",
        "kraken_db": "k_loc",
        "host_ref_fna": "host_fna_loc",
        "host_ref_gtf": "host_gtf_loc",
        "host_map_method": "map_method"
        }
    
    # Remove one param at a time, should error for every required param
    for param in config.keys():
        subset = {k:v for k,v in config.items() if k != param}
        if param not in {"host_map_method"}:
            with pytest.raises(ValueError):
                cc.run_checker(subset)
        else:
            assert cc.run_checker(subset) is None

    