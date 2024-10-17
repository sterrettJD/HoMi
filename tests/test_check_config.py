import pytest
import homi_pipeline.check_config as cc

@pytest.fixture
def config():
    return {
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


def test_config_checker_strings_works(config):
    assert cc.check_strings(config) is None


def test_config_checker_nums_works(config):
    assert cc.check_nums(config) is None
    
    
def test_config_checker_works(config):
    assert cc.run_checker(config) is None


def test_config_checker_valueerror(config):
    # Remove one param at a time, should error for every required param
    for param in config.keys():
        subset = {k:v for k,v in config.items() if k != param}
        if param not in {"host_map_method"}:
            with pytest.raises(ValueError):
                cc.run_checker(subset)
        else:
            assert cc.run_checker(subset) is None


def test_config_checker_num_typeerror(config):
    # Swap each num for a string. Should error
    for param in config.keys():
        current_config = config.copy()
        if type(current_config[param]) == int:
            current_config[param] = "A STRING WHERE IT SHOULDN'T BE"
            with pytest.raises(TypeError):
                cc.check_nums(current_config)


def test_config_checker_string_typeerror(config):
    # Swap each num for a string. Should error
    for param in config.keys():
        current_config = config.copy()
        if type(current_config[param]) == str:
            # Put a numeric type where it shouldn't be
            current_config[param] = 1
            with pytest.raises(TypeError):
                cc.check_strings(current_config)