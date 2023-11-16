import pytest
import pandas as pd

import src.snake_utils as su

def test_get_host_mapping_samples_nocol():
    metadata = pd.DataFrame({"Sample": [1,2,3]})
    out = su.get_host_mapping_samples(metadata)
    assert out == [1,2,3]


def test_get_host_mapping_samples_withcol():
    metadata = pd.DataFrame({"Sample": [1,2,3],
                             "map_host": [False, True, True]})
    out = su.get_host_mapping_samples(metadata)
    assert out == [2,3]


def test_get_host_mapping_samples_raises_error_nobool():
    metadata = pd.DataFrame({"Sample": [1,2,3],
                             "map_host": ["not sure", "maybe I want host", "I definitely want host here"]})
    
    with pytest.raises(ValueError):
        out = su.get_host_mapping_samples(metadata)


def test_get_slurm_extra_noparam():
    config = dict()
    assert su.get_slurm_extra(config, "test_rule") == ""


def test_get_slurm_extra_ruleparam():
    config = {"test_rule_slurm_extra": "qos=long"}
    assert su.get_slurm_extra(config, "test_rule") == "qos=long"


def test_get_slurm_extra_defaultparam():
    config = {"default_slurm_extra": "qos=long"}
    assert su.get_slurm_extra(config, "test_rule") == "qos=long"


def test_get_slurm_extra_priority():
    config = {"test_rule_slurm_extra": "qos=short", # this should be prioritized over the default
              "default_slurm_extra": "qos=long"}
    assert su.get_slurm_extra(config, "test_rule") == "qos=short"


def test_get_kraken_db_default():
    config = dict()
    assert su.get_kraken_db_loc(default="data_location", config=config) == "data_location"


def test_get_kraken_db_nondefault():
    config = {"kraken_db": "NEW_location"}
    assert su.get_kraken_db_loc(default="data_location", config=config) == "NEW_location"