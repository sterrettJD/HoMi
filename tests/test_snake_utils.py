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