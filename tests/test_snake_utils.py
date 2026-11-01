import pytest
import pandas as pd

import homi_pipeline.snake_utils as su

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


def test_get_host_map_method_bbmap():
    config = {"host_map_method": "BBMap"}
    assert su.get_host_map_method(config) == "BBMap"


def test_get_host_map_method_default():
    config = {"Nothing relevant": "is here"}
    assert su.get_host_map_method(config) == "HISAT2"


def test_get_host_map_method_not_implemented():
    config = {"host_map_method": "A method that isn't implemented"}
    with pytest.raises(NotImplementedError):
        out = su.get_host_map_method(config)


def test_get_rule_extra_args_default():
    config = {"Nothing relevant": "is here"}
    assert su.get_rule_extra_args(config, "run_nonpareil") == ""


def test_get_rule_extra_args_nondefault():
    config = {"run_nonpareil_extra": "-X 500"}
    assert su.get_rule_extra_args(config, "run_nonpareil") == "-X 500"


def test_get_taxa_bar_rmd():
    assert str(su.get_taxa_barplot_rmd_path()).endswith("Metaphlan_microshades.Rmd")
    assert str(su.get_taxa_barplot_rmd_path("Kraken")).endswith("Kraken_microshades.Rmd")


def test_get_hostile_aligner():
    config = {"hostile_aligner": "hisat2"}
    assert su.get_hostile_aligner(config) == "hisat2"

    config = {"hostile_aligner": "bowtie2"}
    assert su.get_hostile_aligner(config) == "bowtie2"

    config = {"hostile_aligner": "hi"}
    with pytest.raises(ValueError):
        out = su.get_hostile_aligner(config)

    config = {"doesn't contain": "hostile aligner"}
    with pytest.raises(ValueError):
        out = su.get_hostile_aligner(config)


def test_get_hostile_exts():
    config = {"hostile_aligner": "hisat2"}
    exts = su.get_hostile_exts(config) 
    assert exts == [".1.ht2", ".2.ht2", ".3.ht2", 
                    ".4.ht2", ".5.ht2", ".6.ht2",
                    ".7.ht2", ".8.ht2"]

    config = {"hostile_aligner": "bowtie2"}
    exts = su.get_hostile_exts(config) 
    assert exts== [".1.bt2", ".2.bt2", ".3.bt2", 
                   ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]

    config = {"hostile_aligner": "hi"}
    with pytest.raises(ValueError):
        out = su.get_hostile_exts(config)


@pytest.mark.parametrize("value, expected", [
    (True, True),
    (False, False),
    ("True", True),
    ("true", True),
    ("False", False),
])
def test_is_humann_skipped(value, expected):
    config = {"skip_humann": value}
    assert su.is_humann_skipped(config) == expected


def test_is_humann_skipped_default():
    assert su.is_humann_skipped(dict()) == False


def test_get_reads_breakdown_inputs_with_humann():
    config = {"METADATA": "metadata.csv "}
    inputs = su.get_reads_breakdown_inputs(config, ["s1", "s2"], "proj.f0.0.r0.0")
    assert inputs == {"METADATA": "metadata.csv",
                      "REPORTS": ["proj.f0.0.r0.0.nonhost/s1.report",
                                  "proj.f0.0.r0.0.nonhost/s2.report"],
                      "GENEFAMS": "proj.f0.0.r0.0.nonhost.humann/all_genefamilies.tsv"}


def test_get_reads_breakdown_inputs_skip_humann():
    config = {"METADATA": "metadata.csv", "skip_humann": True}
    inputs = su.get_reads_breakdown_inputs(config, ["s1"], "proj.f0.0.r0.0")
    assert "GENEFAMS" not in inputs
    assert inputs["REPORTS"] == ["proj.f0.0.r0.0.nonhost/s1.report"]


def test_get_genefams_arg():
    assert su.get_genefams_arg(dict(), "proj") == "--genefams_filepath proj.nonhost.humann/all_genefamilies.tsv"
    assert su.get_genefams_arg({"skip_humann": True}, "proj") == ""
