from importlib.resources import files
import shutil
from os import path

def hostile_db_to_path(HOSTILE_DB, parent):
    if HOSTILE_DB in ["human-t2t-hla-argos985", "human-t2t-hla"] :
        return path.join(parent, HOSTILE_DB)
    else:
        return HOSTILE_DB


def get_adapters_path():
    adapters_path = files("homi_pipeline").joinpath("data/adapters.fa")
    if path.exists(adapters_path):
        return adapters_path
    else:
        raise FileNotFoundError(f"Adapters were not found in {adapters_path}.")


def get_nonpareil_rmd_path():
    rmd_path = files("rule_utils").joinpath("nonpareil_curves.Rmd")

    if path.exists(rmd_path):
        return rmd_path
    else:
        raise FileNotFoundError(f"Nonpareil Rmd was not found at {rmd_path}.")


def get_nonpareil_html_path():
    html_path = files("homi_pipeline").joinpath("rule_utils/nonpareil_curves.html")
    return html_path


def get_agg_script_path():
    agg_path = shutil.which("aggregate_metaphlan_bugslists.py")
    if path.exists(agg_path):
        return agg_path
    else:
        raise FileNotFoundError(f"Aggregation script was not found at {agg_path}.")


def get_mphlan_conv_script_path():
    conv_path = shutil.which("convert_mphlan_v4_to_v3.py")

    if path.exists(conv_path):
        return conv_path
    else:
        raise FileNotFoundError(f"MetaPhlan conversion script was not found at {conv_path}.")


def get_taxa_barplot_rmd_path(map_method="Metaphlan"):
    rmd_path = files("homi_pipeline").joinpath(f"rule_utils/{map_method}_microshades.Rmd")
    
    if path.exists(rmd_path):
        return rmd_path
    else:
        raise FileNotFoundError(f"Taxa Barplot Rmd was not found at {rmd_path}.")


def get_func_barplot_rmd_path():
    rmd_path = files("homi_pipeline").joinpath("rule_utils/HUMAnN_microshades.Rmd")

    if path.exists(rmd_path):
        return rmd_path
    else:
        raise FileNotFoundError(f"Functional Barplot Rmd was not found at {rmd_path}.")


def get_gmm_rmd_path():
    rmd_path = files("homi_pipeline").joinpath("rule_utils/Gut_metabolic_modules.Rmd")

    if path.exists(rmd_path):
        return rmd_path
    else:
        raise FileNotFoundError(f"Gut metabolic modules Rmd was not found at {rmd_path}.")


def get_sam2bam_path():
    path_2_script = path.dirname(__file__)
    sam2bam_path = path.join(path_2_script, "rule_utils", "sam2bam.sh")
    # remove ".."
    sam2bam_path = path.normpath(sam2bam_path)

    if path.exists(sam2bam_path):
        return sam2bam_path
    else:
        raise FileNotFoundError(f"Sam2bam script was not found at {sam2bam_path}.")


def get_partition(default, config, rule_name):
    # First, check the config for a rule-specific partition
    config_param = config.get(f"{rule_name}_partition")
    if config_param is not None:
        return config_param
    # Then, check the config for a different default partition name
    # This looks for default_short_partition_name or default_long_partition_name
    # in case users need a different name for these short and long partitions
    config_diff_default = config.get(f"default_{default}_partition_name")
    if config_diff_default is not None:
        return config_diff_default
    
    # If there's nothing in the config, use what's in the snakefile
    return default


def get_mem(default, config, rule_name):
    config_param = config.get(f"{rule_name}_mem_mb")
    if config_param is not None:
        return int(config_param)
    return default


def get_runtime(default, config, rule_name):
    config_param = config.get(f"{rule_name}_runtime")
    if config_param is not None:
        return int(config_param)
    return default


def get_threads(default, config, rule_name):
    config_param = config.get(f"{rule_name}_threads")
    if config_param is not None:
        return int(config_param)
    return default


def get_slurm_extra(config, rule_name):
    config_param = config.get(f"{rule_name}_slurm_extra")
    if config_param is not None:
        return config_param
    
    config_diff_default = config.get(f"default_slurm_extra")
    if config_diff_default is not None:
        return config_diff_default
    
    return ""


def get_host_mapping_samples(metadata, sample_column="Sample"):
    if "map_host" in metadata.columns:
        if metadata["map_host"].dtype != "bool":
            raise ValueError("Please provide only boolean values in the column map_host. "
                             "There cannot be anything in this column other than True/False.")
        return metadata.loc[metadata["map_host"]==True, sample_column].to_list()
    return metadata["Sample"].to_list()


def get_kraken_db_loc(default, config):
    loc = config.get("kraken_db")
    if loc is not None:
        return loc
    return default

def get_tpm_converter_path():
    converter = shutil.which("counts_to_tpm.py")

    if path.exists(converter):
        return converter
    else:
        raise FileNotFoundError(f"TPM conversion script was not found at {converter}.")
    
def get_host_map_method(config):
    implemented_mappers = ["HISAT2", "BBMap"]
    map_method = config.get("host_map_method")
    if map_method is None:
        return "HISAT2"
    if map_method in implemented_mappers:
        return map_method
    raise NotImplementedError(f"Mapping host transcriptomes with {map_method} is not yet an implemented option. Please use an option from {implemented_mappers}")


def get_rule_extra_args(config, rule_name):
    """
    This looks for a config param named <rule_name>_extra.
    If that exists, it's returned. 
    Otherwise, a blank string is returned.
    """
    return config.get(f"{rule_name}_extra", "")


def get_metaphlan_index_name(config):
    """
    Checks the config for the index name. If there isn't one, it returns "latest".
    """
    # See if one is provided
    if config.get("metaphlan_index_name") is not None:
        return config.get("metaphlan_index_name")
    
    return "latest"
    

def read_latest_metaphlan_index_name(config):
    """
    This uses the mpa_latest file to get the database name for metaphlan.
    """
    # If one isn't provided, read it from mpa_latest
    loc = config["metaphlan_bowtie_db"]
    with open(path.join(loc, "mpa_latest")) as f:
        name = f.read()
    return name


def get_R_installation_path():
    r_path = files("homi_pipeline").joinpath("rule_utils/R_packages.R")
    
    if path.exists(r_path):
        return r_path
    else:
        raise FileNotFoundError(f"R package installation script was not found at {r_path}.")
    

def get_read_reports_path():
    reporter = shutil.which("read_reports.py")
    
    if path.exists(reporter):
        return reporter
    else:
        raise FileNotFoundError(f"Read reports script was not found at {reporter}.")