from os import path

def hostile_db_to_path(HOSTILE_DB, parent):
    if HOSTILE_DB in ["human-t2t-hla-argos985", "human-t2t-hla"] :
        return path.join(parent, HOSTILE_DB)
    else:
        return HOSTILE_DB


def get_adapters_path():
    path_2_script = path.dirname(__file__)
    adapters_path = path.join(path_2_script, "..", "data", "adapters.fa")
    # remove ".."
    adapters_path = path.normpath(adapters_path)

    if path.exists(adapters_path):
        return adapters_path
    else:
        raise FileNotFoundError(f"Adapters were not found in {adapters_path}.")


def get_nonpareil_rmd_path():
    path_2_script = path.dirname(__file__)
    rmd_path = path.join(path_2_script, "rule_utils", "nonpareil_curves.Rmd")
    # remove ".."
    rmd_path = path.normpath(rmd_path)

    if path.exists(rmd_path):
        return rmd_path
    else:
        raise FileNotFoundError(f"Nonpareil Rmd was not found at {rmd_path}.")


def get_nonpareil_html_path():
    path_2_script = path.dirname(__file__)
    html_path = path.join(path_2_script, "rule_utils", "nonpareil_curves.html")
    # remove ".."
    html_path = path.normpath(html_path)

    return html_path


def get_agg_script_path():
    path_2_script = path.dirname(__file__)
    agg_path = path.join(path_2_script, "rule_utils", "aggregate_metaphlan_bugslists.py")
    # remove ".."
    agg_path = path.normpath(agg_path)

    if path.exists(agg_path):
        return agg_path
    else:
        raise FileNotFoundError(f"Aggregation script was not found at {agg_path}.")


def get_mphlan_conv_script_path():
    path_2_script = path.dirname(__file__)
    conv_path = path.join(path_2_script, "rule_utils", "convert_mphlan_v4_to_v3.py")
    # remove ".."
    conv_path = path.normpath(conv_path)

    if path.exists(conv_path):
        return conv_path
    else:
        raise FileNotFoundError(f"MetaPhlan conversion script was not found at {conv_path}.")


def get_taxa_barplot_rmd_path():
    path_2_script = path.dirname(__file__)
    rmd_path = path.join(path_2_script, "rule_utils", "Metaphlan_microshades.Rmd")
    # remove ".."
    rmd_path = path.normpath(rmd_path)

    if path.exists(rmd_path):
        return rmd_path
    else:
        raise FileNotFoundError(f"Taxa Barplot Rmd was not found at {rmd_path}.")


def get_func_barplot_rmd_path():
    path_2_script = path.dirname(__file__)
    rmd_path = path.join(path_2_script, "rule_utils", "HUMAnN_microshades.Rmd")
    # remove ".."
    rmd_path = path.normpath(rmd_path)

    if path.exists(rmd_path):
        return rmd_path
    else:
        raise FileNotFoundError(f"Functional Barplot Rmd was not found at {rmd_path}.")


def get_gmm_rmd_path():
    path_2_script = path.dirname(__file__)
    rmd_path = path.join(path_2_script, "rule_utils", "Gut_metabolic_modules.Rmd")
    # remove ".."
    rmd_path = path.normpath(rmd_path)

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


def get_host_mapping_samples(metadata, sample_column="Sample"):
    if "map_host" in metadata.columns:
        if metadata["map_host"].dtype != "bool":
            raise ValueError("Please provide only boolean values in the column map_host. "
                             "There cannot be anything in this column other than True/False.")
        return metadata.loc[metadata["map_host"]==True, sample_column].to_list()
    return metadata["Sample"].to_list()
