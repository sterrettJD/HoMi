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
    