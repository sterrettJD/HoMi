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