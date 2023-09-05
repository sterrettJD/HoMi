from os import path

def hostile_db_to_path(HOSTILE_DB):
    if HOSTILE_DB in ["human-t2t-hla-argos985", "human-t2t-hla"] :
        return path.join("data", HOSTILE_DB)
    else:
        return HOSTILE_DB