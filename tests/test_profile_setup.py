import pytest
import os
from cookiecutter.main import cookiecutter

import src.profile_setup as ps

def test_check_name_loc(tmpdir):
    cookiecutter("gh:Snakemake-Profiles/slurm", 
                 output_dir=tmpdir, 
                 no_input=True)
    
    assert ps.check_profile_named_slurm(tmpdir) is None

    os.rename(os.path.join(tmpdir, "slurm"),
              os.path.join(tmpdir, "notslurm"))
    with pytest.raises(ValueError):
        ps.check_profile_named_slurm(tmpdir)


def test_check_conda(tmpdir):
    cookiecutter("gh:Snakemake-Profiles/slurm", 
                 output_dir=tmpdir, 
                 no_input=True)
    
    with pytest.raises(ValueError):
        ps.check_use_conda_slurm(tmpdir)