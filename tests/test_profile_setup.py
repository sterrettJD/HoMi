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


def test_ntasks_mod(tmpdir):
    cookiecutter("gh:Snakemake-Profiles/slurm", 
                 output_dir=tmpdir, 
                 no_input=True)
    
    slurm_utils_path = os.path.join(tmpdir, "slurm", "slurm_utils.py")
    with open(slurm_utils_path, 'r') as file:
        content = file.read()
    
    # check it's there by default
    assert "options[\"cpus-per-task\"] = job_properties[\"threads\"]" in content

    # This function should replace the cpus-per-task with ntasks
    ps.convert_slurm_profile_tasks(tmpdir)

    with open(slurm_utils_path, 'r') as file:
        content = file.read()
    
    # check it's been updated
    assert "options[\"ntasks\"] = job_properties[\"threads\"]" in content
    assert "options[\"cpus-per-task\"] = job_properties[\"threads\"]" not in content
