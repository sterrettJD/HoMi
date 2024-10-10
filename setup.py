from setuptools import setup, find_packages
from pathlib import Path

__author__ = 'sterrettJD'
__version__ = "1.0.8"

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
      name="HoMi-pipeline",
      version=__version__,
      python_requires='<3.12.0', # snakemake f-strings break in 3.12 (as of Nov 15, 2023)
      setup_requires=['pytest-runner'],
      tests_require=['pytest', 'cookiecutter', 'pandas'],
      install_requires=['snakemake<8.0.0', # Temporarily keeping to snakemake 7 because 8 overhauled cluster job submission, and HoMi will need to be updated
                        'pulp==2.7.0', # higher versions of pulp don't work with snakemake https://github.com/snakemake/snakemake/issues/2606
                        'pandas', 'numpy', 
                        'biopython', 
                        'argparse', 'pyyaml', 
                        'cookiecutter'],
      scripts=['src/HoMi.py', 'src/snake_utils.py', 'src/HoMi_cleanup.py', 'src/profile_setup.py', 'src/check_config.py'],
      packages=find_packages(where="src"),
      package_dir={"": "src"},
      include_package_data=True,
      package_data={"": ["snakefile"], "data": ["adapters.fa"]},
      description="Pipeline for analysis of host-microbiome dual transcriptome data",
      long_description=long_description,
      long_description_content_type='text/markdown',
      author="John Sterrett",
      author_email='john.sterrett@colorado.edu',
      url="https://github.com/sterrettJD/HoMi",
      download_url=f"https://github.com/sterrettJD/HoMi/tarball/{__version__}"
)