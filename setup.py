from setuptools import setup, find_packages
from src.init import __version__ as version

__author__ = 'sterrettJD'
__version__ = version

setup(
      name="HoMi",
      version=__version__,
      python_requires='<3.12.0', # snakemake f-strings break in 3.12 (as of Nov 15, 2023)
      setup_requires=['pytest-runner'],
      tests_require=['pytest', 'cookiecutter', 'pandas'],
      install_requires=['snakemake', 'pandas', 'biopython', 'numpy', 'argparse', 'pyyaml', 'cookiecutter'],
      scripts=['src/HoMi.py', 'src/snake_utils.py', 'src/HoMi_cleanup.py', 'src/profile_setup.py'],
      packages=find_packages(),
      description="Pipeline for analysis of host-microbiome dual transcriptome data",
      author="John Sterrett",
      author_email='john.sterrett@colorado.edu',
      url="https://github.com/sterrettJD/HoMi",
      download_url=f"https://github.com/sterrettJD/HoMi/tarball/{__version__}"
)