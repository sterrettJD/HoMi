from setuptools import setup, find_packages
from HoMi.init import __version__ as version

__author__ = 'sterrettJD'
__version__ = version

setup(
      name="HoMi",
      version=__version__,
      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      install_requires=['snakemake', 'pandas', 'subprocess'],
      scripts=['HoMi/HoMi.py'],
      packages=find_packages(),
      description="Pipeline for analysis of host-microbiome dual transcriptome data",
      author="John Sterrett",
      author_email='john.sterrett@colorado.edu',
      url="https://github.com/sterrettJD/HoMi",
      download_url=f"https://github.com/sterrettJD/HoMi/tarball/{__version__}"
)