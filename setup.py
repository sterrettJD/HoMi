from setuptools import setup, find_packages
from pathlib import Path

__author__ = 'sterrettJD'
__version__ = '1.0.9'

this_directory = Path(__file__).parent
long_description = (this_directory / 'README.md').read_text()

setup(
      name='homi-pipeline',
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
      scripts=['src/homi_pipeline/HoMi.py', 
               'src/homi_pipeline/snake_utils.py', 'src/homi_pipeline/HoMi_cleanup.py', 
               'src/homi_pipeline/profile_setup.py', 'src/homi_pipeline/check_config.py',
               'src/homi_pipeline/rule_utils/combine_kreports.py', 'src/homi_pipeline/rule_utils/aggregate_metaphlan_bugslists.py',
               'src/homi_pipeline/rule_utils/convert_mphlan_v4_to_v3.py', 'src/homi_pipeline/rule_utils/counts_to_tpm.py',
               'src/homi_pipeline/rule_utils/read_reports.py'],
      packages=find_packages(where='src'),
      package_dir={'': 'src'},
      include_package_data=True,
      package_data={'homi_pipeline': ['snakefile',
                                      'data/adapters.fa',
                                    'rule_utils/sam2bam.sh',
                                   'rule_utils/Gut_metabolic_modules.Rmd', 
                                   'rule_utils/HUMAnN_to_phyloseq_helper.R', 'rule_utils/HUMAnN_microshades.Rmd', 
                                   'rule_utils/Kraken_microshades.Rmd',
                                   'rule_utils/metaphlan_to_phyloseq_helper.R','rule_utils/Metaphlan_microshades.Rmd', 
                                   'rule_utils/nonpareil_curvces.Rmd',
                                   'rule_utils/bracken_to_phyloseq_helper.R',
                                   'rule_utils/R_packages.R',
                                   'conda_envs/*.yaml']
                    },
      description='Pipeline for analysis of host-microbiome dual transcriptome data',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='John Sterrett',
      author_email='john.sterrett@colorado.edu',
      url='https://github.com/sterrettJD/HoMi',
      download_url=f'https://github.com/sterrettJD/HoMi/tarball/{__version__}'
)