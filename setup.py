import os
from setuptools import setup, find_packages

def getversion():
    """Fetches version information from VERSION file"""
    with open(os.path.join('meam2nn', 'VERSION')) as version_file:
        version = version_file.read().strip()
    return version

def getreadme():
    with open('README.rst') as readme_file:
        return readme_file.read()
   
setup(name = 'meam2nn',
      version = getversion(),
      description = 'Tools for managing parameters for 2NN-MEAM potentials',
      long_description = getreadme(),
      classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Physics'
      ],
      keywords = [
        'atom', 
        'atomic', 
        'atomistic', 
        'molecular dynamics', 
        'interatomic potential',
      ], 
      url = 'https://github.com/lmhale99/meam2nn',
      author = 'Lucas Hale',
      author_email = 'lucas.hale@nist.gov',
      packages = find_packages(),
      install_requires = [
        'numpy', 
        'pandas',
      ],
      include_package_data = True,
      zip_safe = False)