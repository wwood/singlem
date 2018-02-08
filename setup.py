from setuptools import setup, find_packages
import sys, os, shutil
import subprocess
import itertools

with open('README.md') as readme_file:
    readme = readme_file.read()

exec(open('singlem/version.py').read()) # loads __version__

def recursive_find(directory):
    """List the files in a directory recursively, sort of like Unix 'find'"""
    file_list = []
    for root, _, files in os.walk(directory):
        for f in files:
            file_list.append(os.path.join(root,f))
            if len(file_list) > 300:
                raise Exception("Too many files added to the recursive list")
    return file_list

# See https://stackoverflow.com/questions/20298729/pip-installing-data-files-to-the-wrong-place
# for details on how to get them working.
spkg_data_files = list([f.replace('singlem/data/','') for f in recursive_find('singlem/data')])

setup(
    name='singlem',
    version=__version__,
    description='Find de-novo operational taxonomic units (OTUs) from metagenome data',
    long_description=readme,
    url="https://github.com/wwood/SingleM",
    author='Ben Woodcroft',
    license='GPL3+',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        # Indicate who your project is intended for
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7'
    ],
    keywords="metagenomics bioinformatics",
    packages=find_packages(exclude=['contrib','docs']),
    install_requires=('graftm >= 0.10.0',
                      'extern >= 0.0.4',
                      'tempdir >= 0.6',
                      'biopython >= 1.64',
                      'dendropy >=0.4.0',
                      'pandas >= 0.19.2',
                      'biom-format >= 2.1.6',
                      'orator >= 0.9.7',
                      'squarify >= 0.3.0',
                      'matplotlib >= 2.0.2'
    ),
    setup_requires=['nose >= 1.0'],
    test_suite='nose.collector',
    scripts=['bin/singlem'],
    package_data = {'singlem.data': spkg_data_files}
)
