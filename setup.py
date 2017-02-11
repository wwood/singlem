from setuptools import setup, find_packages
import sys, os, shutil
import subprocess
import itertools

with open('README.md') as readme_file:
    readme = readme_file.read()

exec(open('singlem/version.py').read()) # loads __version__

# create archives during build
if 'sdist' in sys.argv:
    # Get a list of the DBs that are in 'db' of the root directory.
    from singlem.singlem import HmmDatabase # import here so it is not needed during install
    pkgs = HmmDatabase()
    archive_list = []
    for pkg in itertools.chain(*[pkgs.protein_packages(), pkgs.nucleotide_packages()]):
        # Archive the GraftM gpkg
        gpkg = pkg.graftm_package_path()
        archive = gpkg+".tar.gz"
        subprocess.check_call(["graftM", "archive", "--create", "--graftm_package",
                               gpkg, "--archive", archive, '--force'])
        # The CONTENTS file comes courtesy of the Python MANIFEST.in
        archive_list.append(archive)

if 'build' in sys.argv or 'bdist_wheel' in sys.argv:
    # Extract each of the GraftM gpkgs in 'db' into 'singlem/db'
    output_db_dir = 'singlem/db'
    if os.path.exists(output_db_dir):
        raise Exception("Please delete or move '%s' as this is the staging area for building the databases." % output_db_dir)
    os.mkdir(output_db_dir)
    spkgs = os.listdir('db')
    for spkg in spkgs:
        if spkg in ('.gitattributes'): continue
        gpkg = spkg.replace('.gpkg.spkg','')
        spkg_path = os.path.join(output_db_dir,spkg)
        os.mkdir(spkg_path)
        subprocess.check_call(
            ["graftM", "archive", "--extract",
             "--graftm_package", os.path.join(spkg_path,gpkg),
             "--archive", os.path.join(
                 'db', spkg, gpkg+".tar.gz"),
             "--force"])
        shutil.copyfile(os.path.join('db',spkg,'CONTENTS.json'),
                        os.path.join('singlem','db',spkg,'CONTENTS.json'))

def recursive_find(directory):
    """List the files in a directory recursively, sort of like Unix 'find'"""
    file_list = []
    for root, _, files in os.walk(directory):
        for f in files:
            file_list.append(os.path.join(root,f))
            if len(file_list) > 300:
                raise Exception("Too many files added to the recursive list")
    return file_list

spkg_data_files = list([f.replace('singlem/','') for f in recursive_find('singlem/db')])

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
                      'seqmagick >= 0.5.0',
                      'dendropy >=0.4.0'),
    setup_requires=['nose >= 1.0'],
    test_suite='nose.collector',
    scripts=['bin/singlem'],
    package_data = {'singlem': spkg_data_files}
)
