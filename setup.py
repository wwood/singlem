from setuptools import setup
from os.path import dirname, join
import io

with open('README.md') as readme_file:
    readme = readme_file.read()


def get_version(relpath):
    """Read version info from a file without importing it"""
    for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
        if "__version__" in line:
            version_dict = eval(line.split("=")[1])
            return version_dict["singlem"]


setup(
    name='singlem',
    version=get_version("singlem/version.py"),
    description='Novelty-inclusive microbial community profiling of shotgun metagenomes',
    long_description=readme,
    long_description_content_type='text/markdown',
    url="https://github.com/wwood/SingleM",
    author='Ben Woodcroft',
    license='GPL3+',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',
        # Indicate who your project is intended for
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
    ],
    keywords="metagenomics bioinformatics",
    # Exclude test (and test data) since they takes up too much space.
    packages=['singlem','singlem.biolib_lite'],
    data_files=[(".", ["README.md", "LICENCE.txt"])],
    include_package_data=True,
    install_requires=(
        'biopython ~= 1.84',
        'extern ~= 0.4.0',
        'graftm ~= 0.15.1',
        'squarify ~= 0.4.0',
        'sqlalchemy ~= 2.0.0',
        'pandas ~= 2.2.0',
        'bird_tool_utils ~= 0.4.1',
        'pyranges ~= 0.1.0',
        'polars ~= 1.1.0',
        'tqdm ~= 4.66.0',
        'pyarrow ~= 16.1.0',
        'zenodo_backpack ~= 0.3.0',
    ),
    entry_points = {
        'console_scripts': ['singlem = singlem.main:main',
                            'lyrebird = singlem.lyrebird:main']
    },
)
