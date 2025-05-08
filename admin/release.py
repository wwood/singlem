#!/usr/bin/env python3

import io
from os.path import dirname, join
import extern


def get_version(relpath):
    """Read version info from a file without importing it"""
    for line in io.open(join(dirname(__file__), relpath), encoding="cp437"):
        if "__version__" in line:
            if '"' in line:
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]


if __name__ == "__main__":
    version = get_version('singlem/version.py')
    print("version is {}".format(version))
    raise Exception("Is that version right? Check because version.py changed format.")

    yes_no = input(
        "Did you run the non-CI tests first, to make sure everything is OK (y/n)? \n\nmqsub -t 8 --hours 6 -- pytest --run-expensive test/test_outside_ci.py\n\n"
    )
    if yes_no != "y":
        raise Exception("Please run the non-CI tests first")

    print("building docs")
    extern.run("python3 build_docs.py")

    print(
        "Checking if repo is clean. If this fails it might be because the docs have changed from the previous command here?"
    )
    extern.run('if [[ $(git diff --shortstat 2> /dev/null | tail -n1) != "" ]]; then exit 1; fi')

    extern.run('git tag v{}'.format(version))
    print("Now run 'git push && git push --tags' and GitHub actions will build and upload to PyPI".format(version))
    print('You have to run ./build.sh from the docker directory to build the docker image, once the tag is on GitHub')
    raise Exception("Docker not setup yet for lyrebird")

    print("Once pushed to BioConda, also run https://github.com/wwood/singlem-installation to verify deployment and installation instructions")
