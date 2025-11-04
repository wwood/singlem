#!/usr/bin/env python3

import extern
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please provide the version number as an argument e.g. 0.19.0")
        sys.exit(1)
    version = sys.argv[1]

    yes_no = input(
        "Did you run the non-CI tests first, to make sure everything is OK (y/n)? \n\nmqsub -t 8 --hours 6 -- pytest --run-expensive test/test_outside_ci.py\n\n"
    )
    if yes_no != "y":
        raise Exception("Please run the non-CI tests first")

    yes_no = input(
        "Did you update CHANGELOG.md?\n\n"
    )
    if yes_no != "y":
        raise Exception("Please update the CHANGELOG.md file")
    
    print("version is {}".format(version))

    print("building docs")
    extern.run("pixi run python3 admin/build_docs.py --version {}".format(version))

    print(
        "Checking if repo is clean. If this fails it might be because the docs have changed from the previous command here? If so you need to remove the git tag with 'git tag -d v{}'".format(version)
    )
    extern.run('if [[ $(git diff --shortstat 2> /dev/null | tail -n1) != "" ]]; then exit 1; fi')

    # Generate the version file based on the git tag
    extern.run("pixi run -e dev bash -c 'SETUPTOOLS_SCM_PRETEND_VERSION={} python -m setuptools_scm --force-write-version-files'".format(version))

    print("Committing the version file")
    extern.run('git commit -a -m "v{}"'.format(version))

    print("Tagging the release as v{}".format(version))
    extern.run('git tag v{}'.format(version))
    
    print("Now run 'git push && git push --tags' and GitHub actions will build and upload to PyPI".format(version))
    print('You have to run "pixi run bash ./build.sh" from the docker directory to build the docker image, once the tag is on GitHub')
    print("Once pushed to BioConda, also run https://github.com/wwood/singlem-installation to verify deployment and installation instructions")
