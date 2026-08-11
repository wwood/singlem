#!/usr/bin/env python3

import argparse
import io
import logging
import os
from os.path import dirname, join

import extern


def get_version(relpath):
    """Read version info from a file without importing it."""
    for line in io.open(join(dirname(__file__), '..', relpath), encoding="cp437"):
        if "__version__" in line:
            if '"' in line:
                return line.split('"')[1]
            if "'" in line:
                return line.split("'")[1]


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--quiet', help='only output errors', action='store_true')
    parent_parser.add_argument('--version', help='not with v e.g. 0.19.0', required=True)
    args = parent_parser.parse_args()

    logging.basicConfig(
        level=logging.ERROR if args.quiet else logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')
    os.chdir(dirname(dirname(__file__)))

    for template_path, output_path in [
        ('docs/Installation.md.in', 'docs/Installation.md'),
        ('docs/SKILL.md.in', 'docs/SKILL.md'),
    ]:
        logging.info('Updating [RELEASE_TAG] in %s', output_path)
        with open(template_path) as template:
            rendered = template.read().replace('[RELEASE_TAG]', args.version)
        with open(output_path, 'w') as output:
            output.write(rendered)

    extern.run('pixi run -e docs python bin/generate_docs_reference.py')
    extern.run('pixi run -e docs zensical build --clean')
    extern.run('pixi run -e docs python bin/add_docs_redirects.py')
