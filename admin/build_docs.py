#!/usr/bin/env python3

import extern
import logging
import argparse
import io
from os.path import dirname, join
import os

def remove_before(marker, string_to_process):
    splitter = '\n# ' + marker + '\n'
    if splitter not in string_to_process:
        raise Exception("Marker '{}' not found in string".format(marker))
    return splitter + string_to_process.split(splitter)[1]


def get_version(relpath):
    """Read version info from a file without importing it"""
    for line in io.open(join(dirname(__file__), '..', relpath), encoding="cp437"):
        if "__version__" in line:
            if '"' in line:
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")
    # version
    parent_parser.add_argument('--version', help='not with v e.g. 0.19.0', required=True)

    args = parent_parser.parse_args()

    # Setup logging
    debug = True
    if args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.DEBUG
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Change to parent directory, i.e. the root of the repo
    os.chdir(dirname(dirname(__file__)))

    # Update [RELEASE_TAG] in installation.md
    version = args.version
    logging.info("Updating [RELEASE_TAG] in Installation.md to {}".format(version))
    with open('docs/Installation.md.in') as f:
        installation = f.read()
    installation = installation.replace('[RELEASE_TAG]', version)
    with open('docs/Installation.md', 'w') as f:
        f.write(installation)
    logging.info("Done updating [RELEASE_TAG] in Installation.md to {}".format(version))

    subdir_and_commands = [
        ['tools', ['data','pipe','appraise','summarise','renew','supplement','prokaryotic_fraction',
                   ['lyrebird','data'], ['lyrebird','pipe']]],
        ['advanced', ['makedb','query','condense','seqs','create','metapackage','regenerate',
                      ['lyrebird','condense'], ['lyrebird','renew']]]
    ]

    for subdir, commands in subdir_and_commands:
        for subcommand in commands:
            if isinstance(subcommand, list):
                exe, subcommand = subcommand
            else:
                exe = 'singlem'
            cmd_stub = "pixi run {} {} --full-help-roff |pandoc - -t markdown-multiline_tables-simple_tables-grid_tables -f man |sed 's/\\\\\\[/[/g; s/\\\\\\]/]/g; s/^: //'".format(exe, subcommand)
            man_usage = extern.run(cmd_stub)

            if exe == 'singlem':
                subcommand_prelude = 'docs/preludes/{}_prelude.md'.format(subcommand)
                final_path = 'docs/{}/{}.md'.format(subdir, subcommand)
                title = 'SingleM'
            else:
                subcommand_prelude = 'docs/preludes/{}_{}_prelude.md'.format(exe, subcommand)
                final_path = 'docs/{}/{}_{}.md'.format(subdir, exe, subcommand)
                title = exe.capitalize()
            if os.path.exists(subcommand_prelude):
                # Remove everything before the options section
                splitters = {
                    'pipe': 'COMMON OPTIONS',
                    'prokaryotic_fraction': 'OPTIONS',
                    'data': 'OPTIONS',
                    'summarise': 'TAXONOMIC PROFILE INPUT',
                    'makedb': 'REQUIRED ARGUMENTS',
                    'appraise': 'INPUT OTU TABLE OPTIONS',
                    'seqs': 'OPTIONS',
                    'metapackage': 'OPTIONS',
                    'supplement': 'OPTIONS',
                }
                logging.info("For ROFF for command {}, removing everything before '{}'".format(
                    subcommand, splitters[subcommand]))
                man_usage = remove_before(splitters[subcommand], man_usage)

                with open(final_path, 'w') as f:
                    f.write('---\n')
                    f.write('title: {} {}\n'.format(title, subcommand))
                    f.write('---\n')
                    f.write('# {} {}\n'.format(exe, subcommand))

                    with open(subcommand_prelude) as f2:
                        f.write(f2.read())

                    f.write(man_usage)
            else:
                man_usage = remove_before('DESCRIPTION', man_usage)
                with open(final_path, 'w') as f:
                    f.write('---\n')
                    f.write('title: {} {}\n'.format(title, subcommand))
                    f.write('---\n')
                    f.write('# {} {}\n'.format(exe, subcommand))

                    f.write(man_usage)

    extern.run("doctave build")