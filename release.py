#!/usr/bin/env python3

import extern
import logging
import argparse

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    # parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    args = parent_parser.parse_args()

    # Setup logging
    debug = True
    if args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.DEBUG
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    for subcommand in ['data','pipe','appraise','makedb','query','condense','summarise','renew',]:
        extern.run("bin/singlem {} --full-help-roff |pandoc - -t markdown-multiline_tables-simple_tables-grid_tables -f man |sed 's/\\\\\\[/[/g; s/\\\\\\]/]/g' |cat <(sed s/SUBCOMMAND/{}/ prelude) - >docs/usage/{}.md".format(subcommand,subcommand,subcommand))

    extern.run("doctave build")

