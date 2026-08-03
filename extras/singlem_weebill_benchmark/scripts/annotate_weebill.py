#!/usr/bin/env python3
import argparse

from singlem.metapackage import Metapackage
from singlem.weebill import WeebillProfiler


parser = argparse.ArgumentParser()
parser.add_argument("--profile", required=True)
parser.add_argument("--metapackage", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

metapackage = Metapackage.acquire(args.metapackage)
WeebillProfiler().annotate(args.profile, metapackage, args.output)
