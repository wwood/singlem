#!/usr/bin/env python2.7
import argparse
import os

from collections import OrderedDict


parser = argparse.ArgumentParser()
parser.add_argument('--arb_metadata', help='metadata file to fix', required=True)
args = parser.parse_args()


with open(args.arb_metadata) as f:
    current_record = None

    for line in f:
        l = line.strip()
        if l == 'BEGIN':
            current_record = OrderedDict()

        elif l == 'END':
            assembly_name = current_record['ncbi_genbank_assembly_accession']
            print "BEGIN"
            print 'db_name=%s' % assembly_name
            del current_record['db_name']
            for k, v in current_record.items():
                print "%s=%s" % (k,v)
            print 'END'
            print

        elif l != '':
            k, v = l.split('=')
            current_record[k] = v
