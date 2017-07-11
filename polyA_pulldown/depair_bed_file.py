#!/usr/bin/python
'''
Depairs the paired read of appropriately oriented adapter seq.
'''

import argparse
import sys


def filter_read_names(filename, names):
    '''
    Parses file and filters read column by given set.
    '''
    with open(filename, 'rU') as handle:
        for line in handle:
            if line.strip('\n').split('\t')[3] not in names:
                sys.stdout.write(line)


def parse_for_read_names(filename):
    '''
    Parses file for names column to return a set of names.
    '''
    read_names = set([])
    with open(filename, 'rU') as handle:
        for line in handle:
            read_name = line.strip('\n').split('\t')[3]
            read_names.add(line.strip('\n').split('\t')[3])
    return read_names


def main():
    '''
    Arg parsing and central dispatch.
    '''
    # Arg parsing
    parser = argparse.ArgumentParser(description="Depairs BED files")
    parser.add_argument("bed", metavar="BED", help="Main BED file.")
    parser.add_argument("pair", metavar="PAIRED_BED",
                        help="Paired bed to filter.")
    args = parser.parse_args()
    # Function dispatch
    filter_read_names(args.pair, parse_for_read_names(args.bed))


if __name__ == "__main__":
    main()
