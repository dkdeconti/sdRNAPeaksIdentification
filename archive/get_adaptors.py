#! /usr/bin/python

import re
import sys


def main(sa):
    filename = sa[0]
    csv_filename = sa[1]
    with open(csv_filename) as handle:
        for line in handle:
            arow = line.strip('\n').split(',')
            if re.search(arow[1], filename):
                print arow[-3]
                break


if __name__ == "__main__":
    main(sys.argv[1:])
