#! /usr/bin/python

import math
import pysam
import sys


def annot_line(line, samfile):
    arow = line.strip('\n').split('\t')
    chrom = arow[1]
    start = int(arow[2])
    end = int(arow[3])
    reads = [read.is_reverse for read in samfile.fetch(chrom, start, end)]
    posstrand = reads.count(False)
    negstrand = reads.count(True)
    if posstrand == 0:
        strand_bias = "Inf"
    elif negstrand == 0:
        strand_bias = "-Inf"
    else:
        strand_bias = math.log(posstrand/float(negstrand), 2)
    out = arow[:4] + [str(strand_bias), str(posstrand), str(negstrand)] \
          + arow[7:]
    return '\t'.join(out)


def fix_header(line):
    arow = line.strip('\n').split('\t')
    out = arow[:4] + ["strandratio", "posstrand", "negstrand"] + arow[7:]
    return '\t'.join(out)


def parse_file(filename, samfile):
    header = True
    with open(filename) as handle:
        for line in handle:
            if header:
                sys.stdout.write(fix_header(line) + '\n')
                header = False
            else:
                sys.stdout.write(annot_line(line, samfile) + '\n')


def main(sa):
    filename = sa[0]
    samfilename = sa[1]
    samfile = pysam.AlignmentFile(samfilename)
    parse_file(filename, samfile)
    samfile.close()
    pass


if __name__ == "__main__":
    main(sys.argv[1:])
