#! /usr/bin/python

import sys


def parse_header(handle):
    header = handle.readline().strip('\n').strip(',')
    if header[0] == "#":
        while len(header) == 0 or header[0] == "#":
            header = handle.readline().strip('\n').strip(',')
    sys.stdout.write('#' + header + '\n')


def parse_line(line):
    arow = line.strip('\n').strip(',').split('\t')
    chrom = arow[0]
    start = arow[1]
    end = arow[2]
    peakid = arow[3]
    fc = arow[6]
    p = arow[7]
    q = arow[8]
    strand = arow[10]
    unused = "."
    #return [chrom, start, end, peakid, unused, strand, fc, p, q]
    return [chrom, start, end, peakid, unused, strand]


def parse_body(handle):
    sys.stdout.write('\t'.join(["#chr",
                                "start",
                                "end", 
                                "peak_id",
                                "unused",
                                "strand"]) + '\n')
    for line in handle:
        out = parse_line(line)
        sys.stdout.write('\t'.join(out) + '\n')


def convert_file(filename):
    with open(filename, 'rU') as handle:
        #parse_header(handle)
        #if handle:
        parse_body(handle)


def main(sa):
    filename = sa[0]
    convert_file(filename)


if __name__ == "__main__":
    main(sys.argv[1:])
