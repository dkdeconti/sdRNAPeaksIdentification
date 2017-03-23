#!/usr/bin/python

import pysam
import sys


def skip_header(peaks_file):
    '''
    Skips the header region of the file until the line that starts with 'chr'.
    :param peaks_file:
    :return peaks_file:
    '''
    in_header = True
    while in_header:
        line = peaks_file.readline()
        sys.stdout.write(line)
        if line[:3] == "chr":
            in_header = False
    return peaks_file


def contains_depth(bam_file, chrom, begin, end):
    '''
    Checks if bam_file has any reads mapped to the given region.
    :param bam_file:
    :param chrom:
    :param begin:
    :param end:
    :return has_depth:
    '''
    has_depth = True
    iter = bam_file.fetch(chrom, begin, end)
    try:
        _ = next(iter)
    except StopIteration:
        has_depth = False
    finally:
        return has_depth


def filter_peaks(peaks_file, bam_file):
    '''
    Outputs line from peaks_file to stdout if there is no depth in the
     bam_file for that given region.
    :param peaks_file:
    :param bam_file: pysam.calignmentfile.AlignmentFile
    :return:
    '''
    peaks_file = skip_header(peaks_file)
    for line in peaks_file:
        arow = line.strip('\n').split('\t')
        chrom = arow[0]
        begin = int(arow[1])
        end = int(arow[2])
        if not contains_depth(bam_file, chrom, begin, end):
            sys.stdout.write(line)


def main(sa):
    peaks_filename = sa[-2]
    bam_filename = sa[-1]

    with open(peaks_filename, 'rU') as peaks_file, \
            pysam.AlignmentFile(bam_filename, 'rb') as bam_file:
        filter_peaks(peaks_file, bam_file)


if __name__ == "__main__":
    main(sys.argv[1:])
