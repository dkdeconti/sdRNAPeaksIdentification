#! /usr/bin/python

import FileHandler
import pysam
import numpy
import sys


def build_depth_dict(chrom, start, end, samfile):
    depth_dict = {}
    for pileupcolumn in samfile.pileup(chrom, start, end):
        pos = pileupcolumn.reference_pos
        depth = len(pileupcolumn.pileups)
        depth_dict[pos] = depth
    return depth_dict


def is_tight(peak, samfile, cutoff=.2):
    chrom = peak[0]
    start = peak[1]
    end = peak[2]
    summit = peak[3]
    depth_dict = build_depth_dict(chrom, start, end, samfile)
    reads = list(samfile.fetch(chrom, start, end))
    mean_read_len = numpy.mean(reads)



def main(sa):
    bam_filename = sa[0]
    peak_filename = sa[1]
    samfile = pysam.AlignmentFile(bam_filename, 'rb')
    csv_dict = FileHandler.CSVFile(peak_filename).to_dict()
    samfile.close()


if __name__ == "__main__":
    main(sys.argv[1:])
