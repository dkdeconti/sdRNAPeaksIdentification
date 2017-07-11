#! /usr/bin/python

import argparse
#import pysam
import re
import sys


def parse_annot(annot_file):
    '''
    Parses annotation file to matrix with separate header.
    '''
    annot_list = []
    with open(annot_file, 'rU') as handle:
        header = handle.readline().strip('\n').split('\t')
        for line in handle:
            if re.match("ValueError", line) or re.match("{}", line):
                continue
            annot_list.append(line.strip('\n').split('\t'))
    return header, annot_list


def parse_peaks(peakfile):
    '''
    Parses peak files to map by peak id.
    '''
    peaks_map = {}
    with open(peakfile, 'rU') as handle:
        for line in handle:
            arow = line.strip('\n').split('\t')
            peak_id = arow[3]
            fc = arow[6]
            pval = arow[7]
            qval = arow[8]
            peaks_map[peak_id] = (fc, pval, qval)
    return peaks_map


def re_annot(annots, peaks_map):
    '''
    Adds new annotations to annot file.
    '''
    #new_annots = []
    for annot in annots:
        peak_id = annot[0]
        fc, pval, qval = peaks_map[peak_id]
        #new_annots.append(annot + [fc, pval, qval])
        sys.stdout.write('\t'.join(annot + [fc, pval, qval]) + '\n')


def main():
    '''
    Parses CLI args and central dispatch.
    '''
    # Arg parsing
    parser = argparse.ArgumentParser(description="Addition annotation")
    parser.add_argument("annot", metavar="ANNOT",
                        help="Previously annotated files")
    parser.add_argument("peaks", metavar="MACSPEAK",
                        help="original MACS2 peak output")
    #parser.add_argument("bam", metavar="BAM",
    #                    help="BAM file of non-control/input")
    args = parser.parse_args()
    # Central dispatch
    #samfile = pysam.AlignmentFile(args.bam)
    header, annot_list = parse_annot(args.annot)
    peak_map = parse_peaks(args.peaks)
    sys.stdout.write('\t'.join(header + ["fc", "pval", "qval"]) + '\n')
    re_annot(annot_list, peak_map)
    # Wrap up files
    #samfile.close()


if __name__ == "__main__":
    main()
