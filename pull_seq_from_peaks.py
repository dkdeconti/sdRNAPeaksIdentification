#!/usr/bin/python

import pysam
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO


def get_seq(samfile, chrom, begin, end):
    seq_counts = {}
    for read in samfile.fetch(chrom, begin, end):
        if read.query_sequence in seq_counts:
            seq_counts[read.query_sequence] += 1
        else:
            seq_counts[read.query_sequence] = 1
    return max(seq_counts, key=seq_counts.get)


def read_bed(filename):
    bed_list = []
    with open(filename, 'rU') as handle:
        for line in handle:
            arow = line.strip('\n').split()
            chrom = arow[0]
            try:
                begin = int(arow[1])
                end = int(arow[2])
            except ValueError as e:
                sys.stderr.write("ValueError: " + str(e) + "\n")
                sys.exit()
            bed_list.append([chrom, begin, end])
    return bed_list


def main(sa):
    target_region_filename = sa[0]
    bam_filename = sa[1]
    list_of_regions = read_bed((target_region_filename))
    seq_list = []
    with pysam.AlignmentFile(bam_filename, 'rb') as samfile:
        for chrom, begin, end in list_of_regions:
            seq_list.append(get_seq(samfile, chrom, begin, end))
    targets = []
    targets_rc = []
    for i, s in enumerate(seq_list):
        s = Seq(s, IUPAC.unambiguous_dna)
        rc = s.reverse_complement()
        targets.append(SeqRecord(s, id="-".join([str(j) for j in list_of_regions[i]])))
        targets_rc.append(SeqRecord(rc, id="-".join([str(j) for j in list_of_regions[i]] + ["rc"])))
    targets_rc = [s.reverse_complement() for s in targets]
    SeqIO.write(targets, "target_seq.fa", "fasta")
    SeqIO.write(targets_rc, "reverse_complement.fa", "fasta")


if __name__ == "__main__":
    main(sys.argv[1:])
