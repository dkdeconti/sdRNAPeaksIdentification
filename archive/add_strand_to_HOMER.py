#! /usr/bin/python

import math
import pysam
import sys
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


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
    seq_counts = {}
    seq_pos = {}
    for read in samfile.fetch(chrom, start, end):
        if read.query_sequence in seq_counts:
            seq_counts[read.query_sequence] += 1
        else:
            seq_counts[read.query_sequence] = 1
        seq_pos[read.query_sequence] = read.reference_start
    try:
        seq = max(seq_counts, key=seq_counts.get)
        seq_rc = str(Seq(seq, IUPAC.unambiguous_dna).reverse_complement())
        seq_start = seq_pos[seq]
        seq_end = seq_start + len(seq)
    except ValueError as e:
        print "ValueError:", e
        print seq_counts
        seq = "NA"
        seq_rc = "NA"
        seq_start = "NA"
        seq_end = "NA"
    except KeyError as e:
        print "KeyError:", e
        print seq_pos
        seq = "NA"
        seq_rc = "NA"
        seq_start = "NA"
        seq_end = "NA"
    out = arow[:4] + [str(strand_bias), str(posstrand), str(negstrand)] \
          + arow[7:] + [str(seq_start), str(seq_end), seq, seq_rc]
    return '\t'.join(out)


def fix_header(line):
    arow = line.strip('\n').split('\t')
    out = arow[:4] + ["strandratio", "posstrand", "negstrand", "rnaBegin", 
                      "rnaEnd", "RNAseq", "RNAseqRC"] + arow[7:]
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
    with pysam.AlignmentFile(samfilename, 'rb') as samfile:
        #samfile = pysam.AlignmentFile(samfilename)
        parse_file(filename, samfile)
        #samfile.close()


if __name__ == "__main__":
    main(sys.argv[1:])
