#! /usr/bin/python

#import FileHandler
import pysam
#import matplotlib.pyplot as plt
import numpy
import sys


class CSVFile(object):
    def __init__(self, filename):
        self.filename = filename


    def to_dict(self):
        csv_dict = {}
        with open(self.filename, 'rU') as handle:
            header, colnames = self.parse_header(handle)
            # transform colnames?
            csv_dict["header"] = header + [colnames]
            for line in handle:
                arow = line.strip('\n').split('\t')
                chrom = arow[0]
                start = int(arow[1])
                end = int(arow[2])
                summit = int(arow[4])
                v = [chrom, start, end, summit, line]
                if chrom in csv_dict:
                    csv_dict[chrom].append(v)
                else:
                    csv_dict[chrom] = [v]
        return csv_dict



    def parse_header(self, handle):
        header = []
        line = handle.readline()
        while line[0] in ('#', ','):
            header.append(line)
            line = handle.readline()
        return header, line


def build_depth_dict(chrom, start, end, samfile):
    depth_dict = {}
    for pileupcolumn in samfile.pileup(chrom, start, end):
        pos = pileupcolumn.reference_pos
        depth = len(pileupcolumn.pileups)
        depth_dict[pos] = depth
    return depth_dict


def build_depth_vec(chrom, start, end, samfile):
    depth_dict = {}
    for pileupcolumn in samfile.pileup(chrom, start, end):
        if pileupcolumn.pos > end or pileupcolumn.pos < start:
            continue
        depth_dict[pileupcolumn.pos] = len(pileupcolumn.pileups)
    depth = [0 for _ in range(end-start+1)]
    for k in depth_dict:
        try:
            depth[k-start] = depth_dict[k]
        except IndexError:
            print start, end, k
            print len(depth), k-start
            sys.stderr.write("IndexError")
            sys.exit()
    return depth


def build_model_depth(peak_vect):
    summit_depth = max(peak_vect)
    model_vect = [0 for i in range(len(peak_vect))]
    for i, v in enumerate(peak_vect):
        if v > summit_depth - 2:
            model_vect[i] == v
    return model_vect


def calc_nonmodel_val(peak, model):
    diff = list(numpy.subtract(peak, model))
    return sum(diff)/float(len(peak))


def is_tight(peak, samfile):
    chrom = peak[0]
    start = peak[1]
    end = peak[2]
    summit = peak[3]
    peak_depth = build_depth_vec(chrom, start, end, samfile)
    model_depth = build_model_depth(peak_depth)
    nonmodel = calc_nonmodel_val(peak_depth, model_depth)
    if max(peak_depth) < 5:
        return False
    else:
        return nonmodel


#def plot_hist(x):
#    x = [i for i in x if i < 4]
#    plt.hist(x, bins=50, normed=True)
#    plt.savefig("data.png")
#    plt.show()
#    plt.close()


def main(sa):
    bam_filename = sa[0]
    peak_filename = sa[1]
    samfile = pysam.AlignmentFile(bam_filename, 'rb')
    csv_dict = CSVFile(peak_filename).to_dict()
    nonmodel = []
    for line in csv_dict["header"]:
        print line.strip('\n')
    for k, v in csv_dict.items():
        if k == 'header': continue
        for peak in v:
            nonmodel = is_tight(peak, samfile)
            if nonmodel and nonmodel < 0.5:
                print peak[4].strip('\n')
    #plot_hist(nonmodel)
    samfile.close()


if __name__ == "__main__":
    main(sys.argv[1:])
