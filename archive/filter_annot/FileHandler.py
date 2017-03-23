import math
import sys

class BedFile(object):
    def to_dict(self, filename):
        # for parsing
        d = {}
        with open(filename, 'rU') as handle:
            for line in handle:
                if line[0] == "#":
                    continue
                k, v = self.parse_bed_line(line)
                d = self.fill_dict(k, v, d)
        return d

    def fill_dict(self, k, v, d):
        if k in d:
            d[k].append(v)
        else:
            d[k] = [v]
        return d


    def parse_bed_line(self, line):
        arow = line.strip('\n').split('\t')
        chrom = arow[0]
        start = int(arow[1])
        end = int(arow[2])
        return chrom, (start, end)


class PeaksFile(BedFile):
    def to_dict(self, filename):
        # for parsing
        d = {}
        with open(filename, 'rU') as handle:
            header = True
            for line in handle:
                if line[0] == "#" or line[0] == ",":
                    continue
                if header:
                    header = False
                    continue
                k, v = self.parse_peaks_line(line)
                d = self.fill_dict(k, v, d)
        return d


    def parse_peaks_line(self, line):
        arow = line.strip('\n').strip(',').split('\t')
        chrom = arow[0]
        start = int(arow[1])
        end = int(arow[2])
        depth = float(arow[5])
        foldchange = float(arow[7])
        q = float(arow[8])
        name = arow[9]
        return chrom, [start, end, depth, foldchange, q, name]


    @staticmethod
    def print_filtered_annot_peaks(filename, names):
        with open(filename, 'rU') as handle:
            line = handle.readline().strip('\n')
            while line[0] == "#" or line[0] == ",":
                sys.stdout.write(line + '\n')
                handle.readline()
            sys.stdout.write(line + '\n')  # print column headers
            for line in handle:
                arow = line.strip('\n').split('\t')
                peak_name = arow[0]
                if peak_name in names:
                    sys.stdout.write(line)

    @staticmethod
    def print_filtered_annot_peaks_w_match(filename, names):
        with open(filename, 'rU') as handle:
            line = handle.readline().strip('\n')
            while line[0] == "#" or line[0] == ",":
                sys.stdout.write(line + '\n')
                line = handle.readline().strip('\n')
            arow = line.split('\t')
            header = [arow[0]] + \
                ["chrom", "start", "end", "depth", "foldenrichment",
                 "qvalue", "strandratio", "+", "-", "brca1 intersects",
                 "ptt intersects"] + \
                arow[7:]
            sys.stdout.write('\t'.join(header) + '\n')  # print column headers
            for line in handle:
                arow = line.strip('\n').split('\t')
                peak_name = arow[0]
                if peak_name in names:
                    stat_annot = names[peak_name]
                    chrom = arow[1]
                    start = arow[2]
                    end = arow[3]
                    depth = stat_annot[2]
                    foldenrich = stat_annot[3]
                    qvalue = stat_annot[4]
                    posstrand = stat_annot[6][0]
                    negstrand = stat_annot[6][1]
                    if posstrand == 0:
                        strand_bias = "Inf"
                    elif negstrand == 0:
                        strand_bias = "-Inf"
                    else:
                        strand_bias = math.log(posstrand/float(negstrand), 2)
                    if stat_annot[7]:
                        brca = ','.join([str(t) for t in stat_annot[7]])
                    else:
                        brca = "None"
                    if stat_annot[8]:
                        ptt = ','.join([str(t) for t in stat_annot[8]])
                    else:
                        ptt = "None"
                    annot = arow[7:]
                    out_line = [peak_name,
                                chrom,
                                start,
                                end,
                                depth,
                                foldenrich,
                                qvalue,
                                strand_bias,
                                posstrand,
                                negstrand,
                                brca,
                                ptt] + annot
                    out_line = '\t'.join([str(foo) for foo in out_line]) + '\n'
                    sys.stdout.write(out_line)
