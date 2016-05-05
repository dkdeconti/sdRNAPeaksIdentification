#! /usr/bin/python

import sys


class CSVFile(object):
    def __init__(self, filename):
        self.filename = filename


    def to_dict(self):
        csv_dict = {}
        with open(self.filename, 'rU') as handle:
            header, colnames = self.parse_header(handle)
            # transform colnames?
            csv_dict["heaer"] = header + [colnames]
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
