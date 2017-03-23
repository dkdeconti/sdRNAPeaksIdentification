#! /usr/bin/python

import FileHandler as fh
import PeaksFoo as pf
import pysam
import sys


'''
Filter annotations given

Input:
    peak file
    annotated peak file
    Elodie BRCA1 peaks
    PTT candidate peaks


'''


def main(sa):
    peaks_filename = sa[0]
    annot_filename = sa[1]
    brca1_filename = sa[2]
    pttpk_filename = sa[3]
    samplebam_filename = sa[4]

    '''
    peaks_dict = fh.PeaksFile().to_dict(peaks_filename)
    flt_peaks_dict = pf.filter_peaks_dict(peaks_dict)
    comparison_dict = pf.merge_dict(fh.BedFile().to_dict(brca1_filename),
                                    fh.BedFile().to_dict(pttpk_filename))
    flt_peaks_dict = pf.intersect_peaks(flt_peaks_dict,
                                        comparison_dict,
                                        buf=1000)
    peak_names = pf.transform_to_set(flt_peaks_dict)
    fh.PeaksFile.print_filtered_annot_peaks(annot_filename,
                                            peak_names)
    '''


    peaks_dict = fh.PeaksFile().to_dict(peaks_filename)
    flt_peaks_dict = pf.filter_peaks_dict(peaks_dict)
    brca_dict = fh.BedFile().to_dict(brca1_filename)
    pttpk_dict = fh.BedFile().to_dict(pttpk_filename)
    sample_bam = pysam.AlignmentFile(samplebam_filename, 'rb')

    annot_peaks_dict = pf.intersect_peaks_w_extra_annot(flt_peaks_dict,
                                                        brca_dict,
                                                        buf=1000)
    annot_peaks_dict = pf.intersect_peaks_w_extra_annot(annot_peaks_dict,
                                                        pttpk_dict,
                                                        buf=1000)
    annot_peaks_dict = pf.remove_non_match_peaks(annot_peaks_dict)
    annot_peaks_dict = pf.annot_strands(annot_peaks_dict, sample_bam)
    peak_names = pf.transform_to_dict_w_annot(annot_peaks_dict)
    fh.PeaksFile.print_filtered_annot_peaks_w_match(annot_filename,
                                                    peak_names)


if __name__ == "__main__":
    main(sys.argv[1:])
