'''
Detects small RNA peaks with annotated data.
'''

from collections import defaultdict
import argparse
import ConfigParser
import itertools
import os
import pysam
import re
import subprocess
import sys


def add_strands(beds, bam):
    '''
    Adds strand information to beds from BAM file.
    '''
    stranded_beds = []
    samfile = pysam.AlignmentFile(bam)
    for bed in beds:
        chrom = bed[0]
        start = int(bed[1])
        end = int(bed[2])
        reads = [read.is_reverse for read in samfile.fetch(chrom, start, end)]
        pos = reads.count(False)
        neg = reads.count(True)
        if pos < neg:
            strand = '-'
        else:
            strand = '+'
        stranded_beds.append(bed[:5] + [strand] + bed[6:] + 
                             [str(pos), str(neg)])
    return stranded_beds


def annotate_w_homer(bed, genome, config, dir_map):
    '''
    Annotates bed file with HOMER.
    '''
    homer_beds = defaultdict(list)
    homer = config.get("Binaries", "HOMER")
    for samplename, contrasts in bed.items():
        for positive, negative in contrasts:
            for bed in (positive, negative):
                homer_basename = re.sub(r'.bed',
                                        config.get('Suffix', 'homer'),
                                        bed.split('/')[-1])
                homer_bed = '/'.join([dir_map["homerdir"], homer_basename])
                path = config.get('Binaries', 'HOMER_dir')
                cmd = ' '.join(["PATH=$PATH:%s;" % path, "cut -f1-6", bed,
                                "|", homer, "-", genome, '>', homer_bed])
                if not os.path.exists(homer_bed):
                    subprocess.call(cmd, shell=True)
                else:
                    sys.stderr.write(cmd + '\n')
                homer_beds[samplename].append((homer_bed, bed))
    return homer_beds


def align_reads(fastq_map, genome, config, dir_map):
    '''
    Aligns fastq files, assumed no RG values.
    '''
    star = config.get('Binaries', 'star')
    samtools = config.get('Binaries', 'samtools')
    num_threads = config.get('Options', 'star_threads')
    ref_genome = config.get(genome, 'ref_genome')
    genome_dir = config.get(genome, 'star_dir')
    # inverts map to link the fastq to its samplename to map in another dict
    inverted_fastq_map = dict((v, k) for k in fastq_map for v in fastq_map[k])
    fastqs = list(itertools.chain.from_iterable(fastq_map.values()))
    fastq_pairs = map_fastq_pairs(fastqs)
    bams = defaultdict(list)
    for first, second in fastq_pairs.items():
        sample_name = inverted_fastq_map[first]
        bam_basename = re.sub(r'_R[1-2]_.final.clipped.fastq.gz',
                              '', first).split('/')[-1]
        basename = '/'.join([dir_map["indbamdir"], bam_basename])
        bam = '.'.join([basename, 'bam'])
        cmd1 = ' '.join([star, '--runThreadN', num_threads,
                         '--genomeDir', genome_dir,
                         '--readFilesIn', first, second,
                         '--readFilesCommand zcat',
                         '--outFileNamePrefix', basename])
        cmd2 = ' '.join(['mv', basename + 'Aligned.out.sam',
                         basename + '.sam'])
        cmd3 = ' '.join([samtools, 'view -bht', ref_genome, basename + '.sam',
                         '|', samtools, 'sort', '>', bam])
        cmd4 = ' '.join([samtools, "index", bam])
        bams[sample_name].append(bam)
        if not os.path.exists(bam):
            for cmd in (cmd1, cmd2, cmd3, cmd4):
                subprocess.call(cmd, shell=True)
    return bams


def call_macs(bams_map, contrasts, config, dir_map, optionals=""):
    '''
    Calls MACS2.
    '''
    python = config.get('Binaries', 'macs_python')
    macs = config.get('Binaries', 'macs')
    macs_dir = dir_map["macsdir"]
    peaks = defaultdict(list)
    for experiment, controls in contrasts.items():
        exp_bam = bams_map[experiment]
        for control in controls:
            positive = "%s_vs_ctrl-%s.positive%s" %\
                       (experiment, control, optionals)
            negative = "%s_vs_ctrl-%s.negative%s" %\
                       (experiment, control, optionals)
            positive_peaks = os.path.join(macs_dir, positive)
            negative_peaks = os.path.join(macs_dir, negative)
            positive_bed = positive_peaks + "_peaks.narrowPeak"
            negative_bed = negative_peaks + "_peaks.narrowPeak"
            ctrl_bam = bams_map[control]
            cmd1 = ' '.join([python, macs, "callpeak -t", exp_bam,
                             "-c", ctrl_bam, "--outdir", macs_dir,
                             "--name", positive,
                             "-g hs --keep-dup all --nomodel"])
            cmd2 = ' '.join([python, macs, "callpeak -t", ctrl_bam,
                             "-c", exp_bam, "--outdir", macs_dir,
                             "--name", negative,
                             "-g hs --keep-dup all --nomodel"])
            if not os.path.exists(positive_bed):
                subprocess.call(cmd1, shell=True)
            else:
                sys.stderr.write(cmd1 + '\n')
            if not os.path.exists(negative_bed):
                subprocess.call(cmd2, shell=True)
            else:
                sys.stderr.write(cmd2 + '\n')
            peaks[experiment].append((positive_bed, negative_bed,
                                      exp_bam, ctrl_bam))
    return peaks


def clip_fastqs(fastq_map, adapters_map, config, dir_map):
    '''
    Clips barcodes from fastqs with fastx.
    '''
    clipped_fastq_map = defaultdict(list)
    fastx = config.get('Binaries', 'fastx')
    suffix = config.get('Suffix', 'clipped')
    for samplename, fastqs in fastq_map.items():
        for fastq in fastqs:
            adapter = adapters_map[samplename]
            clipped_fastq = '/'.join([dir_map["clipfdir"],
                                      re.sub(r'\.fastq.gz',
                                             suffix + '.fastq.gz',
                                             fastq.split('/')[-1])])
            cmd = ' '.join(["zcat", fastq, "|", fastx, "-a", adapter, "|",
                            "gzip >", clipped_fastq])
            sys.stderr.write(cmd + '\n')
            if not os.path.exists(clipped_fastq):
                subprocess.call(cmd, shell=True)
            clipped_fastq_map[samplename].append(clipped_fastq)
    return clipped_fastq_map


def filter_for_squeezed_peaks(peaks, config, dir_map):
    '''
    Filters for narrow peaks in bed file.
    '''
    bedtools = config.get('Binaries', 'bedtools')
    filt_beds = defaultdict(list)
    for samplename, contrasts in peaks.items():
        for positive, negative, pos_bam, neg_bam in contrasts:
            bed_out = []
            for bed, bam in ((positive, pos_bam), (negative, neg_bam)):
                suffix_pos = re.search(r'.narrowPeak', bed).start()
                bed_basename = bed.split('/')[-1][:suffix_pos]
                filt_bed_name = '/'.join([dir_map["sizefiltmacsdir"],
                                          bed_basename + ".size_filt.bed"])
                tmp = '/'.join([dir_map["outdir"], "bedtools_coverage.tmp"])
                cmd = ' '.join([bedtools, 'coverage', '-hist',
                                '-a', bed, '-b', bam, '>', tmp])
                if not os.path.exists(filt_bed_name):
                    subprocess.call(cmd, shell=True)
                    write_bed(add_strands(parse_bedtools_hist(tmp), bam),
                              filt_bed_name)
                else:
                    sys.stderr.write(cmd + '\n')
                bed_out.append(filt_bed_name)
            filt_beds[samplename].append(bed_out)
    return filt_beds


def get_adapters(adapters_file):
    '''
    Maps fastq and adapter from TSV.
    '''
    adapters_map = {}
    with open(adapters_file, 'rU') as handle:
        for line in handle:
            arow = line.strip('\n').split(',')
            samplename = arow[1]
            adapter = arow[5]
            adapters_map[samplename] = adapter
    return adapters_map


def map_fastq_from_tsv(fastq_tsv):
    '''
    Parses tsv into map of sample name and lane specific fastq.
    '''
    fastq_map = defaultdict(list)
    with open(fastq_tsv, 'rU') as handle:
        for line in handle:
            arow = line.strip('\n').split('\t')
            sample = arow[0]
            fastq = arow[1]
            fastq_map[sample].append(fastq)
    return fastq_map


def map_fastq_pairs(fastqs):
    '''
    Gets fastqs from fastq dir and returns dict (k=first, v=second).
    '''
    first_read_fastqs = [f for f in fastqs if re.search("_R1_", f)]
    second_read_fastqs = [f for f in fastqs if re.search("_R2_", f)]
    basename_second_map = {f.replace("_R2_", "") : f
                           for f in second_read_fastqs}
    mapped_pairs = {f : basename_second_map[f.replace("_R1_", "")]
                    for f in first_read_fastqs
                    if f.replace("_R1_", "") in basename_second_map}
    first_only = {f : "" for f in first_read_fastqs
                  if f.replace("_R1_", "") not in basename_second_map}
    second_only = {f : "" for f in second_read_fastqs
                   if re.sub('_R2_', '_R1_', f) not in mapped_pairs}
    mapped_pairs.update(first_only)
    mapped_pairs.update(second_only)
    return mapped_pairs


def mask_bam(bam_map, mask_bed, suffix, config):
    '''
    Masks bam from given bed with bedtools.
    '''
    bedtools = config.get('Binaries', 'bedtools')
    masked_map = {}
    for samplename, bam in bam_map.items():
        masked_bam = re.sub(r'\.bam', suffix + '.bam', bam)
        cmd = ' '.join([bedtools, 'intersect -abam', bam, '-b', mask_bed,
                        '-wa -v >', masked_bam])
        if not os.path.exists(masked_bam):
            subprocess.call(cmd, shell=True)
        masked_map[samplename] = masked_bam
    return masked_map


def merge_bams(indv_bams_map, config, dir_map):
    '''
    Merges the individual bams into a sample bam.
    '''
    merged_bams = {}
    samtools = config.get('Binaries', 'samtools')
    for samplename, bams in indv_bams_map.items():
        input_bams = ' '.join(bams)
        output_bam = '/'.join([dir_map["bamdir"],
                               samplename + config.get('Suffix', 'bam')])
        cmd1 = ' '.join([samtools, "merge", output_bam, input_bams])
        cmd2 = ' '.join([samtools, "index", output_bam])
        if not os.path.exists(output_bam):
            subprocess.call(cmd1, shell=True)
            subprocess.call(cmd2, shell=True)
        merged_bams[samplename] = output_bam
    return merged_bams


def merge_annot_to_bed(homer_beds, config, dir_map):
    '''
    Merges data from HOMER annotated beds with MACS beds.
    '''
    final_map = defaultdict(list)
    header = ["chrom", "start", "end", "peak_name", "score", "strand",
              "fold_change", "-log10_pval", "-log10_qval", "dist_to_summit",
              "positive_reads", "negative_reads", "annotation",
              "detail_annotation", "distance_to_nearest_RefSeq_TSS",
              "Nearest_TSS-Native_ID", "Nearest_TSS-Entrez_Gene_ID",
              "Nearest_TSS-Unigene_ID", "Nearest_TSS-RefSeq_ID",
              "Nearest_TSS-Ensembl_ID", "Nearest_TSS-Gene_Symbol",
              "Nearest_TSS-Gene_Aliases", "Nearest_TSS-Gene_Description"]
    for samplename, beds in homer_beds.items():
        for homer, peak in beds:
            final_basename = re.sub(r'.bed', '.final.bed', homer.split('/')[-1])
            final_bed = '/'.join([dir_map["finaldir"], final_basename])
            homer_map = parse_bed_by_peakname(homer, 0, skip_header=True)
            peak_map = parse_bed_by_peakname(peak, 3, skip_header=False)
            if not os.path.exists(final_bed):
                with open(final_bed, 'w') as handle:
                    handle.write('\t'.join(header) + '\n')
                    for peak_name, homer_row in homer_map.items():
                        peak_row = peak_map[peak_name]
                        out_row = peak_row + homer_row[7:]
                        handle.write('\t'.join(out_row) + '\n')
            final_map[samplename].append(final_bed)
    return final_map


def parse_contrasts(filename):
    '''
    Parses two column TSV for contrasts.
    '''
    contrasts_map = defaultdict(list)
    with open(filename, 'rU') as handle:
        for line in handle:
            arow = line.strip('\n').split('\t')
            experiment = arow[0]
            control = arow[1]
            contrasts_map[experiment].append(control)
    return contrasts_map


def parse_bedtools_hist(hist_filename):
    '''
    Parses bedtools coverage --hist for regions 28-36bp & >=70% coverage.
    '''
    filtered_bed_regions = []
    peaks = {}
    with open(hist_filename, 'rU') as handle:
        for line in handle:
            if re.match("all", line):
                continue
            arow = line.strip('\n').split('\t')
            peak_id = arow[3]
            if peak_id not in peaks or int(arow[10]) > int(peaks[peak_id][10]):
                peaks[peak_id] = arow
    filtered_bed_regions = [v[:10] for v in peaks.values()
                            if 28 <= int(v[11]) <= 36]
    return filtered_bed_regions


def parse_bed_by_peakname(bed_filename, name_idx, skip_header=False):
    '''
    Parses bed file to return dict os split rows by peak name.
    '''
    row_map = {}
    with open(bed_filename, 'rU') as handle:
        if skip_header:
            handle.readline()
        for line in handle:
            arow = line.strip('\n').split('\t')
            peak_name = arow[name_idx]
            row_map[peak_name] = arow
    return row_map


def write_bed(beds, filename):
    '''
    Writes list of bed rows to file.
    '''
    with open(filename, 'w') as handle:
        for bed in beds:
            handle.write('\t'.join(bed) + '\n')


def setup_dir(cur_dir, out_dir_name):
    '''
    Sets up map of directories in output directory.
    '''
    if cur_dir == '.':
        cur_dir = os.getcwd()
    if not os.path.isdir(cur_dir):
        sys.stderr.write("Error: project directory path does not exist.\n")
        sys.exit()
    out_dir = '/'.join([cur_dir, out_dir_name])
    clipf_dir = '/'.join([out_dir, "clipped_fastqs"])
    bam_dir = '/'.join([out_dir, "bam_files"])
    indbam_dir = '/'.join([bam_dir, "individual_bam_files"])
    macs_dir = '/'.join([out_dir, 'macs_peaks'])
    size_filt_macs_dir = '/'.join([macs_dir, 'size_filtered_macs_dir'])
    homer_dir = '/'.join([out_dir, "homer_annot_dir"])
    final_dir = '/'.join([out_dir, "final_annot"])
    for folder in [out_dir, clipf_dir, bam_dir, indbam_dir, macs_dir,
                   size_filt_macs_dir, homer_dir, final_dir]:
        try:
            os.makedirs(folder)
        except OSError as err:
            sys.stderr.write("%s\n" % err)
            sys.stderr.write("Error: %s directory already exists.\n" % folder)
            sys.stderr.write("Skipping...\n")
    return {"bamdir": bam_dir,
            "clipfdir": clipf_dir,
            "outdir": out_dir,
            "projdir": cur_dir,
            "indbamdir": indbam_dir,
            "macsdir": macs_dir,
            "sizefiltmacsdir": size_filt_macs_dir,
            "homerdir": homer_dir,
            "finaldir": final_dir}


def main():
    '''
    Parses CLI args and central dispath of functionality.
    '''
    # Argument parsing
    parser = argparse.ArgumentParser(description="MACS-based sRNA peak id")
    parser.add_argument("fastqtsv", metavar="FASTQTSV",
                        help="TSV of sample name and FASTQ")
    parser.add_argument("contrasts", metavar="ContrastFile",
                        help="Two column TSV of pairs")
    parser.add_argument("-d", "--projectdir",
                        metavar="PROJECT_DIR",
                        help="project directory; default: .")
    parser.add_argument("-a", "--adapters",
                        metavar="ADAPTERS",
                        help="adapters file")
    parser.add_argument("-o", "--outdir",
                        metavar="OUTPUT_DIR",
                        help="output dir; relative to -d; default: sRNA_peaks")
    parser.add_argument("-g", "--genome",
                        metavar="GENOME",
                        help="genome; [hg19,]; default: hg19")
    parser.add_argument("--noRepeats",
                        metavar="REAPEATS_BED",
                        help="Mask repeat regions")
    parser.add_argument("--nomiR",
                        metavar="MIR_bed",
                        help="Mask miRs")
    parser.set_defaults(projectdir='.', outdir='sRNA_peaks', genome='hg19',
                        adapters=None, noRepeats=None, nomiR=None)
    args = parser.parse_args()
    # Set up globally used maps
    config = ConfigParser.ConfigParser()
    script_dir = os.path.dirname(os.path.realpath(__file__))
    config.read(os.path.join(script_dir, "config"))
    dir_map = setup_dir(args.projectdir, args.outdir)
    # Central dispatch
    contrasts_map = parse_contrasts(args.contrasts)
    fastq_maps = map_fastq_from_tsv(args.fastqtsv)
    if args.adapters:
        adapters_map = get_adapters(args.adapters)
        clipped_fastq_maps = clip_fastqs(fastq_maps, adapters_map, config,
                                         dir_map)
        fastq_maps = clipped_fastq_maps
    indv_bams_maps = align_reads(fastq_maps, args.genome, config, dir_map)
    merged_bams_map = merge_bams(indv_bams_maps, config, dir_map)
    if args.noRepeats:
        merged_bams_map = mask_bam(merged_bams_map, args.noRepeats,
                                   config.get('Suffix', 'no_repeats'), config)
    if args.nomiR:
        merged_bams_map = mask_bam(merged_bams_map, args.nomiR,
                                   config.get('Suffix', 'no_mir'), config)
    optionals = []
    if args.noRepeats:
        optionals.append(config.get('Suffix', 'no_repeats'))
    if args.nomiR:
        optionals.append(config.get('Suffix', 'no_mir'))
    optionals = ''.join(optionals)
    peaks_map = call_macs(merged_bams_map, contrasts_map, config, dir_map,
                          optionals)
    f_peaks_map = filter_for_squeezed_peaks(peaks_map, config, dir_map)
    homer_map = annotate_w_homer(f_peaks_map, args.genome, config, dir_map)
    _ = merge_annot_to_bed(homer_map, config, dir_map)
    # bedtools closest may be a very useful feature
    # may be worthwhile to use histogram coverage from bedtools coverage
    # to filter peaks


if __name__ == "__main__":
    main()
