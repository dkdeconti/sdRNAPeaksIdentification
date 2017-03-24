'''
Detects small RNA peaks with annotated data.
'''

from collections import defaultdict
import argparse
import ConfigParser
import itertools
import os
import re
import subprocess
import sys


def align_reads(fastq_map, genome, config, dir_map):
    '''
    Aligns fastq files, assumed no RG values.
    '''
    star = config.get('Binaries', 'star')
    samtools = config.get('Binaries', 'samtools')
    num_threads = config.get('Options', 'star_threads')
    ref_genome = config.get(genome, 'ref_genome')
    genome_dir = config.get(genome, 'star_dir')

    inverted_fastq_map = dict((v, k) for k in fastq_map for v in fastq_map[k])
    fastqs = list(itertools.chain.from_iterable(fastq_map.values()))
    fastq_pairs = map_fastq_pairs(fastqs)
    bams = defaultdict(list)
    for first, second in fastq_pairs.items():
        sample_name = inverted_fastq_map[first]
        bam_basename = re.sub(r'_R[1-2]_.final.fastq.gz',
                              '.bam', first).split('/')[-1]
        basename = '/'.join([dir_map["indbamdir", bam_basename]])
        bam = '.'.join([basename, 'bam'])
        cmd1 = ' '.join([star, '--runThreadN', num_threads,
                         '--genomeDir', genome_dir,
                         '--readFilesIn', first, second,
                         '--readFilesCommand zcat',
                         '--outFileNamePrefix', basename])
        cmd2 = ' '.join(['mv', basename + '_Aligned.out.sam',
                         basename + '.sam'])
        cmd3 = ' '.join([samtools, '-bht', ref_genome, basename + '.sam',
                         '|', samtools, 'sort', '>', bam])
        bams[sample_name].append(bam)
        if not os.path.exists(bam):
            for cmd in (cmd1, cmd2, cmd3):
                subprocess.call(cmd, shell=True)
    return bams


def call_macs(bams_map, contrasts, config, dir_map):
    '''
    Calls MACS.
    '''
    for experiment, controls in contrasts.items():
        exp_bam = bams_map[experiment]
        for control in controls:
            ctrl_bam = bams_map[control]


def clip_fastqs(fastq_map, adapters_map, config, dir_map):
    '''
    Clipes barcodes from fastqs with fastx.
    '''
    clipped_fastq_map = defaultdict(list)
    fastx = config.get('Binaries', 'fastx')
    suffix = config.get('Suffix', 'clipped')
    for samplename, fastqs in fastq_map.items():
        for fastq in fastqs:
            adapter = adapters_map[fastq]
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
        masked_bam = re.sub(r'\.bam', suffix, bam)
        cmd = ' '.join([bedtools, '-abam', bam, '-b', mask_bed,
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
    pileup_dir = '/'.join([out_dir, "pileups"])
    vcf_dir = '/'.join([out_dir, "vcfs"])
    report_dir = '/'.join([out_dir, "report_html"])
    coverage_dir = '/'.join([out_dir, "coverage"])
    for folder in [out_dir, bam_dir, indbam_dir, pileup_dir,
                   coverage_dir, vcf_dir, report_dir]:
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
            "pileupdir": pileup_dir,
            "vcfdir": vcf_dir,
            "coveragedir": coverage_dir,
            "reportdir": report_dir}


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
                                   config.get('Suffix', 'no-repeats'), config)
    if args.nomiR:
        merged_bams_map = mask_bam(merged_bams_map, args.nomiR,
                                   config.get('Suffix', 'no_mir'), config)
    #optionals = args.noMir + args.noRepeats
    # TODO add optional arguments addition to name of macs to track
    # filtering, etc.
    call_macs(merged_bams_map, contrasts_map, config, dir_map, optionals)



if __name__ == "__main__":
    main()
