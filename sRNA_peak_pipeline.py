'''
Detects small RNA peaks with annotated data.
'''

from collections import defaultdict
import argparse
import ConfigParser
import os
import re
import subprocess
import sys


def align_reads(fastq_map, genome, config, dir_map):
    '''
    Aligns fastq files, assumed no RG values.
    '''
    bwa = config.get('Binaries', 'bwa')
    samtools = config.get('Binaries', 'samtools')
    ref_genome = config.get(genome, 'ref_genome')
    num_threads = config.get('Options', 'bwa_threads')
    inverted_fastq_map = dict((v, k) for k in fastq_map for v in fastq_map[k])
    fastqs = list(itertools.chain.from_iterable(fastq_map.values()))
    fastq_pairs = map_fastq_pairs(fastqs)
    bams = defaultdict(list)
    for first, second in fastq_pairs.items():
        sample_name = inverted_fastq_map[first]
        bam_basename = re.sub(r'_R[1-2]_[0-9]+\.fastq.gz',
                              '', first).split('/')[-1]
        bam = dir_map["indbamdir"] + '/' + bam_basename + \
              config.get('Suffix', 'bam')
        cmd = ' '.join([bwa, 'mem -t', num_threads,
                        ref_genome, first, second, '|',
                        samtools, "view -bht", ref_genome, "|",
                        samtools, "sort", ">", bam])
        if not os.path.exists(bam):
            subprocess.call(cmd, shell=True)
        bams[sample_name].append(bam)
    return bams


def clip_fastqs(fastq_map, adapters_map, config, dir_map):
    '''
    Clipes barcodes from fastqs with fastx.
    '''
    fastx = config.get('Binaries', 'fastx')
    for samplename, fastqs in fastq_map.items():
        for fastq in fastqs:
            adapter = adapters_map[fastq]
            clipped_fastq = '/'.join(dir_map[""])
            cmd1 = ' '.join(["zcat", fastq, "|", fastx, "-a", adapter, "|",
                             "gzip >", clipped_fastq])



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


def setup_dir(cur_dir, out_dir_name):
    '''
    Sets up map of directories in output directory.
    '''
    dir_map = {}
    return dir_map


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
    parser.add_argument("-o", "--outdir",
                        metavar="OUTPUT_DIR",
                        help="output dir; relative to -d; default: sRNA_peaks")
    parser.set_defaults(projectdir='.', outdir='sRNA_peaks', genome='hg19')
    args = parser.parse_args()
    # Set up globally used maps
    config = ConfigParser.ConfigParser()
    script_dir = os.path.dirname(os.path.realpath(__file__))
    config.read(os.path.join(script_dir, "config"))
    dir_map = setup_dir(args.projectdir, args.outdir)
    # Central dispatch
    fastq_maps = map_fastq_from_tsv(args.fastqtsv, config, dir_map)
    #clipped_fastq_maps = clip_fastqs(args.fastqtsv, config, )
    indv_bams_maps = align_reads(fastq_maps, args.genome, config, dir_map)
    merged_bams_map = merge_bams(indv_bams_maps, config, dir_map)




if __name__ == "__main__":
    main()
