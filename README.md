# sdRNAPeaksIdentification
Various scripts for handling peak calling of short RNA.

## Requirements
* biopython (< v1.78)
* pysam (>= v0.12.0.1)
* matplotlib

## Usage

In it's simplest use, the program requires a tab delimited file for matching the sample names to the fastqs and a tab delimited file for the contrasts of samples. 

```bash
python sdRNA_peak_pipeline.py \
  tab-delim-file-for-matching-samplename-to-fastq.tsv \
  two-column-paired-contrasts.tsv;  
```

## CLI help output

```bash
usage: sdRNA_peak_pipeline.py [-h] [-d PROJECT_DIR] [-a ADAPTERS]
                              [-o OUTPUT_DIR] [-g GENOME]
                              [--noRepeats REAPEATS_BED] [--nomiR MIR_bed]
                              FASTQTSV ContrastFile

MACS-based sRNA peak id

positional arguments:
  FASTQTSV              TSV of sample name and FASTQ
  ContrastFile          Two column TSV of pairs

optional arguments:
  -h, --help            show this help message and exit
  -d PROJECT_DIR, --projectdir PROJECT_DIR
                        project directory; default: .
  -a ADAPTERS, --adapters ADAPTERS
                        adapters file
  -o OUTPUT_DIR, --outdir OUTPUT_DIR
                        output dir; relative to -d; default: sRNA_peaks
  -g GENOME, --genome GENOME
                        genome; [hg19,]; default: hg19
  --noRepeats REAPEATS_BED
                        Mask repeat regions
  --nomiR MIR_bed       Mask miRs
```
