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

## Help output

```bash

```
