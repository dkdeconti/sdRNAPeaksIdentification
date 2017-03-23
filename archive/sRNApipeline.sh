samples=$1
adapters=$2
contrasts=$3

db_path="/cccbstore-rc/projects/cccb/pipelines/sRNApipeline"

# Clips adapters off fastq
${db_path}/preprocessing_sRNA_reads.sh \
    $samples \
    $adapters \
    ${db_path}/repeat_elements.bed \
    ${db_path}/hg19_miR.gff3

${db_path}/macs.sh $contrasts

${db_path}/postprocessing_sRNA_peaks.sh $contrasts
