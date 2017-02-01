##############################################################################
# Reads go in -> bam files come out
# 1. fastx clipper
# 2. star align
# 3. sam -> bam
# 4. filter repeat reginos
# 5. filter miR
##############################################################################
# cli args
samples=$1;
adapters=$2;
repeats=$3;
miRDB=$4;

# binary paths
apps_dir="/cccbstore-rc/projects/cccb/apps";
bin_libs="/cccbstore-rc/projects/cccb/pipelines/sRNApipeline";
fastx="${apps_dir}/fastx_toolkit-0.0.14/bin/fastx_clipper";
star="${apps_dir}/STAR_2.5.0/bin/Linux_x86_64_static/STAR"
bedtools="${apps_dir}/bedtools2-2.22.1/bin/intersectBed"
java="${apps_dir}/jre1.8.0_74/bin/java"
picardtools="${apps_dir}/picard-tools-2.1.1/picard.jar"

# required files paths
db_dir="/cccbstore-rc/projects/db/gatk/hg19"
genome_dir=/cccbstore-rc/projects/db/gatk/hg19/STAR

# make directory for clipped fastq
if [ ! -d ./fastq_clipped ]; then
    mkdir ./fastq_clipped;
fi;

# Identify adapter sequence and use fastx_clipper to clipp adapter sequence
# gzips output fastq
while read line; do
    fastq=./fastq/${line}_R1_.final.fastq.gz;
    clipped=./fastq_clipped/${line}_R1_.final.clipped.fastq;
    if [ ! -f $clipped ]; then
        adpt=`python ${bin_libs}/get_adaptors.py ./fastq/$fastq $adapters`;
        zcat $fastq | $fastx -a $adpt -o $clipped;
        gzip $clipped;
    fi;
done < $samples;

if [ ! -d ./star_align ]; then
    mkdir star_align;
fi

# Aligns fastq to genome with STAR
# Removes extraneous STAR output
# Renames STAR sam file
while read line; do
    fastq=./fastq_clipped/${line}_R1_.final.clipped.fastq.gz;
    sam=./star_align/${line}_R1_.final.clipped;    
    if [ ! -f ${sam}.sam ]; then
        $star --runThreadN 24 \
          --genomeDir $genome_dir \
          --readFilesIn ${fastq} \
          --readFilesCommand zcat \
          --outFileNamePrefix ${sam}_;
        rm -r ${sam}_Log.final.out ${sam}_Log.out \
          ${sam}_SJ.out.tab ${sam}_Log.final.out \
          ${sam}_Log.progress.out ${sam}_STARtmp;
        mv ${sam}_Aligned.out.sam ${sam}.sam
    fi;
done < $samples;

# Converts SAM to BAM with picardtools
mkdir ./tmp_dir
tmp_dir="./tmp_dir"
while read line; do
    sam=./star_align/${line}_R1_.final.clipped.sam;
    bam=./star_align/${line}_R1_.final.clipped.bam;
    if [ ! -f $bam ]; then
        $java -jar $picardtools SamFormatConverter \
          I=$sam O=$bam TMP_DIR=$tmp_dir;
    fi;
    #rm $sam
done < $samples;

# filters reads from repeat regions
while read line; do
    bam=./star_align/${line}_R1_.final.clipped.bam;
    repeat_filt=./star_align/${line}_R1_.final.clipped.noRepeats.bam;
    if [ ! -f $repeat_filt ]; then
        $bedtools -abam $bam -b $repeats -wa -v > $repeat_filt;
    fi;
done < $samples;

while read line; do
    #bam=./star_align/${line}_R1_.final.clipped.bam;
    repeat_filt=./star_align/${line}_R1_.final.clipped.noRepeats.bam;
    miR_filt=./star_align/${line}_R1_.final.clipped.noRepeats.noMiR.bam;
    if [ ! -f $miR_filt ]; then
        $bedtools -abam $repeat_filt -b $miRDB -wa -v > $miR_filt;
    fi;
done < $samples;

# filters reads from miR
while read line; do
    bam=./star_align/${line}_R1_.final.clipped.noRepeats.noMiR.bam;
    sorted_bam=./star_align/${line}_R1_.final.clipped.noRepeats.noMiR.sort.bam;
    if [ ! -f $sorted_bam ]; then
        $java -jar $picardtools SortSam \
          I=$bam O=$sorted_bam SORT_ORDER=coordinate \
          TMP_DIR=$tmp_dir;
        $java -jar $picardtools BuildBamIndex \
          I=$sorted_bam;
    fi;
done < $samples;

