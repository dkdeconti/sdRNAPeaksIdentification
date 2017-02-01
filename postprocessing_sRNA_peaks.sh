# cli args
#samples=$1;
contrasts=$1;

# binary paths
apps_dir="/cccbstore-rc/projects/cccb/apps";
bin_libs="/cccbstore-rc/projects/cccb/pipelines/sRNApipeline/";
fastx="${apps_dir}/fastx_toolkit-0.0.14/bin/fastx_clipper";
star="${apps_dir}/STAR_2.5.0/bin/Linux_x86_64_static/STAR";
bedtools="${apps_dir}/bedtools2-2.22.1/bin/intersectBed";
java="${apps_dir}/jre1.8.0_74/bin/java";
picardtools="${apps_dir}/picard-tools-2.1.1/picard.jar";

# required files paths
db_dir="/cccbstore-rc/projects/db/gatk/hg19";
genome_dir="/cccbstore-rc/projects/db/gatk/hg19/STAR";
genome="/cccbstore-rc/projects/db/gatk/hg19/ucsc.hg19.fasta"

# set up for macs2
#scl enable python27 bash;
#source /ifs/labs/cccb/projects/cccb/apps/MACS2-2.1.1/bin/activate;
#macs="/ifs/labs/cccb/projects/cccb/apps/MACS2-2.1.1/bin/macs2";
py="${bin_libs}/virt_env/bin/python"
filter_py="${bin_libs}/peak_squeeze.py"
homer_apps="/ifs/labs/cccb/projects/cccb/apps/HOMER_v4.8/bin"
PATH=$PATH:${homer_apps}
macs_dir="./macs_output"
#bam_suffix="_R1_.final.clipped.noRepeats.noMiR.sort.bam"

for bam_suffix in "_R1_.final.clipped.noRepeats.noMiR.sort.bam" \
                  "_R1_.final.clipped.sort.bam"; do 

    # sets name addition
    nameAdd="${bam_suffix##*.clipped}"
    nameAdd="${nameAdd%.sort.bam}"    

    # filters for peaks that are narrow
    while read line; do
        experiment=`python -c "s='${line}'; print s.strip('\n').split()[0]"`;
        control=`python -c "s='${line}'; print s.strip('\n').split()[1]"`;
        name="${experiment}-vs-${control}";
        positive="${name}-${nameAdd}-positive_peaks";
        negative="${name}-${nameAdd}-negative_peaks";
        #bam_suffix="_R1_.final.clipped.noRepeats.noMiR.sort.bam"
        bam=./star_align/${experiment}${bam_suffix};
        $py $filter_py $bam ${macs_dir}/${positive}.narrowPeak > \
          ${macs_dir}/${positive}.widthFilt.narrowPeak;
        $py $filter_py $bam ${macs_dir}/${negative}.narrowPeak > \
          ${macs_dir}/${negative}.widthFilt.narrowPeak;
    done < $contrasts;

    ## Annotate peaks with HOMER
    #while read line; do
    #    experiment=`python -c "s='${line}'; print s.strip('\n').split()[0]"`;
    #    control=`python -c "s='${line}'; print s.strip('\n').split()[1]"`;
    #    name="${experiment}-vs-${control}-${nameAdd}";
    #    # All peaks
    #    positive="${macs_dir}/${name}-positive_peaks.narrowPeak";
    #    negative="${macs_dir}/${name}-negative_peaks.narrowPeak";
    #    if [ ! -f $positive ]; then
    #        # convert to HOMER compatible format
    #        $py ${bin_libs}/macs21csv_to_HOMERbed.py $positive > \
    #          ${positive%.narrowPeak}.bed;
    #        $py ${bin_libs}/macs21csv_to_HOMERbed.py $negative > \
    #          ${negative%.narrowPeak}.bed;
    #        # invoke HOMER
    #        ${homer_apps}/annotatePeaks.pl ${positive%.narrowPeak}.bed > \
    #          ${positive%.narrowPeak}.homerAnnot.bed;
    #        ${homer_apps}/annotatePeaks.pl ${negative%.narrowPeak}.bed > \
    #          ${negative%.narrowPeak}.homerAnnot.bed;
    #    fi;
    #    # Width filtered peaks
    #    positive="${macs_dir}/${name}-positive_peaks.widthFilt.narrowPeak";
    #    negative="${macs_dir}/${name}-negative_peaks.widthFilt.narrowPeak";
    #    if [ ! -f $positive ]; then
    #        # convert to HOMER compatible format
    #        $py ${bin_libs}/macs21csv_to_HOMERbed.py $positive > \
    #          ${positive%.narrowPeak}.bed;
    #        $py ${bin_libs}/macs21csv_to_HOMERbed.py $negative > \
    #          ${negative%.narrowPeak}.bed;
    #        # invoke HOMER
    #        ${homer_apps}/annotatePeaks.pl ${positive%.narrowPeak}.bed > \
    #          ${positive%.narrowPeak}.homerAnnot.bed;
    #        ${homer_apps}/annotatePeaks.pl ${negative%.narrowPeak}.bed > \
    #          ${negative%.narrowPeak}.homerAnnot.bed;
    #    fi;
    #done < $contrasts;

    # Add miscellaneous information to all peaks.
    # Add strand information
    # Convert to 
    for f in ${macs_dir}/*.narrowPeak; do
        # No idea why this needed to be put here...
        # virtual env won't load the actual copy here.
        #py="/cccbstore-rc/projects/cccb/projects/2015/9/EH_08312015_882/version_3/bin/virt_env/bin/python"
        peak_wStrand="${f%.narrowPeak}.wStrand.narrowPeak"
        strand=${f%_peaks.*}
        strand=${strand##*-}
        bam_suffix="_R1_.final.clipped.sort.bam"
        if [ $strand == "positive" ]; then
            bam=${f%-vs-*};
            bam=${bam%-};
            bam=./star_align/${bam##*/}${bam_suffix};
        else
            bam=${f##*-vs-};
            bam=${bam%%-*};
            bam=./star_align/${bam%-negative*}${bam_suffix};
        fi;
        bed="${f%.narrowPeak}.bed"
        annot="${f%.narrowPeak}.homerAnnot.bed"
        annotWStrand="${annot%.bed}.wStrandNum.bed"
        # add strand information to peaks
        $py ${bin_libs}/add_strand_to_peaks.py $f $bam > $peak_wStrand;
        # convert macs output to Homer capable bed
        $py ${bin_libs}/macs21csv_to_HOMERbed.py $peak_wStrand > $bed;
        # annotate peaks with homer
        ${homer_apps}/annotatePeaks.pl $bed hg19 > $annot;
        # Add strand data to homer annotations
        echo "$py ${bin_libs}/add_strand_to_HOMER.py $annot $bam > $annotWStrand;"
        $py ${bin_libs}/add_strand_to_HOMER.py $annot $bam > $annotWStrand;
    done;
done;

