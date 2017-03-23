# cli args
#samples=$1;
contrasts=$1;

# binary paths
apps_dir="/cccbstore-rc/projects/cccb/apps";
bin_libs="/cccbstore-rc/projects/cccb/pipelines/sRNApipeline";
fastx="${apps_dir}/fastx_toolkit-0.0.14/bin/fastx_clipper";
star="${apps_dir}/STAR_2.5.0/bin/Linux_x86_64_static/STAR";
bedtools="${apps_dir}/bedtools2-2.22.1/bin/intersectBed";
java="${apps_dir}/jre1.8.0_74/bin/java";
picardtools="${apps_dir}/picard-tools-2.1.1/picard.jar";

# required files paths
db_dir="/cccbstore-rc/projects/db/gatk/hg19";
genome_dir="/cccbstore-rc/projects/db/gatk/hg19/STAR";

# set up for macs2
#scl enable python27 bash;
#source /ifs/labs/cccb/projects/cccb/apps/MACS2-2.1.1/bin/activate;
py="/ifs/labs/cccb/projects/cccb/apps/MACS2-2.1.1/bin/python"
macs="/ifs/labs/cccb/projects/cccb/apps/MACS2-2.1.1/bin/macs2";
macs_dir="./macs_output"
if [ ! -d $macs_dir ]; then
    mkdir $macs_dir;
fi;

for bam_suffix in "_R1_.final.clipped.noRepeats.noMiR.sort.bam" \
                  "_R1_.final.clipped.sort.bam"; do
    nameAdd=${bam_suffix##*.clipped}
    nameAdd=${nameAdd%.sort.bam}
    # reads contrast file to run macs2
    while read line; do
        experiment=`python -c "s='${line}'; print s.strip('\n').split()[0]"`;
        control=`python -c "s='${line}'; print s.strip('\n').split()[1]"`;
        #bam_suffix="_R1_.final.clipped.noRepeats.noMiR.sort.bam";
        exp_bam="./star_align/${experiment}${bam_suffix}";
        ctrl_bam="./star_align/${control}${bam_suffix}";
        $py $macs callpeak \
	        -t $exp_bam \
	        -c $ctrl_bam \
            -f BAM \
            --outdir $macs_dir \
	        --name ${experiment}-vs-${control}-${nameAdd}-positive \
            -g hs \
	        --keep-dup all \
	        --nomodel;
        $py $macs callpeak \
	        -t $exp_bam \
	        -c $ctrl_bam \
	        -f BAM \
            --outdir $macs_dir \
	        --name ${experiment}-vs-${control}-${nameAdd}-negative \
	        -g hs \
	        --keep-dup all \
	        --nomodel;
    done < $contrasts;
done;

