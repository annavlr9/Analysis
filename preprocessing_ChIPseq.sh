
# working dir
wdir="/mnt/beegfs/alopez/Gemma_project/ChIPseq/data_sept_2023"
cd $wdir


echo "############################################################## Step 1: Run fastqc in raw samples ######################################################################"
	module load fastqc-0.11.9-gcc-11.2.0-dd2vd2m
	fastqc *
echo "############################################################## Step 1: DONE ######################################################################"




#echo "############################################################## Step 2: Trim reads with Trim Galore (if necessary) ######################################################################"
	module load TrimGalore/0.6.6
	trim_galore --paired --length 35 --fastqc $read1 $read2

	for i in /mnt/beegfs/alopez/Gemma_project/ChIPseq/data_sept_2023/*.fq.gz
	do
        	trim_galore --length 35 --fastqc "${i}"	
	done
#echo "############################################################## Step 2: DONE ######################################################################"




echo "############################################################## Step3: Alignment with bowtie2 ######################################################################"
	index_human="/mnt/beegfs/public/references/index/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as"

	module load Bowtie2/2.4.4.1-GCC-11.2.0

	for i in /mnt/beegfs/alopez/Gemma_project/ChIPseq/data_sept_2023/*.fq.gz
	do
        	bowtie2 -p 10 -x $index_human -U "${i}" -S "${i%.fq.gz}".bt2.sam 2> "${i%.fq.gz}".bt2.log
	done
echo "############################################################## Step 3: DONE ######################################################################"




echo "############################################################## Step 4: Samtools filtering  ######################################################################"
	module load SAMtools/1.13-foss-2021b

	echo -e "### Samtools SAM to BAM ###"
	for i in *.bt2.sam; do sample=`echo $i | awk -F "." '{print $1}'`; samtools view -bS $i > ${sample}.unsorted.bam; done


	echo -e "### Samtools filter BAM ###"
	for i in *.unsorted.bam; do sample=`echo $i | awk -F "." '{print $1}'`; samtools view -F 2304 -b -q 10 $i > ${sample}.clean.bam; done

	for i in *.clean.bam; do sample=`echo $i | awk -F "." '{print $1}'`; samtools sort $i -o ${sample}.sorted.bam; done

	# Create indeces for all the bam files for visualization and QC
	#samtools index ${sample}.sorted.bam
	for i in *.sorted.bam; do sample=`echo $i | awk -F "." '{print $1}'`; samtools index $i; done
echo "############################################################## Step 4: DONE ######################################################################"




echo "############################################################## Step 5: Picard MarkDuplicates  ######################################################################"
	module load picard/2.26.3-Java-11

	for i in *.sorted.bam; do sample=`echo $i | awk -F "." '{print $1}'`; java -jar $EBROOTPICARD/picard.jar MarkDuplicates -INPUT $i -OUTPUT ${sample}.rmdup.bam -METRICS_FILE ${sample_name}.stats -REMOVE_DUPLICATES true ASSUME_SORTED true -VERBOSITY WARNING ; done
echo "############################################################## Step 5: DONE  ######################################################################"




echo "############################################################## Step 6: Remove blacklisted regions started ######################################################################"
	module load BEDTools/2.30.0-GCC-11.2.0

 	for i in *.rmdup.bam ; do sample=`echo $i | awk -F "." '{print $1}'`; bedtools intersect -nonamecheck -v -abam ${sample}.rmdup.bam -b  ${BLACKLIST} > ${sample}.filtered.bam


	BLACKLIST="/mnt/beegfs/alopez/bin/ENCFF356LFX.bed"
	bedtools intersect -v -a ${sample}.peaks.narrowPeak -b ${BLACKLIST} > ${sample}.peaks.narrowPeak.filt
echo "############################################################## Step 6: DONE ######################################################################"





echo "############################################################## Step 7: MACS2 Peak Calling  ######################################################################"
	module load MACS2/2.2.5-foss-2021b-Python-3.8.5	
	for i in *.filtered.bam; do sample=`echo $i | awk -F "." '{print $1}'`; macs2 callpeak -t $sample.filtered.bam -n $sample.filtered.peaks -g hs --keep-dup auto --call-summits --verbose 3 ; done
echo "############################################################## Step 7: DONE ######################################################################"







echo "############################################################## Step 8: Get consensus of the replicates ######################################################################"
	#cat *.narrowPeak.filt | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > merged_replicates_peaks.168h.bed

	# Add a column to know to which timepoint each peak corresponds to
	#awk 'BEGIN{ FS = OFS = "\t" } { print $0, (NR==1? "Timepoint" : "168h") }' merged_replicates_peaks.168h.bed > merged_replicates_peaks.168h_timepoint.bed
echo "############################################################## Step 8: DONE ######################################################################"




echo "############################################################## Step 9: Run multiQC ######################################################################"
	module load multiqc/1.12
	#python -m multiqc .
echo "############################################################## Step 9: DONE ######################################################################"





echo "############################################################## Step 10: OBTAIN BIGWIG FILES ######################################################################"
	module load deepTools/3.5.1-foss-2021b

	human_effGenomeSize=2913022398 #from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
	#mouse_effGenomeSize=2652783500
	PPN=10

	for i in *.rmdup.bam; do sample=`echo $i | awk -F "." '{print $1}'`; bamCoverage  --numberOfProcessors $PPN --binSize 10 --smoothLength 50 --normalizeUsing RPKM --effectiveGenomeSize $human_effGenomeSize --bam $i -o 	$i.sorted.RMDUP.RPKM.smooth50.bw ; done
echo "############################################################## Step 10: DONE ######################################################################"






