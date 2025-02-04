sample=$1
cd $sample

read1=$(ls *_1.fq.gz)
read2=$(ls *_2.fq.gz)

echo "read1=$read1"
echo "read2=$read2"


#module load STAR/2.7.8a-GCC-11.2.0
#STAR --runMode genomeGenerate --genomeDir /mnt/beegfs/alopez/bin/annotation/GRCh38 --genomeFastaFiles /mnt/beegfs/alopez/bin/annotation/GRCh38.primary_assembly.genome.fa --sjdbGTFfile /mnt/beegfs/alopez/bin/annotation/gencode.v43.annotation.gtf --sjdbOverhang 50 --outFileNamePrefix GRCh38.
#echo -e "index done"



echo -e "##################################################### 1) QC with FastQC started #####################################################"
module load FastQC/0.11.9-Java-11
fastqc *
echo " ##################################################### 1) QC with FastQC finished #####################################################"



echo -e " ##################################################### 3) STAR alignment started #####################################################"
module load STAR/2.7.8a-GCC-11.2.0

STAR --genomeDir /mnt/beegfs/alopez/bin/annotation/GRCh38 --readFilesIn $read1 $read2 --runThreadN 8 --outFileNamePrefix $sample. --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
echo -e " ##################################################### 3) STAR alignment done #####################################################"



echo -e " ##################################################### 4) samtools index started #####################################################"
module load SAMtools/1.19.2-foss-2021b

samtools index $sample.Aligned.sortedByCoord.out.bam
samtools index $sample.Aligned.sortedByCoord.out.bam
echo -e " ##################################################### 4) samtools index done #####################################################"





#echo -e " ##################################################### 7) QUANTIFY READS WITH FEATURECOUNTS #####################################################"
#module load Subread/2.0.3-GCC-11.2.0
#featureCounts -p --countReadPairs -a /mnt/beegfs/alopez/bin/annotation/gencode.v43.annotation.gtf -o a2.4_rawCounts.txt *Aligned.sortedByCoord.out.bam
#echo -e " ##################################################### 7) featureCounts done #####################################################"

