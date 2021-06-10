

#!/bin/bash

##################################################################

echo ""
echo "============================================================="
echo "[Starting Bash Script] " ${__file} && date 
echo ""


declare -a arrayOfSampleBaseNames=(	"9WT24-R3_S1" 
								"8WT24-R2_S2" 
								"3Tub24-R3_S3" 
								"1Tub24-R1_S4" 
								"7WT24-R1_S5" 
								"2Tub24-R2_S6" 
								)

for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
do
	echo "  [START working on sample] " ${sampleBaseName} && date 


	## ( map each lane to see if any lane effects exist)
	STAR --runThreadN 20 --genomeDir ./genomeDirectory --readFilesIn ./raw_fastq/${sampleBaseName}_L001_R1_001.fastq ./raw_fastq/${sampleBaseName}_L001_R2_001.fastq --sjdbGTFfile ./mouse_annotation/Mus_musculus.GRCm38.101.gtf --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate

	cp -av Log.final.out starMapJob_${sampleBaseName}_L001.log
	cp -av Aligned.sortedByCoord.out.bam tmp/bamFile_${sampleBaseName}_L001.bam
	# cat starMapJob_R1_L001.log

	STAR --runThreadN 20 --genomeDir ./genomeDirectory --readFilesIn ./raw_fastq/${sampleBaseName}_L002_R1_001.fastq ./raw_fastq/${sampleBaseName}_L002_R2_001.fastq --sjdbGTFfile ./mouse_annotation/Mus_musculus.GRCm38.101.gtf --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate
	cp -av Log.final.out starMapJob_${sampleBaseName}_L002.log
	cp -av Aligned.sortedByCoord.out.bam tmp/bamFile_${sampleBaseName}_L002.bam
	# cat starMapJob_R1_L002.log

	STAR --runThreadN 20 --genomeDir ./genomeDirectory --readFilesIn ./raw_fastq/${sampleBaseName}_L003_R1_001.fastq ./raw_fastq/${sampleBaseName}_L003_R2_001.fastq --sjdbGTFfile ./mouse_annotation/Mus_musculus.GRCm38.101.gtf --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate
	cp -av Log.final.out starMapJob_${sampleBaseName}_L003.log
	cp -av Aligned.sortedByCoord.out.bam tmp/bamFile_${sampleBaseName}_L003.bam
	# cat starMapJob_R1_L003.log

	STAR --runThreadN 20 --genomeDir ./genomeDirectory --readFilesIn ./raw_fastq/${sampleBaseName}_L004_R1_001.fastq ./raw_fastq/${sampleBaseName}_L004_R2_001.fastq --sjdbGTFfile ./mouse_annotation/Mus_musculus.GRCm38.101.gtf --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate
	cp -av Log.final.out starMapJob_${sampleBaseName}_L004.log
	cp -av Aligned.sortedByCoord.out.bam tmp/bamFile_${sampleBaseName}_L004.bam
	# cat starMapJob_R1_L004.log



	## (  then merge with merge_bam
	samtools merge tmp/bamFile_${sampleBaseName}_LALL.bam tmp/bamFile_${sampleBaseName}_L001.bam tmp/bamFile_${sampleBaseName}_L002.bam tmp/bamFile_${sampleBaseName}_L003.bam tmp/bamFile_${sampleBaseName}_L004.bam

	samtools sort -@ 8 -o tmp/bamFile_${sampleBaseName}_LALL.sorted.bam tmp/bamFile_${sampleBaseName}_LALL.bam 

	htseq-count tmp/bamFile_${sampleBaseName}_LALL.sorted.bam ./mouse_annotation/Mus_musculus.GRCm38.101.gtf --format=bam --order=pos > ${sampleBaseName}_counts.txt

	# sort -k2,2gr 7WT24-R1_S5_counts.txt | head -20


	echo "  [DONE working on sample] " ${sampleBaseName} && date 
done

## _TODO_ - make generic loop to do
	# generate PCGcounts
	# convert to FPKM

# merge fpkm files



echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
