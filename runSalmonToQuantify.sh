

#!/bin/bash

# (about)  --- use salmon for quasialignment    (salmon gives) TPM = normalized within each sample
#   overall - for DESEQ2; using salmon pipeline, merge each sample by lanes so end up with 6 total samples....then have 6 salmon calls, merge outputs into a matrix where you have ~50,000 rows of geneID by 6 cols of sampleID.
#   goal: get normalization factor and generate two visualizations  for samples
#   -  PCA plot = see clustering of absolute values
#   -  MDS (multi-dimensional scaling) plot = based on distance b/w samples





##################################################################

echo ""
echo "============================================================="
echo "[Starting Bash Script] " ${__file} && date 
echo ""

# prepare files for index creation (run this inside mouse_annotation folder)
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz

# grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
# sed -i.bak -e 's/>//g' decoys.txt
# cat gencode.vM25.transcripts.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome.fa.gz

# make index once
salmon index --threads 40 \
			 -t mouse_annotation/gentrome.fa.gz \
			 -d mouse_annotation/decoys.txt \
			 -i transcripts_index \
			 --gencode

# (use samples where lanes have been merged)
declare -a arrayOfSampleBaseNames=(	
								"7WT24-R1_S5"
								"8WT24-R2_S2" 		# _TODO__uncomment in samples 1 through 6 when pipeline ready for 6 sample calls total
								"9WT24-R3_S1" 		# _TODO__uncomment in samples 1 through 6 when pipeline ready for 6 sample calls total
								"1Tub24-R1_S4" 		# _TODO__uncomment in samples 1 through 6 when pipeline ready for 6 sample calls total
								"2Tub24-R2_S6" 		# _TODO__uncomment in samples 1 through 6 when pipeline ready for 6 sample calls total
								"3Tub24-R3_S3" 		# _TODO__uncomment in samples 1 through 6 when pipeline ready for 6 sample calls total
								)



for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
do
	echo "  [START working on sample using SALMON] " ${sampleBaseName} && date 

	salmon quant --threads 40 \
				 -i transcripts_index \
				 --libType A \
				 	-1 raw_fastq/${sampleBaseName}_LALL_R1_001.fastq \
				 	-2 raw_fastq/${sampleBaseName}_LALL_R2_001.fastq \
				 --gcBias \
				 --validateMappings -o transcripts_quant_${sampleBaseName}
done




echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
