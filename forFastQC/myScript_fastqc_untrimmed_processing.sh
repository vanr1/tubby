#!/bin/bash

# (about generating summary for fastqc) 
##################################################################

echo ""
echo "============================================================="
echo "[Starting Bash Script] " ${__file} && date 
echo ""


# making a variable for our main pipeline path
TUBBY_PIPELINE=/data2/han_lab/richardvan/tubby_pipeline

set -e  #exits the script when an error occurs for safety

####################################################################################
# (use samples where lanes have been merged) 
declare -a arrayOfSampleBaseNames=(	
								# "7WT24-R1_S5"	#batch1
								# "8WT24-R2_S2" 		
								# "9WT24-R3_S1" 		
								# "1Tub24-R1_S4" 		
								# "2Tub24-R2_S6" 		
								# "3Tub24-R3_S3" 		

								# "14Nora_S15"	#batch2
								# "15Nora_S9" 		
								# "17Nora_S17" 	
								# "21Nora_S4" 		
								# "23Nora_S2" 	
								# "24Nora_S13" 

								# "10WT3m-R1_S16"				#batch3
								# "11WT3m-R2_S14"	
								"12WT3m-R3_S8"
								# "13Tub24-H1_S17"
								# "18Tub3m-H5_S2"	
								"19WT24-H4_S11"
								# "22WT3m-H1_S10"			
								)

# list all samples we will be looking at, make sure all their paired-end read files exist
echo ""
echo "Listing all samples we will be looking at ..."
for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
do
	if [ -e ${TUBBY_PIPELINE}/data/untrimmed_fastq/${sampleBaseName}_LALL_R1_001.fastq ] && \
	   [ -e ${TUBBY_PIPELINE}/data/untrimmed_fastq/${sampleBaseName}_LALL_R2_001.fastq ]
	then
	    echo "  [CHECKING if files exist for] " ${sampleBaseName} "-- YES"
	else
	    echo "  [CHECKING if files exist for] " ${sampleBaseName} "-- Not Found!"
		echo "Previous file CANNOT BE FOUND, fix it! ..."
	    exit 1
	fi
done
echo "All files exist so can proceed to next step ..."
echo ""
####################################################################################




cd ${TUBBY_PIPELINE}/data/untrimmed_fastq/
mkdir -p ${TUBBY_PIPELINE}/results/fastqc_untrimmed_reads

for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
do

	echo "Running FastQC for 'first pair' of paired-end reads ..."
	fastqc ./${sampleBaseName}_LALL_R1_001.fastq    # makes qc files, html, and zip;  
	echo "Running FastQC for 'second pair' of paired-end reads ..."
	fastqc ./${sampleBaseName}_LALL_R2_001.fastq    # makes qc files, html, and zip;  
        #  don't use multi-thread since not well implemented (slows down val too much)


	echo "Saving FastQC results ..."
	mv -v ./${sampleBaseName}*.zip  ${TUBBY_PIPELINE}/results/fastqc_untrimmed_reads
	mv -v ./${sampleBaseName}*.html ${TUBBY_PIPELINE}/results/fastqc_untrimmed_reads
done




echo "Unzipping files in the results folder ..."
cd ${TUBBY_PIPELINE}/results/fastqc_untrimmed_reads

for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
    do
    unzip ./${sampleBaseName}*.zip
    done

echo "Saving summary to docs (contains all summaries of FastQC per sample) ..."
cat */summary.txt > ${TUBBY_PIPELINE}/docs/fastqc_untrimmed_summaries.txt

echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
