#!/bin/bash

# (about finding the tubby mutation) 
   # (literature) Kleyn, 1996 - Identification and Characterization of the Mouse Obesity Gene tubby: A Member of a Novel Gene Family
   # (rest of nucleotides) search for "ACGGCAATGACC" at this link
   			# https://www.ncbi.nlm.nih.gov/nuccore/NC_000073.7?report=fasta&from=108610087&to=108633666

##################################################################

echo ""
echo "============================================================="
echo "[Starting Bash Script] " ${__file} && date 
echo ""

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


mkdir -p ${TUBBY_PIPELINE}/results/mutscan/

# generate results of mutscan on all samples
echo "[RV] generate results of mutscan for each sample "

for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
do
	# use string split to grab only the sample name
	#	(ex)	IN="bla@some.com;john@home.com"
	#			arrIN=(${IN//;/ })
	#			echo ${arrIN[1]}                  # Output: john@home.com

	echo "  [calling ./mutscan on the sample] " ${sampleBaseName}



	# (proof of concept command)
	./bin/mutscan 	-1 ${TUBBY_PIPELINE}/data/untrimmed_fastq/${sampleBaseName}_LALL_R1_001.fastq \
					-2 ${TUBBY_PIPELINE}/data/untrimmed_fastq/${sampleBaseName}_LALL_R2_001.fastq \
					-m ${TUBBY_PIPELINE}/bin/mutation_tubby.csv \
					-t 40 \
					-h ${TUBBY_PIPELINE}/results/mutscan/${sampleBaseName}_mutscan_htmlReport.html \
					--standalone \
					>  ${TUBBY_PIPELINE}/results/mutscan/${sampleBaseName}_mutscan_readsFound.txt 

done


#### [post command] to generate useful info for docs
# [richardvan@atty tubby_pipeline]$ for filename in results/mutscan/*.txt
# > do
# >   grep ^mutscan $filename
# >   grep ^- $filename | wc -l
# >   echo " "
# > done




echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""