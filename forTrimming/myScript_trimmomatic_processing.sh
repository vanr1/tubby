#!/bin/bash

# (about trimming) no adapter; just regular trimming 
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
                        # "7WT24-R1_S5"   #batch1
                        # "8WT24-R2_S2"      
                        # "9WT24-R3_S1"      
                        # "1Tub24-R1_S4"     
                        # "2Tub24-R2_S6"     
                        # "3Tub24-R3_S3"     

                        # "14Nora_S15" #batch2
                        # "15Nora_S9"     
                        # "17Nora_S17"    
                        # "21Nora_S4"     
                        # "23Nora_S2"  
                        # "24Nora_S13" 

                        "10WT3m-R1_S16"          #batch3
                        "11WT3m-R2_S14" 
                        "12WT3m-R3_S8"
                        "13Tub24-H1_S17"
                        "18Tub3m-H5_S2" 
                        "19WT24-H4_S11"
                        "22WT3m-H1_S10"       
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

echo "Running Trimmomatic ..."

for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
do
   echo "  [working with sample] $sampleBaseName"
   java -jar ${TUBBY_PIPELINE}/bin/trimmomatic-0.39.jar PE \
               -threads 40 \
               ${sampleBaseName}_LALL_R1_001.fastq ${sampleBaseName}_LALL_R2_001.fastq \
               ${sampleBaseName}_LALL_1.trim.fastq ${sampleBaseName}_LALL_1un.trim.fastq \
               ${sampleBaseName}_LALL_2.trim.fastq ${sampleBaseName}_LALL_2un.trim.fastq \
               SLIDINGWINDOW:4:20 MINLEN:25 \
               2>&1 | tee ${sampleBaseName}_trimmomatic_processing.log
   echo "  [done with sample] $sampleBaseName"
   echo "  "
done

mkdir -p ${TUBBY_PIPELINE}/results/logs/from_trimmomatic

echo "Moving output logs that tells (percent of reads discarded per sample, percent of reads we keep in both pairs) ..."
mv -v *.log ${TUBBY_PIPELINE}/results/logs/from_trimmomatic

mkdir -p ${TUBBY_PIPELINE}/data/trimmed_fastq/
mv -v ${TUBBY_PIPELINE}/data/untrimmed_fastq/*.trim* ${TUBBY_PIPELINE}/data/trimmed_fastq/


echo "Saving summary to docs (contains grep of the relevant lines) ..."
cd ${TUBBY_PIPELINE}/results/logs/from_trimmomatic
grep Surviving *processing.log > ${TUBBY_PIPELINE}/docs/trimmomatic_summaries.txt


echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
