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

                        "10WT3m-R1_S16"            #batch3
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
   if [ -e ${TUBBY_PIPELINE}/data/trimmed_fastq/${sampleBaseName}_LALL_1.trim.fastq ] && \
      [ -e ${TUBBY_PIPELINE}/data/trimmed_fastq/${sampleBaseName}_LALL_2.trim.fastq ]
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

mkdir -p ${TUBBY_PIPELINE}/results/bam/

echo "Running Samtools ..."


for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
do
   echo "  [STARTING on sample] " ${sampleBaseName} && date 

   # making a bunch of variables for files we want to use
   sam=${TUBBY_PIPELINE}/results/sam/${sampleBaseName}.aligned.sam
   bam=${TUBBY_PIPELINE}/results/bam/${sampleBaseName}.aligned.bam
   sorted_bam=${TUBBY_PIPELINE}/results/bam/${sampleBaseName}.aligned.sorted.bam
   flagstat_bam=${TUBBY_PIPELINE}/docs/bam_flagstat_${sampleBaseName}.txt

   # the bioinformatics commands themselves here
   echo "    converting..."  
   samtools view --threads 40 -S -b $sam > $bam      #wierd synthax, but this is how sam->bam conversion is done
   echo "    sorting..."  
   samtools sort --threads 40 -o $sorted_bam $bam    #the sorted bam file will be smaller than the unsorted one bc compression is better on sorted
   samtools index $sorted_bam   
   echo "    generating flagstat to docs..."  
   samtools flagstat $sorted_bam > $flagstat_bam

   #untested
         #    Fastest way to count number of reads

         # From http://left.subtree.org/2012/04/13/counting-the-number-of-reads-in-a-bam-file/#comment-403

         # count
         # #number of reads
         # samtools idxstats in.bam | awk '{s+=$3+$4} END {print s}'
         # #number of mapped reads
         # samtools idxstats in.bam | awk '{s+=$3} END {print s}'

   echo "  [DONE working on sample] " ${sampleBaseName} 
   echo "   "
done

echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
