#!/bin/bash

# (about htseq) generate counts with a reliable library
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

mkdir -p ${TUBBY_PIPELINE}/results/counts/

echo "Running htseq-count ..."


for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
do
   echo "  [STARTING on sample] " ${sampleBaseName} && date 

   # making a bunch of variables for files we want to use
   gff=${TUBBY_PIPELINE}/data/ref_genome/mouse_ensemble_GRCm38p6.gff3
   sorted_bam=${TUBBY_PIPELINE}/results/bam/${sampleBaseName}.aligned.sorted.bam

   # the bioinformatics commands themselves here
   echo "    generating counts based samples sorted by position ..."  
   
                                                            #notice using .gff3 file so it is "Name" not "gene_id"(default); 
                                                            #notice --type needs to be 'gene' for instead of the default of 'exon'
   htseq-count --format bam --order pos --type gene --idattr Name \
                  $sorted_bam $gff > ${TUBBY_PIPELINE}/results/counts/${sampleBaseName}.htseq.counts




   echo "  [DONE working on sample] " ${sampleBaseName} 
   echo "   "
done


   # sort -k2,2gr 7WT24-R1_S5_counts.txt | head -20



# echo "Moving output logs that tells (percent of reads discarded per sample, percent of reads we keep in both pairs) ..."
# mv -v *.log ${TUBBY_PIPELINE}/results/trimmomatic_processing_logs/

# mkdir -p ${TUBBY_PIPELINE}/data/trimmed_fastq/
# mv -v ${TUBBY_PIPELINE}/data/untrimmed_fastq/*.trim* ${TUBBY_PIPELINE}/data/trimmed_fastq/


# echo "Saving summary to docs (contains grep of the relevant lines) ..."
# cd ${TUBBY_PIPELINE}/results/trimmomatic_processing_logs/
# grep Surviving *processing.log > ${TUBBY_PIPELINE}/docs/trimmomatic_summaries.txt


echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
