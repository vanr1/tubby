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

mkdir -p ${TUBBY_PIPELINE}/results/bcf/
mkdir -p ${TUBBY_PIPELINE}/results/vcf/

echo "Running bcftools ..."


for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
do
   echo "  [STARTING on sample] " ${sampleBaseName} && date 

   # making a bunch of variables for files we want to use
   genome=${TUBBY_PIPELINE}/data/ref_genome/mouse_ensemble_GRCm38p6.fasta
   sorted_bam=${TUBBY_PIPELINE}/results/bam/${sampleBaseName}.aligned.sorted.bam
   raw_bcf=${TUBBY_PIPELINE}/results/bcf/${sampleBaseName}_raw.bcf
   variants=${TUBBY_PIPELINE}/results/bcf/${sampleBaseName}_variants.vcf
   final_variants=${TUBBY_PIPELINE}/results/vcf/${sampleBaseName}_final_variants.vcf 

   # the bioinformatics commands themselves here
   echo "    calculating read coverage of positions in the genome..."  
   bcftools mpileup --threads 80 -O b -o $raw_bcf -f $genome $sorted_bam
   echo "    detecting the SNPs..."  
   bcftools call --threads 80 --ploidy 1 -m -v -o $variants $raw_bcf 
   echo "    filter and report the SNP variants in variant calling format (VCF)"
   vcfutils.pl varFilter $variants > $final_variants


   echo "    -- consider what information can get from the VCF file in results/vcf/"



   echo "  [DONE working on sample] " ${sampleBaseName} 
   echo "   "
done



echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
