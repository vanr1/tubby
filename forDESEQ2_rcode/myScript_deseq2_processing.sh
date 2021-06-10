#!/bin/bash

# (about preprocessing) remember to run this in ATTY as only it can handle Rscript calls
##################################################################

echo ""
echo "============================================================="
echo "[Starting Bash Script] " ${__file} && date 
echo ""

# making a variable for our main pipeline path
TUBBY_PIPELINE=/data2/han_lab/richardvan/tubby_pipeline

set -e  #exits the script when an error occurs for safety


                           # declare -a arrayOfSampleBaseNames=( 

                           #                         "1Tub24-R1_S4"       # first batch (oct 2020)
                           #                         "2Tub24-R2_S6"       
                           #                         "3Tub24-R3_S3"    
                           #                         "7WT24-R1_S5"
                           #                         "8WT24-R2_S2"     
                           #                         "9WT24-R3_S1"        

                           #                         "14Nora_S15"      # second batch (dec 2020)
                           #                         "15Nora_S9"       
                           #                         "17Nora_S17"   
                           #                         "21Nora_S4"       
                           #                         "23Nora_S2"    
                           #                         "24Nora_S13" 

                           #                                        # third batch (october 2020)
                           #                         )

echo "Renaming the raw htseq outputs into a new folder with better naming ..."

## need to create the input_sets accordingly with what sequences are available and renaming to make more sense
mkdir -p ${TUBBY_PIPELINE}/data/deseq2_input_sets/
mkdir -p ${TUBBY_PIPELINE}/data/renamed_counts

# here copy over the files to their respective folders; so can use different input sets
   # tub 24w 
cp -v ${TUBBY_PIPELINE}/results/counts/1Tub24-R1_S4.htseq.counts \
 ${TUBBY_PIPELINE}/data/renamed_counts/tub_24w_r1.txt
cp -v ${TUBBY_PIPELINE}/results/counts/2Tub24-R2_S6.htseq.counts \
 ${TUBBY_PIPELINE}/data/renamed_counts/tub_24w_r2.txt
cp -v ${TUBBY_PIPELINE}/results/counts/3Tub24-R3_S3.htseq.counts \
 ${TUBBY_PIPELINE}/data/renamed_counts/tub_24w_r3.txt
   # wt 24w 
cp -v ${TUBBY_PIPELINE}/results/counts/7WT24-R1_S5.htseq.counts \
 ${TUBBY_PIPELINE}/data/renamed_counts/wt_24w_r1.txt
cp -v ${TUBBY_PIPELINE}/results/counts/8WT24-R2_S2.htseq.counts \
 ${TUBBY_PIPELINE}/data/renamed_counts/wt_24w_r2.txt
cp -v ${TUBBY_PIPELINE}/results/counts/9WT24-R3_S1.htseq.counts \
 ${TUBBY_PIPELINE}/data/renamed_counts/wt_24w_r3.txt
   # _toOrganizeLater_onceHaveAllSamples_HereComingInRandomly
# cp -v ${TUBBY_PIPELINE}/results/counts/14Nora_S15.htseq.counts \
#  ${TUBBY_PIPELINE}/data/renamed_counts/tub_24w_h2.txt
# cp -v ${TUBBY_PIPELINE}/results/counts/15Nora_S9.htseq.counts \
#  ${TUBBY_PIPELINE}/data/renamed_counts/tub_24w_h3.txt
# cp -v ${TUBBY_PIPELINE}/results/counts/17Nora_S17.htseq.counts \
#  ${TUBBY_PIPELINE}/data/renamed_counts/tub_3m_h2.txt
# cp -v ${TUBBY_PIPELINE}/results/counts/21Nora_S4.htseq.counts \
#  ${TUBBY_PIPELINE}/data/renamed_counts/wt_24w_h3.txt
# cp -v ${TUBBY_PIPELINE}/results/counts/23Nora_S2.htseq.counts \
#  ${TUBBY_PIPELINE}/data/renamed_counts/wt_3m_h2.txt
# cp -v ${TUBBY_PIPELINE}/results/counts/24Nora_S13.htseq.counts \
#  ${TUBBY_PIPELINE}/data/renamed_counts/wt_3m_h3.txt

echo "Creating deseq2_input_sets ..."

# first one just includes everything from the renamed_counts folder; need to label untreated vs treated
#  (untreated = wt ; treated = tub)
mkdir -p ${TUBBY_PIPELINE}/data/deseq2_input_sets/wt_vs_tub_all
rm -f ${TUBBY_PIPELINE}/data/deseq2_input_sets/wt_vs_tub_all/*
cp -v  ${TUBBY_PIPELINE}/data/renamed_counts/* \
 ${TUBBY_PIPELINE}/data/deseq2_input_sets/wt_vs_tub_all
cd ${TUBBY_PIPELINE}/data/deseq2_input_sets/wt_vs_tub_all
for filename in *wt*.txt
do
   name=$(basename $filename .txt)
   mv $filename ${name}_untreated.txt
done
for filename in *tub*.txt
do
   name=$(basename $filename .txt)
   mv $filename ${name}_treated.txt
done



echo ""
echo "Checking deseq2_input_sets ..."
num_files=$(ls ${TUBBY_PIPELINE}/data/deseq2_input_sets/wt_vs_tub_all/* | wc -w )
echo "  --  wt_vs_tub_all.txt has the following " $num_files " total files:"  
ls ${TUBBY_PIPELINE}/data/deseq2_input_sets/wt_vs_tub_all
echo ""


echo "Running myScript_deseq2_processing.r ..."


Rscript ${TUBBY_PIPELINE}/scripts/myScript_deseq2___commands.r

# for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
# do
#    echo "  [STARTING on input_set (out of six analysis total] " ${sampleBaseName} && date 


#    echo "  [DONE working on input_set] " ${sampleBaseName} 
#    echo "   "
# done




echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
