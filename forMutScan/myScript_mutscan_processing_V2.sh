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
TUBBY_PIPELINE_RENAMED=/data2/han_lab/richardvan/tubby_pipeline/data/renamed_fastq
TUBBY_PIPELINE_TMP=/data2/han_lab/richardvan/tubby_pipeline/tmp
set -e  #exits the script when an error occurs for safety

####################################################################################
# (use samples where lanes have been merged) 

samples=(sample01 sample02 sample03 sample04 sample05 sample06 \
         sample07 sample08 sample09 sample10 sample11 sample12 \
         sample13 sample14 sample15 sample16 sample17 sample18 \
         sample19 sample20 sample21 sample22 sample23 sample24 \
        )

declare -A renamed_fastq=(
    [sample01]="tb_24d_ret1"
    [sample02]="tb_24d_ret2"
    [sample03]="tb_24d_ret3"
    [sample04]="tb_3mo_ret1"
    [sample05]="tb_3mo_ret2"
    [sample06]="tb_3mo_ret3"
    [sample07]="wt_24d_ret1"
    [sample08]="wt_24d_ret2"
    [sample09]="wt_24d_ret3"
    [sample10]="wt_3mo_ret1"
    [sample11]="wt_3mo_ret2"
    [sample12]="wt_3mo_ret3"
    [sample13]="tb_24d_hyp1"
    [sample14]="tb_24d_hyp2"
    [sample15]="tb_24d_hyp3"
    [sample16]="tb_3mo_hyp1"
    [sample17]="tb_3mo_hyp2"
    [sample18]="tb_3mo_hyp3"
    [sample19]="wt_24d_hyp1"
    [sample20]="wt_24d_hyp2"
    [sample21]="wt_24d_hyp3"
    [sample22]="wt_3mo_hyp1"
    [sample23]="wt_3mo_hyp2"
    [sample24]="wt_3mo_hyp3"
    )


# for sample in ${samples[@]}
# do
#     # echo $sample raw is ${raw_fastq[$sample]} 
#     # echo $sample renamed is ${renamed_fastq[$sample]} 
#     cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L001_R1_001.fastq.gz \
#            ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L1_R1.fastq.gz 

#     cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L002_R1_001.fastq.gz \
#            ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L2_R1.fastq.gz 

#     cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L003_R1_001.fastq.gz \
#            ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L3_R1.fastq.gz 

#     cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L004_R1_001.fastq.gz \
#            ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L4_R1.fastq.gz  

#     cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L001_R2_001.fastq.gz \
#            ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L1_R2.fastq.gz 

#     cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L002_R2_001.fastq.gz \
#            ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L2_R2.fastq.gz 

#     cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L003_R2_001.fastq.gz \
#            ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L3_R2.fastq.gz 

#     cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L004_R2_001.fastq.gz \
#            ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L4_R2.fastq.gz 
# done


# list all samples we will be looking at, make sure all their paired-end read files exist
echo ""
echo "Listing all samples we will be looking at ..."
for sample in ${samples[@]}
do
	if [ -e ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L1_R1.fastq.gz ] && \
	   [ -e ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L1_R2.fastq.gz ]
	then
	    echo "  [CHECKING if files exist for] " ${renamed_fastq[$sample]} "-- YES"
	else
	    echo "  [CHECKING if files exist for] " ${renamed_fastq[$sample]} "-- Not Found!"
		echo "Previous file CANNOT BE FOUND, fix it! ..."
	    exit 1
	fi
done


echo "All files exist so can proceed to next step ..."
echo ""
####################################################################################


mkdir -p ${TUBBY_PIPELINE}/results_mutscan/

# generate results of mutscan on all samples
echo "[RV] generate results of mutscan for each sample "

for sample in ${samples[@]}
do
	# use string split to grab only the sample name
	#	(ex)	IN="bla@some.com;john@home.com"
	#			arrIN=(${IN//;/ })
	#			echo ${arrIN[1]}                  # Output: john@home.com

	echo "  [combining lanes 1, 2, 3, and 4 for R1] " ${sample}
	cat ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L1_R1.fastq.gz \
		  ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L2_R1.fastq.gz \
		  ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L3_R1.fastq.gz \
		  ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L4_R1.fastq.gz \
		  > ${TUBBY_PIPELINE_TMP}/mutscan_combined_lanes_R1.fastq.gz
	echo "  [combining lanes 1, 2, 3, and 4 for R2] " ${sample}
	cat ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L1_R2.fastq.gz \
		  ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L2_R2.fastq.gz \
		  ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L3_R2.fastq.gz \
		  ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L4_R2.fastq.gz \
		  > ${TUBBY_PIPELINE_TMP}/mutscan_combined_lanes_R2.fastq.gz

	echo "  [calling ./mutscan on the sample] " ${sample} "-- all lanes -- "
	./bin/mutscan 	-1 ${TUBBY_PIPELINE_TMP}/mutscan_combined_lanes_R1.fastq.gz \
					-2 ${TUBBY_PIPELINE_TMP}/mutscan_combined_lanes_R2.fastq.gz \
					-m ${TUBBY_PIPELINE}/bin/mutation_tubby.csv \
					-t 50 \
					-h ${TUBBY_PIPELINE}/results_mutscan/MUTSCAN_${renamed_fastq[$sample]}_htmlReport.html \
					--standalone \
					>  ${TUBBY_PIPELINE}/results_mutscan/MUTSCAN_${renamed_fastq[$sample]}_readsFound.txt 
	
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