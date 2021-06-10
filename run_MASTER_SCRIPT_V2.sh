

#!/bin/bash

##################################################################

echo ""
echo "============================================================="
echo "[Starting Bash Script] " ${__file} && date 
echo ""

## <SPECIFIC TO> second set of raw data provided:
	# (atty location) /data2/han_lab/richardvan/tubby_pipeline/raw_fastq_2020november
ATTY_LOCATION="/data2/han_lab/richardvan/tubby_pipeline/raw_fastq_2020november"





## (here have raw dataset in its own folder)
	# commands used inside "raw_fastq_2020november/RNASeq111920JohnNora-211372170/FASTQ_Generation_2020-11-20_14_02_18Z-344723379/"
		# -bash-4.2$ 	mv -v 1*Nora*/* ../../
			# ‘14Nora_L001-ds.6a6388b0f779497ca393818522c2c68c/14Nora_S15_L001_R1_001.fastq.gz’ -> ‘../../14Nora_S15_L001_R1_001.fastq.gz’
			# ‘14Nora_L001-ds.6a6388b0f779497ca393818522c2c68c/14Nora_S15_L001_R2_001.fastq.gz’ -> ‘../../14Nora_S15_L001_R2_001.fastq.gz’
			# ‘14Nora_L002-ds.7e69160a04e74d3cbfbe850f5d5a179b/14Nora_S15_L002_R1_001.fastq.gz’ -> ‘../../14Nora_S15_L002_R1_001.fastq.gz’
			# ‘14Nora_L002-ds.7e69160a04e74d3cbfbe850f5d5a179b/14Nora_S15_L002_R2_001.fastq.gz’ -> ‘../../14Nora_S15_L002_R2_001.fastq.gz’
			# ‘14Nora_L003-ds.d507e2656bbd434c9dc32ddb542d1e22/14Nora_S15_L003_R1_001.fastq.gz’ -> ‘../../14Nora_S15_L003_R1_001.fastq.gz’
			# ‘14Nora_L003-ds.d507e2656bbd434c9dc32ddb542d1e22/14Nora_S15_L003_R2_001.fastq.gz’ -> ‘../../14Nora_S15_L003_R2_001.fastq.gz’
			# ‘14Nora_L004-ds.14ff98c291244dc99bd2031b20621edf/14Nora_S15_L004_R1_001.fastq.gz’ -> ‘../../14Nora_S15_L004_R1_001.fastq.gz’
			# ‘14Nora_L004-ds.14ff98c291244dc99bd2031b20621edf/14Nora_S15_L004_R2_001.fastq.gz’ -> ‘../../14Nora_S15_L004_R2_001.fastq.gz’
			# ‘15Nora_L001-ds.85663ed59da049d38df6bb5703df2533/15Nora_S9_L001_R1_001.fastq.gz’ -> ‘../../15Nora_S9_L001_R1_001.fastq.gz’
			# ‘15Nora_L001-ds.85663ed59da049d38df6bb5703df2533/15Nora_S9_L001_R2_001.fastq.gz’ -> ‘../../15Nora_S9_L001_R2_001.fastq.gz’
			# ‘15Nora_L002-ds.c40bfd989cbc4f9b9eccada9883654be/15Nora_S9_L002_R1_001.fastq.gz’ -> ‘../../15Nora_S9_L002_R1_001.fastq.gz’
			# ‘15Nora_L002-ds.c40bfd989cbc4f9b9eccada9883654be/15Nora_S9_L002_R2_001.fastq.gz’ -> ‘../../15Nora_S9_L002_R2_001.fastq.gz’
			# ‘15Nora_L003-ds.9f386ab485a14ff29eb3be36cd53a798/15Nora_S9_L003_R1_001.fastq.gz’ -> ‘../../15Nora_S9_L003_R1_001.fastq.gz’
			# ‘15Nora_L003-ds.9f386ab485a14ff29eb3be36cd53a798/15Nora_S9_L003_R2_001.fastq.gz’ -> ‘../../15Nora_S9_L003_R2_001.fastq.gz’
			# ‘15Nora_L004-ds.5372c33a3d5d4c65bba396d846dab8ed/15Nora_S9_L004_R1_001.fastq.gz’ -> ‘../../15Nora_S9_L004_R1_001.fastq.gz’
			# ‘15Nora_L004-ds.5372c33a3d5d4c65bba396d846dab8ed/15Nora_S9_L004_R2_001.fastq.gz’ -> ‘../../15Nora_S9_L004_R2_001.fastq.gz’
			# ‘17Nora_L001-ds.2c2213edc2a04f3db272c96840aae971/17Nora_S17_L001_R1_001.fastq.gz’ -> ‘../../17Nora_S17_L001_R1_001.fastq.gz’
			# ‘17Nora_L001-ds.2c2213edc2a04f3db272c96840aae971/17Nora_S17_L001_R2_001.fastq.gz’ -> ‘../../17Nora_S17_L001_R2_001.fastq.gz’
			# ‘17Nora_L002-ds.3b4d1f7de77e4f3b9ca140cab109c815/17Nora_S17_L002_R1_001.fastq.gz’ -> ‘../../17Nora_S17_L002_R1_001.fastq.gz’
			# ‘17Nora_L002-ds.3b4d1f7de77e4f3b9ca140cab109c815/17Nora_S17_L002_R2_001.fastq.gz’ -> ‘../../17Nora_S17_L002_R2_001.fastq.gz’
			# ‘17Nora_L003-ds.603ef8cca5c0411dbd087fe54c3c1c0e/17Nora_S17_L003_R1_001.fastq.gz’ -> ‘../../17Nora_S17_L003_R1_001.fastq.gz’
			# ‘17Nora_L003-ds.603ef8cca5c0411dbd087fe54c3c1c0e/17Nora_S17_L003_R2_001.fastq.gz’ -> ‘../../17Nora_S17_L003_R2_001.fastq.gz’
			# ‘17Nora_L004-ds.96536362aab54c90abf3e62d76a2ccaa/17Nora_S17_L004_R1_001.fastq.gz’ -> ‘../../17Nora_S17_L004_R1_001.fastq.gz’
			# ‘17Nora_L004-ds.96536362aab54c90abf3e62d76a2ccaa/17Nora_S17_L004_R2_001.fastq.gz’ -> ‘../../17Nora_S17_L004_R2_001.fastq.gz’
		# -bash-4.2$ mv -v 2*Nora*/* ../../
			# ‘21Nora_L001-ds.27455709a8834549845a0e46bf733360/21Nora_S4_L001_R1_001.fastq.gz’ -> ‘../../21Nora_S4_L001_R1_001.fastq.gz’
			# ‘21Nora_L001-ds.27455709a8834549845a0e46bf733360/21Nora_S4_L001_R2_001.fastq.gz’ -> ‘../../21Nora_S4_L001_R2_001.fastq.gz’
			# ‘21Nora_L002-ds.107a19e5028c4204b23b54ab216683e8/21Nora_S4_L002_R1_001.fastq.gz’ -> ‘../../21Nora_S4_L002_R1_001.fastq.gz’
			# ‘21Nora_L002-ds.107a19e5028c4204b23b54ab216683e8/21Nora_S4_L002_R2_001.fastq.gz’ -> ‘../../21Nora_S4_L002_R2_001.fastq.gz’
			# ‘21Nora_L003-ds.0f5139347028487e8c776ba2a051d6cb/21Nora_S4_L003_R1_001.fastq.gz’ -> ‘../../21Nora_S4_L003_R1_001.fastq.gz’
			# ‘21Nora_L003-ds.0f5139347028487e8c776ba2a051d6cb/21Nora_S4_L003_R2_001.fastq.gz’ -> ‘../../21Nora_S4_L003_R2_001.fastq.gz’
			# ‘21Nora_L004-ds.8aaa160a75b6499d8bf9e3f5016647ab/21Nora_S4_L004_R1_001.fastq.gz’ -> ‘../../21Nora_S4_L004_R1_001.fastq.gz’
			# ‘21Nora_L004-ds.8aaa160a75b6499d8bf9e3f5016647ab/21Nora_S4_L004_R2_001.fastq.gz’ -> ‘../../21Nora_S4_L004_R2_001.fastq.gz’
			# ‘23Nora_L001-ds.cd722aaed6274a5f87777f54db1a690f/23Nora_S2_L001_R1_001.fastq.gz’ -> ‘../../23Nora_S2_L001_R1_001.fastq.gz’
			# ‘23Nora_L001-ds.cd722aaed6274a5f87777f54db1a690f/23Nora_S2_L001_R2_001.fastq.gz’ -> ‘../../23Nora_S2_L001_R2_001.fastq.gz’
			# ‘23Nora_L002-ds.1757c22c7ae1479785d3d93fa875289f/23Nora_S2_L002_R1_001.fastq.gz’ -> ‘../../23Nora_S2_L002_R1_001.fastq.gz’
			# ‘23Nora_L002-ds.1757c22c7ae1479785d3d93fa875289f/23Nora_S2_L002_R2_001.fastq.gz’ -> ‘../../23Nora_S2_L002_R2_001.fastq.gz’
			# ‘23Nora_L003-ds.b999f3515e2d4d2dbb20ce94aebb3a72/23Nora_S2_L003_R1_001.fastq.gz’ -> ‘../../23Nora_S2_L003_R1_001.fastq.gz’
			# ‘23Nora_L003-ds.b999f3515e2d4d2dbb20ce94aebb3a72/23Nora_S2_L003_R2_001.fastq.gz’ -> ‘../../23Nora_S2_L003_R2_001.fastq.gz’
			# ‘23Nora_L004-ds.9a7a327249d541ca810e5d831259ec68/23Nora_S2_L004_R1_001.fastq.gz’ -> ‘../../23Nora_S2_L004_R1_001.fastq.gz’
			# ‘23Nora_L004-ds.9a7a327249d541ca810e5d831259ec68/23Nora_S2_L004_R2_001.fastq.gz’ -> ‘../../23Nora_S2_L004_R2_001.fastq.gz’
			# ‘24Nora_L001-ds.8b63ef8168864aefaa00d6f7a3f775c1/24Nora_S13_L001_R1_001.fastq.gz’ -> ‘../../24Nora_S13_L001_R1_001.fastq.gz’
			# ‘24Nora_L001-ds.8b63ef8168864aefaa00d6f7a3f775c1/24Nora_S13_L001_R2_001.fastq.gz’ -> ‘../../24Nora_S13_L001_R2_001.fastq.gz’
			# ‘24Nora_L002-ds.7718a3c06b1e4950bc7a7e562a0aac39/24Nora_S13_L002_R1_001.fastq.gz’ -> ‘../../24Nora_S13_L002_R1_001.fastq.gz’
			# ‘24Nora_L002-ds.7718a3c06b1e4950bc7a7e562a0aac39/24Nora_S13_L002_R2_001.fastq.gz’ -> ‘../../24Nora_S13_L002_R2_001.fastq.gz’
			# ‘24Nora_L003-ds.4c625b2d8b5d4f2182bd3990c5eb13ff/24Nora_S13_L003_R1_001.fastq.gz’ -> ‘../../24Nora_S13_L003_R1_001.fastq.gz’
			# ‘24Nora_L003-ds.4c625b2d8b5d4f2182bd3990c5eb13ff/24Nora_S13_L003_R2_001.fastq.gz’ -> ‘../../24Nora_S13_L003_R2_001.fastq.gz’
			# ‘24Nora_L004-ds.d4757799c6d74501b79bf8dc35ef836f/24Nora_S13_L004_R1_001.fastq.gz’ -> ‘../../24Nora_S13_L004_R1_001.fastq.gz’
			# ‘24Nora_L004-ds.d4757799c6d74501b79bf8dc35ef836f/24Nora_S13_L004_R2_001.fastq.gz’ -> ‘../../24Nora_S13_L004_R2_001.fastq.gz’

## (here checking if is pair end reads)
	# commands used inside "raw_fastq_2020november/"
		# -bash-4.2$ 	wc -l *
		   # 31094236 14Nora_S15_L001_R1_001.fastq
		   # 31094236 14Nora_S15_L001_R2_001.fastq
		   # 30700656 14Nora_S15_L002_R1_001.fastq
		   # 30700656 14Nora_S15_L002_R2_001.fastq
		   # 30890052 14Nora_S15_L003_R1_001.fastq
		   # 30890052 14Nora_S15_L003_R2_001.fastq
		   # 30761964 14Nora_S15_L004_R1_001.fastq
		   # 30761964 14Nora_S15_L004_R2_001.fastq
		   # 25950904 15Nora_S9_L001_R1_001.fastq
		   # 25950904 15Nora_S9_L001_R2_001.fastq
		   # 25665328 15Nora_S9_L002_R1_001.fastq
		   # 25665328 15Nora_S9_L002_R2_001.fastq
		   # 25866996 15Nora_S9_L003_R1_001.fastq
		   # 25866996 15Nora_S9_L003_R2_001.fastq
		   # 25737708 15Nora_S9_L004_R1_001.fastq
		   # 25737708 15Nora_S9_L004_R2_001.fastq
		   # 37286720 17Nora_S17_L001_R1_001.fastq
		   # 37286720 17Nora_S17_L001_R2_001.fastq
		   # 36972180 17Nora_S17_L002_R1_001.fastq
		   # 36972180 17Nora_S17_L002_R2_001.fastq
		   # 37125540 17Nora_S17_L003_R1_001.fastq
		   # 37125540 17Nora_S17_L003_R2_001.fastq
		   # 36999592 17Nora_S17_L004_R1_001.fastq
		   # 36999592 17Nora_S17_L004_R2_001.fastq
		   # 36937172 21Nora_S4_L001_R1_001.fastq
		   # 36937172 21Nora_S4_L001_R2_001.fastq
		   # 36588068 21Nora_S4_L002_R1_001.fastq
		   # 36588068 21Nora_S4_L002_R2_001.fastq
		   # 36809712 21Nora_S4_L003_R1_001.fastq
		   # 36809712 21Nora_S4_L003_R2_001.fastq
		   # 36642152 21Nora_S4_L004_R1_001.fastq
		   # 36642152 21Nora_S4_L004_R2_001.fastq
		   # 35600036 23Nora_S2_L001_R1_001.fastq
		   # 35600036 23Nora_S2_L001_R2_001.fastq
		   # 35146932 23Nora_S2_L002_R1_001.fastq
		   # 35146932 23Nora_S2_L002_R2_001.fastq
		   # 35353496 23Nora_S2_L003_R1_001.fastq
		   # 35353496 23Nora_S2_L003_R2_001.fastq
		   # 35171000 23Nora_S2_L004_R1_001.fastq
		   # 35171000 23Nora_S2_L004_R2_001.fastq
		   # 34681240 24Nora_S13_L001_R1_001.fastq
		   # 34681240 24Nora_S13_L001_R2_001.fastq
		   # 34328476 24Nora_S13_L002_R1_001.fastq
		   # 34328476 24Nora_S13_L002_R2_001.fastq
		   # 34494636 24Nora_S13_L003_R1_001.fastq
		   # 34494636 24Nora_S13_L003_R2_001.fastq
		   # 34365756 24Nora_S13_L004_R1_001.fastq
		   # 34365756 24Nora_S13_L004_R2_001.fastq



# (combine lanes) so get total for 6 samples

		## (EX) next two lines is from the first dataset provided
				# cat 7WT24-R1_S5_L00*_R1_001.fastq > 7WT24-R1_S5_LALL_R1_001.fastq		
				# cat 7WT24-R1_S5_L00*_R2_001.fastq > 7WT24-R1_S5_LALL_R2_001.fastq		

# echo ""
# echo "   [myMessage] about to merge the lanes of the second provided dataset of rna-seq data "
# echo "   [myMessage] 			notice removing 001 that ends every file "

# echo ""
# cat -v ${ATTY_LOCATION}/14Nora_S15_L00*_R1_001.fastq > ${ATTY_LOCATION}/14Nora_S15_LALL_R1_001.fastq	
# cat -v ${ATTY_LOCATION}/14Nora_S15_L00*_R2_001.fastq > ${ATTY_LOCATION}/14Nora_S15_LALL_R2_001.fastq	

# cat -v ${ATTY_LOCATION}/15Nora_S9_L00*_R1_001.fastq > ${ATTY_LOCATION}/15Nora_S9_LALL_R1_001.fastq	
# cat -v ${ATTY_LOCATION}/15Nora_S9_L00*_R2_001.fastq	 > ${ATTY_LOCATION}/15Nora_S9_LALL_R2_001.fastq	

# cat -v ${ATTY_LOCATION}/17Nora_S17_L00*_R1_001.fastq > ${ATTY_LOCATION}/17Nora_S17_LALL_R1_001.fastq
# cat -v ${ATTY_LOCATION}/17Nora_S17_L00*_R2_001.fastq > ${ATTY_LOCATION}/17Nora_S17_LALL_R2_001.fastq

# cat -v ${ATTY_LOCATION}/21Nora_S4_L00*_R1_001.fastq > ${ATTY_LOCATION}/21Nora_S4_LALL_R1_001.fastq
# cat -v ${ATTY_LOCATION}/21Nora_S4_L00*_R2_001.fastq > ${ATTY_LOCATION}/21Nora_S4_LALL_R2_001.fastq

# cat -v ${ATTY_LOCATION}/23Nora_S2_L00*_R1_001.fastq > ${ATTY_LOCATION}/23Nora_S2_LALL_R1_001.fastq	
# cat -v ${ATTY_LOCATION}/23Nora_S2_L00*_R2_001.fastq > ${ATTY_LOCATION}/23Nora_S2_LALL_R2_001.fastq	

# cat -v ${ATTY_LOCATION}/24Nora_S13_L00*_R1_001.fastq > ${ATTY_LOCATION}/24Nora_S13_LALL_R1_001.fastq
# cat -v ${ATTY_LOCATION}/24Nora_S13_L00*_R2_001.fastq > ${ATTY_LOCATION}/24Nora_S13_LALL_R2_001.fastq
# echo ""
# echo "   [myMessage] There were 6 samples total; each had pair-end reads due to same sizes; merged all 4 lanes into *_LALL_*"
# echo ""






# (use samples where lanes have been merged)
declare -a arrayOfSampleBaseNames=(	
		"raw_fastq_2020november/14Nora_S15"
		"raw_fastq_2020november/15Nora_S9" 		
		"raw_fastq_2020november/17Nora_S17" 	
		"raw_fastq_2020november/21Nora_S4" 		
		"raw_fastq_2020november/23Nora_S2" 	
		"raw_fastq_2020november/24Nora_S13" 		
		)



# # (run FAST_QC) on each of these 6 samples (twice because paired-end reads); this command takes ~10 minutes per sample; so 60 minutes total (run in val node)
#   # (why) we use fastqc to analyse pre-mapped fastq reads for R1 and R2 separetly

# mkdir -p results_2020november/FASTQCresults
# output="/data2/han_lab/richardvan/tubby_pipeline/results_2020november/FASTQCresults"

# for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
# do
# 	echo "  [START processing sample using FASTQC] " ${sampleBaseName} && date 
# 	fastqc -f fastq -o ${output} ${sampleBaseName}_LALL_R1.fastq
# 	fastqc -f fastq -o ${output} ${sampleBaseName}_LALL_R2.fastq


# done
# #   (copy to own cpu to view on browser with) scp -r richardvan@atty.nipm.unlv.edu:/data2/han_lab/richardvan/tubby_pipeline/results_2020november .



# (run SALMON) on each
for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
do
	echo "  [START working on sample using SALMON] " ${sampleBaseName} && date 

	salmon quant --threads 40 \
				 -i transcripts_index \
				 --libType A \
				 	-1 ${sampleBaseName}_LALL_R1_001.fastq \
				 	-2 ${sampleBaseName}_LALL_R2_001.fastq \
				 --gcBias \
				 --validateMappings -o transcripts_quant_${sampleBaseName}
done



echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
