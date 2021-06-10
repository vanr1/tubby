#!/bin/bash

# calcualte gene length command (needs to be done once)
# Rscript ./calculateGeneLengthsFromAnnotation.R

# ./trimPairEndReadsAndRunFastQC.sh 
# ./runSTARcommandToMap.sh 
# ./runHTSEQtoCountReadsInFeatures.sh 2>&1 | tee log_runHTSEQtoCountReadsInFeatures_20201004.log



###########################
# convert all counts to rpkm
###########################
mkdir -p tmp
declare -a arrayOfSampleBaseNames=(	"9WT24-R3_S1" 
									"8WT24-R2_S2" 
									"3Tub24-R3_S3" 
									"1Tub24-R1_S4" 
									"7WT24-R1_S5" 
									"2Tub24-R2_S6" 
								)

for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
do
	date
	echo "  [START working on sample] " ${sampleBaseName} 

	# This says, "find all lines that start with # and delete them, leaving everything else."
	sed '/^__/d' ${sampleBaseName}_counts.txt > tmp/${sampleBaseName}_counts_minusLastLines.txt

	echo "    converting from counts to rpkm" 
	Rscript ./convertCountsToRPKM.r tmp/${sampleBaseName}_counts_minusLastLines.txt ./Mus_musculus.GRCm38.geneLength.txt tmp/${sampleBaseName}_convertedRpkm.txt

	echo "    cutting out only the RPKM column and changing the first line to sample name" 
	cut -f4 tmp/${sampleBaseName}_convertedRpkm.txt > tmp/${sampleBaseName}_rpkmOnly.txt 	
	sed "s/RPKM/${sampleBaseName}/g" tmp/${sampleBaseName}_rpkmOnly.txt > tmp/${sampleBaseName}_rpkm.txt 

	# echo "    cutting out only the counts column and changing the first line to sample name" 
	# cut -f3 tmp/${sampleBaseName}_convertedRpkm.txt > tmp/${sampleBaseName}_countOnly.txt 	
	# sed "s/Count/${sampleBaseName}/g" tmp/${sampleBaseName}_countOnly.txt > tmp/${sampleBaseName}_counts.txt 	# just a longer way of adding the sample name to the top of count file

	echo " " 
done


###########################
# merge only the rpkm values into a table; add headers and stuff
###########################

# create the initial file that is just gene names; cut it from a file above
cut -f1 tmp/9WT24-R3_S1_convertedRpkm.txt > tmp/mergedFile.txt

for sampleBaseName in "${arrayOfSampleBaseNames[@]}"
do
	date
	echo "  [pasting sample of rpkm values] " ${sampleBaseName} 
	paste tmp/mergedFile.txt tmp/${sampleBaseName}_rpkm.txt > tmp/newMergedFile.txt
	mv tmp/newMergedFile.txt tmp/mergedFile.txt 

	echo " " 
done

# move the final merged file to the top folder and name it what you want
cp -av tmp/mergedFile.txt tubby_samples_rpkm.txt


