

#!/bin/bash

##################################################################

echo ""
echo "============================================================="
echo "[Starting Bash Script] " ${__file} && date 
echo ""


# making a variable for our main pipeline path
TUBBY_PIPELINE=/data2/han_lab/richardvan/tubby_pipeline

set -e  #exits the script when an error occurs for safety

mkdir -p ${TUBBY_PIPELINE}/data/star_index/genome_directory
STAR 	--runThreadN 80 --runMode genomeGenerate \
		--genomeDir ${TUBBY_PIPELINE}/data/star_index/genome_directory \
		--genomeFastaFiles ${TUBBY_PIPELINE}/data/ref_genome/mouse_ensemble_GRCm38p6.fasta \
		--sjdbGTFfile ${TUBBY_PIPELINE}/data/ref_genome/mouse_ensemble_GRCm38p6.gff3 \
		--limitGenomeGenerateRAM=83476462976 \
		--sjdbOverhang 74




echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
