

#!/bin/bash

##################################################################

echo ""
echo "============================================================="
echo "[Starting Bash Script] " ${__file} && date 
echo ""

## <SPECIFIC TO> second set of raw data provided:
	# (atty location) /data2/han_lab/richardvan/tubby_pipeline/raw_fastq_2020november
ATTY_LOCATION="/data2/han_lab/richardvan/tubby_pipeline"




#NOTICE NEED TO CHANGE THE DATE IN 3 SPOTS
./nextflow run nf-core/rnaseq -r 3.0 --input samplesheet20210512.csv --outdir './results_20210512_24samples' --genome GRCm38 -profile singularity --pseudo_aligner salmon -resume 2>&1 | tee log_nextflow_results_20210512_24samples.log 



echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""