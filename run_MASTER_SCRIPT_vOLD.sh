

#!/bin/bash

##################################################################

echo ""
echo "============================================================="
echo "[Starting Bash Script] " ${__file} && date 
echo ""

## (here have raw dataset in its own folder




## <SPECIFIC TO> first set of raw data provided:
	# (atty location) /data2/han_lab/richardvan/tubby_pipeline/raw_fastq

# (combine lanes) so get total for 6 samples; 3 replicates per experimental/control group

cat 7WT24-R1_S5_L00*_R1_001.fastq > 7WT24-R1_S5_LALL_R1_001.fastq		#done
cat 7WT24-R1_S5_L00*_R2_001.fastq > 7WT24-R1_S5_LALL_R2_001.fastq		#done

cat 8WT24-R2_S2_L00*_R1_001.fastq > 8WT24-R2_S2_LALL_R1_001.fastq		#done
cat 8WT24-R2_S2_L00*_R2_001.fastq > 8WT24-R2_S2_LALL_R2_001.fastq		#done

cat 9WT24-R3_S1_L00*_R1_001.fastq > 9WT24-R3_S1_LALL_R1_001.fastq		#done
cat 9WT24-R3_S1_L00*_R2_001.fastq > 9WT24-R3_S1_LALL_R2_001.fastq		#done


cat 1Tub24-R1_S4_L00*_R1_001.fastq > 1Tub24-R1_S4_LALL_R1_001.fastq		#done
cat 1Tub24-R1_S4_L00*_R2_001.fastq > 1Tub24-R1_S4_LALL_R2_001.fastq		#done

cat 2Tub24-R2_S6_L00*_R1_001.fastq > 2Tub24-R2_S6_LALL_R1_001.fastq		#done
cat 2Tub24-R2_S6_L00*_R2_001.fastq > 2Tub24-R2_S6_LALL_R2_001.fastq		#done

cat 3Tub24-R3_S3_L00*_R1_001.fastq > 3Tub24-R3_S3_LALL_R1_001.fastq		#done
cat 3Tub24-R3_S3_L00*_R2_001.fastq > 3Tub24-R3_S3_LALL_R2_001.fastq		#done


featureCounts -p -t exon -g gene_id -a mouse_annotation/Mus_musculus.GRCm38.101.gtf -o Aligned.counts.txt out.sam 2>&1



echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
