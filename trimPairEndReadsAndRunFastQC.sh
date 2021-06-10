

#!/bin/bash

##################################################################

echo ""
echo "============================================================="
echo "[Starting Bash Script] " ${__file} && date 
echo ""


# (trim basic parameneters)
java -jar ./trimmomatic-0.39.jar PE ./raw_fastq/7WT24-R1_S5_LALL_R1_001.fastq  ./raw_fastq/7WT24-R1_S5_LALL_R2_001.fastq ./trimmed_fastq/7WT24-R1_S5_LALL_R1_001.forward_paired.fastq ./trimmed_fastq/7WT24-R1_S5_LALL_R1_001.forward_unpaired.fastq ./trimmed_fastq/7WT24-R1_S5_LALL_R1_001.reverse_paired.fastq ./trimmed_fastq/7WT24-R1_S5_LALL_R1_001.reverse_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


# (trim down to 50bp)
# java -jar ./trimmomatic-0.39.jar PE ./raw_fastq/7WT24-R1_S5_LALL_R1_001.fastq  ./raw_fastq/7WT24-R1_S5_LALL_R2_001.fastq ./trimmed_fastq/7WT24-R1_S5_LALL_R1_001.forward_paired.fastq ./trimmed_fastq/7WT24-R1_S5_LALL_R1_001.forward_unpaired.fastq ./trimmed_fastq/7WT24-R1_S5_LALL_R1_001.reverse_paired.fastq ./trimmed_fastq/7WT24-R1_S5_LALL_R1_001.reverse_unpaired.fastq LEADING:10 TRAILING:10 CROP:50 SLIDINGWINDOW:4:15 MINLEN:36

# fastqc ./trimmed_fastq/*.fastq




echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
