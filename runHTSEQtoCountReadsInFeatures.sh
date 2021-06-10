

#!/bin/bash

##################################################################

echo ""
echo "============================================================="
echo "[Starting Bash Script] " ${__file} && date 
echo ""


# making a variable for our main pipeline path
TUBBY_PIPELINE=/data2/han_lab/richardvan/tubby_pipeline

set -e  #exits the script when an error occurs for safety

declare -a arrayOfSampleBaseNames=( 

                        "1Tub24-R1_S4"       # first batch (oct 2020)
                        "2Tub24-R2_S6"       
                        "3Tub24-R3_S3"    
                        "7WT24-R1_S5"
                        "8WT24-R2_S2"     
                        "9WT24-R3_S1"        

                        "14Nora_S15"      # second batch (dec 2020)
                        "15Nora_S9"       
                        "17Nora_S17"   
                        "21Nora_S4"       
                        "23Nora_S2"    
                        "24Nora_S13" 

                                       # third batch (october 2020)
                        )


htseq-count ./Aligned.sortedByCoord.out.bam ./mouse_annotation/Mus_musculus.GRCm38.101.gtf --format=bam --order=pos





echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
