

#!/bin/bash

##################################################################

echo ""
echo "============================================================="
echo "[Starting Bash Script] " ${__file} && date 
echo ""

featureCounts -p -t exon -g gene_id -a mouse_annotation/Mus_musculus.GRCm38.101.gtf -o Aligned.counts.txt out.sam 2>&1



echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
