

#!/bin/bash

##################################################################

echo ""
echo "============================================================="
echo "[Starting Bash Script] " ${__file} && date 
echo ""

STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./genomeDirectory --genomeFastaFiles ./mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile ./mouse_annotation/Mus_musculus.GRCm38.101.gtf --limitGenomeGenerateRAM=99687447936




echo ""
echo "[Ending Bash Script] " ${__file} && date 
echo "============================================================="
echo ""
