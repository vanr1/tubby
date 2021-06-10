#!/bin/bash

### change PatientID to Sample
sed -i 's/PatientID/Sample/g' ReadCounts.tsv

### to change up Daniel's ReadCounts file, move the 'type' column to second row and rename it to 'Type'

# extract last col
less ReadCounts.tsv | rev | cut -f1 | rev > lastCol.txt

# drop last col
# awk -F'\t' 'NF{NF-=1};1' OFS='\t' ReadCounts.tsv > ReadCountsWithoutLastCol.tsv
less ReadCounts.tsv | rev | cut -f2- | rev > ReadCountsWithoutLastCol.tsv

# change 'type' to 'Type'
sed 's/type/Type/g' lastCol.txt > lastColType.txt


paste <(cut -f1 ReadCountsWithoutLastCol.tsv) lastColType.txt <(cut -f2- ReadCountsWithoutLastCol.tsv) > ReadCountsRearranged.tsv


## [debug]
head -n 1 ReadCountsWithoutLastCol.tsv | awk '{print NF; exit}'
head -n 1 ReadCounts.tsv | awk '{print NF; exit}'
head -n 1 ReadCountsRearranged.tsv | awk '{print NF; exit}'
# awk '{print $NF}' ReadCountsRearranged.tsv | tail -5

### (now) ReadCounts.tsv in format wanted


### transpose (in bash)
python -c "import sys; print('\n'.join('\t'.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < ReadCountsRearranged.tsv > ReadCountsRearranged_transpose.tsv

### Subset only protein coding genes into file (~60000 genes -> ~19800 genes)

# get list of protein coding genes
less gencode.v22.annotation.gtf | grep "gene" | grep protein_coding | cut -f9 | cut -f2 -d'"' | sort -u > EnsembleIDsPCG.txt

# extract the protein coding genes from gene lengths file (DONT DO THIS, we still need all geneLengths)
# less geneLengthsForV22.txt | grep -Ff EnsembleIDsPCG.txt > geneLengthsForV22_PCG.txt
# wc -l geneLengthsForV22.txt
# wc -l geneLengthsForV22_PCG.txt

# add in "Sample" and "Type" so don't lose those data
cat <(echo Sample) <(echo Type) EnsembleIDsPCG.txt	> EnsembleIDsPCG_SampleType.txt

# extract the protein coding genes from our counts file
less ReadCountsRearranged_transpose.tsv | grep -Ff EnsembleIDsPCG_SampleType.txt > ReadCountsRearranged_transpose_PCG.tsv


### transpose back (in bash)
python -c "import sys; print('\n'.join('\t'.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < ReadCountsRearranged_transpose_PCG.tsv > ReadCountsPCG.tsv


## [debug]
wc -l ReadCounts.tsv && head -n 1 ReadCounts.tsv | awk '{print NF; exit}'
wc -l ReadCountsPCG.tsv && head -n 1 ReadCountsPCG.tsv | awk '{print NF; exit}'


# run script to convert counts to fpkm
Rscript convertCountsToRPKMonTable.r

