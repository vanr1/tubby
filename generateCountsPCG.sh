#!/bin/bash


## [debug]
wc -l 7WT24-R1_S5_counts.txt && head -n 1 7WT24-R1_S5_counts.txt | awk '{print NF; exit}'

### (now) ReadCounts.txt in format wanted


### transpose (in bash)
# python -c "import sys; print('\n'.join('\t'.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < 7WT24-R1_S5_counts.txt > tmp_counts_transposed

### Subset only protein coding genes into file (~60000 genes -> ~19800 genes)

# get list of protein coding genes
less ./mouse_annotation/Mus_musculus.GRCm38.101.gtf | grep "gene" | grep protein_coding | cut -f9 | cut -f2 -d'"' | sort -u > tmp_EnsembleIDsPCG.txt

# extract the protein coding genes from gene lengths file (DONT DO THIS, we still need all geneLengths)
# less geneLengthsForV22.txt | grep -Ff EnsembleIDsPCG.txt > geneLengthsForV22_PCG.txt
# wc -l geneLengthsForV22.txt
# wc -l geneLengthsForV22_PCG.txt

# extract the protein coding genes from our counts file
less 7WT24-R1_S5_counts.txt  | grep --file tmp_EnsembleIDsPCG.txt > 7WT24-R1_S5_countsPCG.txt




## [debug]
wc -l 7WT24-R1_S5_counts.txt && head -n 1 7WT24-R1_S5_counts.txt | awk '{print NF; exit}'
wc -l 7WT24-R1_S5_countsPCG.txt && head -n 1 7WT24-R1_S5_countsPCG.txt | awk '{print NF; exit}'


# run script to convert counts to fpkm
# Rscript convertCountsToRPKMonTable.r

