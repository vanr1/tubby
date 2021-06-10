
#### transpose
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}' merged_myAnalysis_20210202_18samples.txt


#### find value of salmon TPM file      /Users/richardvan/dev/sandboxForTubby/downloads/results_20210202_18samples_retinaVsHypothalamus/star_salmon
cat <(head -1 salmon.merged.gene_tpm.tsv) <(grep ENSMUSG00000022346 salmon.merged.gene_tpm.tsv) | awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}'

#### find value of normalized counts    /Users/richardvan/dev/sandboxForTubby/forDESEQ2_rcode
## need to add "geneid" as part of header at topleft first
cat <(head -1 data_normalized_counts.txt) <(grep ENSMUSG00000089789 data_normalized_counts.txt) | awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}'

## genes to look at
#   ---             ---         
# rhodopsin         ENSMUSG00000030324
# opsin             ENSMUSG00000058831
# agouti-related    ENSMUSG00000005705
# arginine          ENSMUSG00000037727
# tub               ENSMUSG00000031028
# gadph             ENSMUSG00000057666


#### find value of log transformed normalized counts -- log2(n+1)   /Users/richardvan/dev/sandboxForTubby/forDESEQ2_rcode
## need to add "geneid" as part of header at topleft first
cat <(head -1 data_rlog_counts.txt) <(grep ENSMUSG00000057666 data_rlog_counts.txt) | awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}'





### below is to print rows and columns

wc -l _FILE_ && head -1 _FILE_  | awk -F"\t" '{print NF; exit}'

wc -l salmon.merged.gene_tpm.tsv && head -1 salmon.merged.gene_tpm.tsv | awk -F"\t" '{print NF; exit}'

wc -l data_rlog_counts.txt && head -1 data_rlog_counts.txt | awk -F"\t" '{print NF; exit}'