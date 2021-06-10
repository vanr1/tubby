# _TODO_ check if this below output is the same as the file i'm looking at??

library("DESeq2")
library("edgeR")
library("ggplot2")
library("eulerr")
library("UpSetR")

## Constants
P_ALPHA_CUTOFF = 0.05
P_LOG2FC_THRESHOLD = 1.0


## BioMart Setup
print ("[DEBUG] using BIOMART library to get gene symbol and descriptions")
library( "biomaRt" )
ENSEMBL_BIOMART = useMart( "ensembl", dataset = "mmusculus_gene_ensembl", host="uswest.ensembl.org")
possibleAttribute = listAttributes(ENSEMBL_BIOMART)     # choose attributes of your interest
write.table(
  data.frame(possibleAttribute), 
  file="possible_ENSEMBL_BIOMART_attributes.txt", 
  sep="\t", 
  quote=F, 
  row.names=FALSE)





## values coming out of deseq2.dds
# load("18samples_retinaVsHypothalamus.deseq2.dds.RData")
load("24samples_retinaVsHypothalamus.deseq2.dds.RData")
dim(dds)

## first subset the raw dds by getting rid of genes with NA
numericData = counts(dds)
nonzeroData = numericData[rowSums(numericData==0, na.rm=TRUE)<ncol(numericData), ]
nonzeroRowNames = rownames(nonzeroData)
nonzero_dds = dds[nonzeroRowNames, ]


## *** renaming back to 'dds' (so now this variable has nonzero rows)
dds = nonzero_dds
dim(dds)

### understanding "size factors" a type of normalization factor
                  #<REMINDER> this type of normalization for "gene count comparisons 
                  #           between samples and DE analysis; NOT for within sample comparisons

estimateSizeFactors(dds)
sizeFactors(dds)[1:4]         # function from DEseq2 
                                      # hypothalamus_R1 hypothalamus_R2 hypothalamus_R3 hypothalamus_R4 
                                      # 1.0438106       1.9372076       0.6550678       1.9910981 
      
head(counts(dds))[,1:4]         # function from DEseq2
                  #                     hypothalamus_R1 hypothalamus_R2 hypothalamus_R3 hypothalamus_R4
                  # ENSMUSG00000000001             346            1236             187            1191
                  # ENSMUSG00000000003               0               0               0               0
                  # ENSMUSG00000000028              25              55              17              49
                  # ENSMUSG00000000031              45              87              13              59
                  # ENSMUSG00000000037              18              63              11              81
                  # ENSMUSG00000000049               6               9               2               9

normalized_counts <- counts(dds, normalized= TRUE)   # need to reorder it to make it easier for layer (next line)
normalized_counts <- subset(normalized_counts, select=
         c("retina_R1","retina_R2","retina_R3","retina_R4","retina_R5","retina_R6",
           "retina_R7","retina_R8","retina_R9","retina_R10","retina_R11","retina_R12",
           "hypothalamus_R1","hypothalamus_R2","hypothalamus_R3","hypothalamus_R4","hypothalamus_R5","hypothalamus_R6",
           "hypothalamus_R7","hypothalamus_R8","hypothalamus_R9","hypothalamus_R10","hypothalamus_R11","hypothalamus_R12"
         ))
head(normalized_counts)[,1:4] 

normalized_counts_retina_3months <- subset(normalized_counts, select=
                              c("retina_R4","retina_R5","retina_R6",
                                "retina_R10","retina_R11","retina_R12"
                              ))

head(normalized_counts_retina_3months)[,]  



all_genes <- sapply( strsplit( rownames(normalized_counts), split="\\+" ), "[", 1)
all_genes_retina_3months <- sapply( strsplit( rownames(normalized_counts_retina_3months), split="\\+" ), "[", 1)
genemap <- getBM(  attributes = c("ensembl_gene_id","external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = all_genes,
                   mart = ENSEMBL_BIOMART)
genemap_retina_3months <- getBM(  attributes = c("ensembl_gene_id","external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = all_genes_retina_3months,
                   mart = ENSEMBL_BIOMART)

idx <- match( all_genes, genemap$ensembl_gene_id )
idx_retina_3months <- match( all_genes_retina_3months, genemap$ensembl_gene_id )


gene_name <- genemap$external_gene_name[ idx ]
gene_name_3months <- genemap_retina_3months$external_gene_name[ idx ]
# head(all_genes, 4)
# head(gene_name, 4)

# add genename to first column, genedescription to last
normalized_counts = cbind(gene_name, normalized_counts) 
normalized_counts_retina_3months = cbind(gene_name_3months, normalized_counts_retina_3months) 

head(normalized_counts)[,1:4]       

write.table(
  data.frame("geneid"=rownames(normalized_counts),normalized_counts), 
  file="data_normalized_counts.txt", 
  sep="\t", 
  quote=F, 
  row.names=FALSE)
write.table(
    data.frame("geneid"=rownames(normalized_counts_retina_3months),normalized_counts_retina_3months), 
    file=paste("output/retina_3mo/COUNTS_normalized_retina_3mo_DATA.txt",
               sep=""),
    sep="\t", 
    quote=F, 
    row.names=FALSE)




#### now seeing what rlog does to get the NEXTFLOW PCA plot
head(rlog(dds))
                  # class: DESeqTransform 
                  # dim: 45706 18 
                  # metadata(1): version
                  # assays(1): ''
                  # rownames(45706): ENSMUSG00000000001 ENSMUSG00000000003 ... ENSMUSG00000107390 ENSMUSG00000107392
                  # rowData names(24): baseMean baseVar ... dispFit rlogIntercept
                  # colnames(18): hypothalamus_R1 hypothalamus_R2 ... retina_R8 retina_R9
                  # colData names(3): condition sizeFactor replaceable
dim(dds)
                  # [1] 45706    18

dimnames(dds[1:6,1:4])
                  # [[1]]
                  # [1] "ENSMUSG00000000001" "ENSMUSG00000000003" "ENSMUSG00000000028"
                  # 
                  # [[2]]
                  # [1] "hypothalamus_R1" "hypothalamus_R2" "hypothalamus_R3" "hypothalamus_R4"
                  # ] "hypothalamus_R1" "hypothalamus_R2" "hypothalamus_R3" "hypothalamus_R4"


colData(dds)[2]
                  # DataFrame with 18 rows and 3 columns
                  # condition sizeFactor replaceable
                  # <factor>  <numeric>   <logical>
                  # hypothalamus_R1 hypothalamus   1.043811        TRUE
                  # hypothalamus_R2 hypothalamus   1.937208        TRUE
                  # hypothalamus_R3 hypothalamus   0.655068        TRUE
                  # hypothalamus_R4 hypothalamus   1.991098        TRUE
                  # hypothalamus_R5 hypothalamus   1.803478        TRUE
                  # ...                      ...        ...         ...
                  # retina_R5             retina   0.557446        TRUE
                  # retina_R6             retina   1.242944        TRUE
                  # retina_R7             retina   0.816152        TRUE
                  # retina_R8             retina   1.071291        TRUE
                  # retina_R9             retina   1.009718        TRUE

rowData(dds[1:4])
                  # DataFrame with 4 rows and 23 columns
                  # baseMean   baseVar   allZero dispGeneEst dispGeneIter   dispFit dispersion  dispIter dispOutlier   dispMAP Intercept
                  # <numeric> <numeric> <logical>   <numeric>    <numeric> <numeric>  <numeric> <numeric>   <logical> <numeric> <numeric>
                  # ENSMUSG00000000001  480.1399 11838.850     FALSE   0.0535878           10 0.0889105  0.0577923         7       FALSE 0.0577923   8.97073
                  # ENSMUSG00000000003    0.0000     0.000      TRUE          NA           NA        NA         NA        NA          NA        NA        NA
                  # ENSMUSG00000000028   37.8296   344.993     FALSE   0.0453664            8 0.1514059  0.0628013         9       FALSE 0.0628013   4.66096
                  # ENSMUSG00000000031   27.1040   336.804     FALSE   0.5335006            7 0.1782517  0.4580958        10       FALSE 0.4580958   4.89787


dds[,1:4]$sizeFactor
                  # hypothalamus_R1 hypothalamus_R2 hypothalamus_R3 hypothalamus_R4 
                  # 1.0438106       1.9372076       0.6550678       1.9910981 
                  # <same as rlog(dds) so transforming this data doesn change sizeFactor


counts(dds[1:3,1:4])
                # hypothalamus_R1 hypothalamus_R2 hypothalamus_R3 hypothalamus_R4
                # ENSMUSG00000000001             346            1236             187            1191
                # ENSMUSG00000000003               0               0               0               0
                # ENSMUSG00000000028              25              55              17              49



useThisLibraryToGetDEGsFromInputList <- function(argInputList, argLibrary, argExpCond) {
  print (      "  <# arguments for function call>")
  print (paste("    argInputListNumberOfGenes --- ", length(row.names(argInputList) )))
  print (paste("    argLibrary --- ", argLibrary))
  print (paste("    argExpCond --- ", argExpCond))
  
  
  print ("[DEBUG] using BIOMART library to get gene symbol and descriptions")
  argInputList$ensembl <- sapply( strsplit( rownames(argInputList), split="\\+" ), "[", 1)
  genemap <- getBM(  attributes = c("ensembl_gene_id","external_gene_name","interpro_description", "description","refseq_mrna"),
                     filters = "ensembl_gene_id",
                     values = argInputList$ensembl,
                     mart = ENSEMBL_BIOMART)
  idx <- match( argInputList$ensembl, genemap$ensembl_gene_id )
  argInputList$gene_name <- genemap$external_gene_name[ idx ]
  argInputList$gene_description_interpro <- genemap$interpro_description[ idx ]
  argInputList$gene_description_hgnc <- genemap$description[ idx ]
  argInputList$refseq <- genemap$refseq_mrna[ idx ]
  head(argInputList, 4)
  
  print ("[DEBUG] subset to get top down/up-regulated genes")
  resSig = argInputList [ which(argInputList$padj < P_ALPHA_CUTOFF ), ]             # fdr correction
  dim(resSig)
  # head (resSig[ order( resSig$log2FoldChange ), ] )   # for downregulated genes  (double-check not flipped; how?)
  # head (rev(resSig[ order( resSig$log2FoldChange ), ] ))   # for upregulated genes (double-check not flipped; how?)
  
  print ("[DEBUG] subset again to get down/up-regulated genes that are above the threshold")
  resSig_pos = resSig [ which(resSig$log2FoldChange > P_LOG2FC_THRESHOLD ), ]             # upregulated DEGs have log2FC > 1.5
  resSig_neg = resSig [ which(resSig$log2FoldChange < (-1*P_LOG2FC_THRESHOLD) ), ]        # downregulated DEGs have log2FC < -1.5
  
  
  print ("[DEBUG] save info about all genes to some tables")
  upreg_Log2FC_genes = resSig_pos[ order( resSig_pos$log2FoldChange ,decreasing=TRUE), ] 
  downreg_Log2FC_genes = resSig_neg[ order( resSig_neg$log2FoldChange ,decreasing=FALSE), ]
  
  allLog2FC_genes = argInputList[ order( argInputList$padj ,decreasing=FALSE), ] 
  row.has.na <- apply(allLog2FC_genes, 1, function(x){any(is.na(x))})
  allLog2FC_genes_filtered = allLog2FC_genes[!row.has.na,] 
  
  
  print (      "  <# of genes in certain tables current analysis >")
  print (paste("    original # genes begin --- ", length(row.names(argInputList) )))
  print (paste("    droppingNAs --- ", sum(row.has.na) ))
  print (paste("    total genes left (non-NA) --- ", length(row.names(allLog2FC_genes)) - sum(row.has.na) ))
  print (paste("    DEGs (padj less than alpha ) --- ", length(row.names(resSig)) ))
  print (paste("    upreg   (DEGs with log2FC > threshold) --- ", length(row.names(resSig_pos)) ))
  print (paste("    downreg (DEGs with log2FC < threshold) --- ", length(row.names(resSig_neg)) ))
  
  
  output_dir = paste("output/",
                     argExpCond,
                     sep="")
  
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  
  print ("[DEBUG] about to write the tables ")
  
  write.table(upreg_Log2FC_genes[ , c("gene_name", "ensembl","refseq", "log2FoldChange", "pvalue", "padj", 
                                      "gene_description_interpro", "gene_description_hgnc")],
              file=paste(output_dir,
                         "/",
                         "DEG_",
                         argLibrary,
                         "_",
                         argExpCond,
                         "_up.txt",
                         sep=""),
              sep="\t",
              quote=F,
              row.names=FALSE)
  write.table(downreg_Log2FC_genes[ , c("gene_name", "ensembl","refseq", "log2FoldChange", "pvalue", "padj", 
                                        "gene_description_interpro", "gene_description_hgnc")],
              file=paste(output_dir,
                         "/",
                         "DEG_",
                         argLibrary,
                         "_",
                         argExpCond,
                         "_down.txt",
                         sep=""),
              sep="\t",
              quote=F,
              row.names=FALSE)
  
  # (extra info) for the up/down-regulated tables, here is sorted by padj; most significant genes first
  upreg_Padj_genes = resSig_pos[ order( resSig_pos$padj ,decreasing=FALSE), ] 
  write.table(upreg_Padj_genes[ , c("gene_name", "ensembl","refseq", "log2FoldChange", "pvalue", "padj", 
                                    "gene_description_interpro", "gene_description_hgnc")],
              file=paste(output_dir,
                         "/",
                         "DEG_",
                         argLibrary,
                         "_",
                         argExpCond,
                         "_up.sortedByPadj.txt",
                         sep=""),
              sep="\t",
              quote=F,
              row.names=FALSE)
  downreg_Padj_genes = resSig_neg[ order( resSig_neg$padj ,decreasing=FALSE), ] 
  write.table(downreg_Padj_genes[ , c("gene_name", "ensembl","refseq", "log2FoldChange", "pvalue", "padj", 
                                      "gene_description_interpro", "gene_description_hgnc")],
              file=paste(output_dir,
                         "/",
                         "DEG_",
                         argLibrary,
                         "_",
                         argExpCond,
                         "_down.sortedByPadj.txt",
                         sep=""),
              sep="\t",
              quote=F,
              row.names=FALSE)
  
  # (extra info) full list of genes that also contains NA
  write.table(allLog2FC_genes[ , c("gene_name", "ensembl","refseq", "log2FoldChange", "pvalue", "padj", 
                                   "gene_description_interpro", "gene_description_hgnc")],
              file=paste(output_dir,
                         "/",
                         "DEG_",
                         argLibrary,
                         "_",
                         argExpCond,
                         "_allGenes.sortedByPadj.txt",
                         sep=""),
              sep="\t",
              quote=F,
              row.names=FALSE)
  
  # (extra info) full list of genes filtering out NAs
  write.table(allLog2FC_genes_filtered[ , c("gene_name", "ensembl","refseq", "log2FoldChange", "pvalue", "padj", 
                                            "gene_description_interpro", "gene_description_hgnc")],
              file=paste(output_dir,
                         "/",
                         "DEG_",
                         argLibrary,
                         "_",
                         argExpCond,
                         "_allGenesNoNA.sortedByPadj.txt",
                         sep=""),
              sep="\t",
              quote=F,
              row.names=FALSE)
  
  

  
  print (      "  [END_OF_FUNCTION_CALL] done writing tables for DEGs and generating histogram ")
}









######
#  

print ("[DEBUG] modify dds by adding other factors (refer to master excel sheet for what these are)")
#

# # there can be a better way to handle this, use switch statement?
# dds$age = ifelse(  grepl(".*_R1[012]", colnames(dds)), "3mo",
#                     ifelse(  grepl(".*_R[456]", colnames(dds)), "3mo",
#                              ifelse(  grepl(".*_R[789]", colnames(dds)), "24d",
#                                       ifelse(  grepl(".*_R[123]", colnames(dds)), "24d",
#                                                "_TODO1_"))))
# dds$age = as.factor(dds$age)
# 
# # there can be a better way to handle this, use switch statement?
# dds$genotype = ifelse(  grepl(".*_R1[012]", colnames(dds)), "wt",
#                     ifelse(  grepl(".*_R[456]", colnames(dds)), "tb",
#                              ifelse(  grepl(".*_R[789]", colnames(dds)), "wt",
#                                       ifelse(  grepl(".*_R[123]", colnames(dds)), "tb",
#                                                 "_TODO1_"))))
# dds$genotype = as.factor(dds$genotype)

# (forBatchEffect) 
#     "1" = sequencing run on 2020-08-31
#     "2" = sequencing run on 2020-11-19
#     "3" = sequencing run on 2021-04-26

dds$batch = as.factor(c(1,2,3))
dds$genotype = as.factor(c("wt","tb"))
dds$age = as.factor(c("24d","3mo"))

print ("[DEBUG] relevel; convenient way to make sure that wt is the first level in the treatment factor,   ")
# so the default log2FC calculated as treatment over control and not other way around
dds$genotype = relevel(dds$genotype, ref="wt")

i = 1
for(currentColName in colnames(dds))
{
  # print(currentColName)
  switch(currentColName,
         retina_R1 = {dds$batch[i]    = 1
                      dds$genotype[i] = "tb"
                      dds$age[i]      = "24d"},
         retina_R2 = {dds$batch[i]    = 1
                      dds$genotype[i] = "tb"
                      dds$age[i]      = "24d"},
         retina_R3 = {dds$batch[i]    = 1
                      dds$genotype[i] = "tb"
                      dds$age[i]      = "24d"},
         retina_R4 = {dds$batch[i]    = 3
                      dds$genotype[i] = "tb"
                      dds$age[i]      = "3mo"},
         retina_R5 = {dds$batch[i]    = 3
                      dds$genotype[i] = "tb"
                      dds$age[i]      = "3mo"},
         retina_R6 = {dds$batch[i]    = 3
                      dds$genotype[i] = "tb"
                      dds$age[i]      = "3mo"},
         retina_R7 = {dds$batch[i]    = 1
                      dds$genotype[i] = "wt"
                      dds$age[i]      = "24d"},
         retina_R8 = {dds$batch[i]    = 1
                      dds$genotype[i] = "wt"
                      dds$age[i]      = "24d"},
         retina_R9 = {dds$batch[i]    = 1
                      dds$genotype[i] = "wt"
                      dds$age[i]      = "24d"},
         retina_R10 = {dds$batch[i]   = 1
                      dds$genotype[i] = "wt"
                      dds$age[i]      = "3mo"},
         retina_R11 = {dds$batch[i]   = 1
                      dds$genotype[i] = "wt"
                      dds$age[i]      = "3mo"},
         retina_R12 = {dds$batch[i]   = 1
                      dds$genotype[i] = "wt"
                      dds$age[i]      = "3mo"},
         
         hypothalamus_R1 = {dds$batch[i] = 1
                         dds$genotype[i] = "tb"
                         dds$age[i]      = "24d"},
         hypothalamus_R2 = {dds$batch[i] = 2
                         dds$genotype[i] = "tb"
                         dds$age[i]      = "24d"},
         hypothalamus_R3 = {dds$batch[i] = 2
                         dds$genotype[i] = "tb"
                         dds$age[i]      = "24d"},
         hypothalamus_R4 = {dds$batch[i] = 3
                         dds$genotype[i] = "tb"
                         dds$age[i]      = "3mo"},
         hypothalamus_R5 = {dds$batch[i] = 2
                         dds$genotype[i] = "tb"
                         dds$age[i]      = "3mo"},
         hypothalamus_R6 = {dds$batch[i] = 1
                         dds$genotype[i] = "tb"
                         dds$age[i]      = "3mo"},
         hypothalamus_R7 = {dds$batch[i] = 1
                         dds$genotype[i] = "wt"
                         dds$age[i]      = "24d"},
         hypothalamus_R8 = {dds$batch[i] = 3
                         dds$genotype[i] = "wt"
                         dds$age[i]      = "24d"},
         hypothalamus_R9 = {dds$batch[i] = 2
                         dds$genotype[i] = "wt"
                         dds$age[i]      = "24d"},
         hypothalamus_R10 = {dds$batch[i] = 1
                         dds$genotype[i] = "wt"
                         dds$age[i]      = "3mo"},
         hypothalamus_R11 = {dds$batch[i] = 2
                         dds$genotype[i] = "wt"
                         dds$age[i]      = "3mo"},
         hypothalamus_R12 = {dds$batch[i] = 2
                         dds$genotype[i] = "wt"
                         dds$age[i]      = "3mo"},
  )
  i = i + 1
}


colData(dds)[ ,c("genotype","age","batch")]  # can also get tissue which is called "condition"



######

print ("[DEBUG] subsetting dds to different variables ")
#  
#

## (by all samples) make a copy of dds
dds_all_samples = dds
dds_all_samples$condition = "all"          # note modifying the $condition and $age for downstream print statements to make sense
dds_all_samples$age = "samples" 
dim(dds_all_samples)

## (by tissue) two groups of 12 samples each
dds_retina = dds[,dds$condition == "retina"]
dds_retina$condition = "retina"          # note modifying the $condition and $age for downstream print statements to make sense
dds_retina$age = "only" 
dim(dds_retina)
dds_hypothalamus = dds[,dds$condition == "hypothalamus"]
dds_hypothalamus$condition = "hypothalamus"          # note modifying the $condition and $age for downstream print statements to make sense
dds_hypothalamus$age = "only" 
dim(dds_hypothalamus)


## (by experimental conditions) 4 groups of 6 samples each; four comparisons made between wt (control) vs tb (treatment)
dds_retina_24d = dds[, (dds$condition == "retina") & (dds$age == "24d")]
dim(dds_retina_24d)
dds_retina_3mo = dds[, (dds$condition == "retina") & (dds$age == "3mo")]
dim(dds_retina_3mo)
dds_hypothalamus_24d = dds[, (dds$condition == "hypothalamus") & (dds$age == "24d")]
dim(dds_hypothalamus_24d)
dds_hypothalamus_3mo = dds[, (dds$condition == "hypothalamus") & (dds$age == "3mo")]
dim(dds_hypothalamus_3mo)




as.data.frame( colData(dds_all_samples))
as.data.frame( colData(dds_retina))
as.data.frame( colData(dds_hypothalamus))
as.data.frame( colData(dds_retina_3mo))
as.data.frame( colData(dds_hypothalamus_24d))
as.data.frame( colData(dds_hypothalamus_3mo))


print ("[DEBUG] write a table of normalized counts for retina 3 month that can be used for GSEA ")
normalized_counts <- counts(dds_retina_3mo, normalized= TRUE)   # need to reorder it to make it easier for layer (next line)
head(normalized_counts)[,]  
write.table(
  data.frame("geneid"=rownames(normalized_counts),normalized_counts), 
  file="dds_retina_3mo_normalized_counts.txt", 
  sep="\t", 
  quote=F, 
  row.names=FALSE)



#####
#  Running the DESEQ2 pipeline 
####

print ("[DEBUG] change design to look at it based on genotype; rerun DESeq ")



# (notice) for the design we need to account for "batch effect", which is different sequencing runs
#          this is accounted for as a covariate by adding it to the design formula   "design = ~ batch + condition"
#          (src) https://www.biostars.org/p/403053/
#          (about error) https://bioconductor.org/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#model-matrix-not-full-rank
design(dds_all_samples) <- formula(~ batch + genotype)
dds_all_samples <- DESeq(dds_all_samples)
design(dds_retina) <- formula(~ genotype)    # remove batch correction bc "full model matrix is less than full rank" error appears ()
dds_retina <- DESeq(dds_retina)
design(dds_hypothalamus) <- formula(~ batch + genotype)
dds_hypothalamus <- DESeq(dds_hypothalamus)

design(dds_retina_24d) <- formula(~ genotype) # remove batch correction bc "full model matrix is less than full rank" error appears ()
dds_retina_24d <- DESeq(dds_retina_24d)
design(dds_retina_3mo) <- formula(~ genotype) # remove batch correction bc "full model matrix is less than full rank" error appears ()
dds_retina_3mo <- DESeq(dds_retina_3mo)
design(dds_hypothalamus_24d) <- formula(~ batch + genotype)
dds_hypothalamus_24d <- DESeq(dds_hypothalamus_24d)
design(dds_hypothalamus_3mo) <- formula(~ batch + genotype)
dds_hypothalamus_3mo <- DESeq(dds_hypothalamus_3mo)



#(update) not necessary as DESeq() calls this; 
# refer to page 7 of DESeq2 manual: https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf
# ## perform "between group" normalization by calling estimateSizeFactors of DESeq2 library ; 
# dds_retina_24d = estimateSizeFactors(dds_retina_24d)
# dds_retina_3mo = estimateSizeFactors(dds_retina_3mo)
# dds_hypothalamus_24d = estimateSizeFactors(dds_hypothalamus_24d)
# dds_hypothalamus_3mo = estimateSizeFactors(dds_hypothalamus_3mo)


print ("[DEBUG] before FOR LOOP starts")
# listOfExpCond = list( dds_retina_24d, dds_retina_3mo, dds_hypothalamus_24d, dds_hypothalamus_3mo, dds_all_samples, dds_retina, dds_hypothalamus)
# listOfExpCond = list( dds_hypothalamus_24d, dds_hypothalamus_3mo, dds_retina_24d, dds_retina_3mo)
listOfExpCond = list( dds_retina_3mo)
for(dds_currentExpCond in listOfExpCond)
{
  print ("")
  print (paste("== IN LOOP performing DESeq DEG analysis pipeline for ----- ", dds_currentExpCond$condition[1],"_", dds_currentExpCond$age[1]))

  
  myExperimentalCondition = paste(dds_currentExpCond$condition[1],dds_currentExpCond$age[1], sep="_")
  myOutputDir = paste("output/",
                      myExperimentalCondition,
                      sep="")
  
  ####################
  ## now to use ----->   DESEQ2
  ####################
  
  print ("")
  print ("== + generating DESEQ2 based data...")
  print ("[FUNCTION CALL TO MY CUSTOM] useThisLibraryToGetDEGsFromInputList")
  libraryUsedForFilename = "deseq2"
  res = results(dds_currentExpCond)
  myInputList = res
  useThisLibraryToGetDEGsFromInputList (myInputList, libraryUsedForFilename, myExperimentalCondition)
  print ("")
  
  # debug stuff 
  sum (res$pvalue < 0.01, na.rm=TRUE)
  sum (res$pvalue < 0.05, na.rm=TRUE)
  length (res [ which(res$padj <0.10 ), ] )
  length (res [ which(res$padj <0.05 ), ] )
  res
  table (is.na (res$pvalue))
  
  
  
  
  ####################
  ## now to use   --->    edgeR
  ####################
  
  print ("")
  print ("== + generating edgeR based data...")
  
  #(src) https://github.com/markziemann/de_compare/blob/main/deg_compare.Rmd
  
  # convert data needed from 'dds' object and prepare it for edgeR
  
  #   get the data into another format
  numericData = counts(dds_currentExpCond)
  nonzeroData = numericData[rowSums(numericData==0, na.rm=TRUE)<ncol(numericData), ]
  
  #   need to rename the columns to C1,C2,C3,T1,T2,T3
  newColumnNames = vector("list")
  i = 1
  i_C = 1
  i_T = 1
  
  for(current_genotype in dds_currentExpCond$genotype)
  {
    newColumnName = ""
    if(current_genotype=="wt")
    {
      newColumnName = paste("C",i_C, sep="")
      i_C = i_C + 1
    }
    else
    {
      newColumnName = paste("T",i_T, sep="")
      i_T = i_T + 1
    }
    newColumnNames[[i]] = newColumnName
    # need to add the occurence # too
    
    i = i + 1
  }
  colnames(nonzeroData) = newColumnNames
  
  # setup dataset counts
  xx <- nonzeroData[which(rowMeans(nonzeroData)>-1),]
  str(xx)
  dim(xx)
  
  # setup dataset header
  ss <- as.data.frame(colnames(xx))
  ss$trt <- as.integer(grepl("T",ss[,1]))
  rownames(ss) <- ss[,1]
  ss[,1]=NULL
  design <- model.matrix(~ss$trt)
  rownames(design) <- rownames(ss)
  
  # DEG via edgeR using quasi-likelihood (QL) F-test
  z <- DGEList(counts=xx)
  z <- calcNormFactors(z)
  z <- estimateDisp(z, design,robust=TRUE,prior.df=1)
  fit <- glmQLFit(z, design)
  ql <- glmQLFTest(fit)
  
  dge<-as.data.frame(topTags(ql,n=Inf))
  colnames(dge) = c("log2FoldChange", "logCPM", "F", "pvalue", "padj")    # need to rename the columwns so they are like DESEQ2's format
  # dge$dispersion<-ql$dispersion
  # dge<-merge(dge,ql$fitted.values,by='row.names')
  # rownames(dge)=dge$Row.names
  # dge$Row.names=NULL
  # dge<-dge[order(dge$PValue),]
  # head(dge,10)
  # dge_edgerql <- dge
  # sig <- subset(dge_edgerql,FDR<P_ALPHA_CUTOFF)
  # dge_edgerql_up <- rownames(subset(sig,logFC> P_LOG2FC_THRESHOLD ))
  # dge_edgerql_dn <- rownames(subset(sig,logFC< (-1*P_LOG2FC_THRESHOLD)))
  # length(dge_edgerql_up)
  # length(dge_edgerql_dn)
  
  
  print ("[FUNCTION CALL TO MY CUSTOM] useThisLibraryToGetDEGsFromInputList")
  libraryUsedForFilename = "edgeR_"
  myInputList = dge
  useThisLibraryToGetDEGsFromInputList (myInputList, libraryUsedForFilename, myExperimentalCondition)
  print ("")
  
  
  
  ####################
  ## now to plot   --->    venn diagram of deseq2 vs edgeR up/down-regulated genes
  ####################
  print ("== + check if can generate venn diagram...")
  
  deg_deseq2 = res
  deg_edger  = dge
  
  print ("[DEBUG] subset to get top down/up-regulated genes")
  resSig_deseq2 = deg_deseq2 [ which(deg_deseq2$padj < P_ALPHA_CUTOFF ), ]             # fdr correction
  resSig_edger = deg_edger [ which(deg_edger$padj < P_ALPHA_CUTOFF ), ]             # fdr correction
  dim(resSig_deseq2)
  dim(resSig_edger)
  
  
  print ("[DEBUG] subset again to get down/up-regulated genes that are above the threshold")
  resSig_pos_deseq2 = resSig_deseq2 [ which(resSig_deseq2$log2FoldChange > P_LOG2FC_THRESHOLD ), ]             # upregulated DEGs have log2FC > 1.5
  resSig_pos_edger = resSig_edger [ which(resSig_edger$log2FoldChange > P_LOG2FC_THRESHOLD ), ]             # upregulated DEGs have log2FC > 1.5
  resSig_neg_deseq2 = resSig_deseq2 [ which(resSig_deseq2$log2FoldChange < (-1*P_LOG2FC_THRESHOLD) ), ]        # downregulated DEGs have log2FC < -1.5
  resSig_neg_edger = resSig_edger [ which(resSig_edger$log2FoldChange < (-1*P_LOG2FC_THRESHOLD) ), ]        # downregulated DEGs have log2FC < -1.5
  
  print ("[DEBUG] save info about all genes to some tables")
  highestLog2FC_genes_deseq2 = resSig_pos_deseq2[ order( resSig_pos_deseq2$log2FoldChange ,decreasing=TRUE), ] 
  highestLog2FC_genes_edger = resSig_pos_edger[ order( resSig_pos_edger$log2FoldChange ,decreasing=TRUE), ] 
  lowestLog2FC_genes_deseq2 = resSig_neg_deseq2[ order( resSig_neg_deseq2$log2FoldChange ), ]
  lowestLog2FC_genes_edger = resSig_neg_edger[ order( resSig_neg_edger$log2FoldChange ), ]
  
  genenames_deseq2_up = rownames(highestLog2FC_genes_deseq2)
  genenames_edger_up = rownames(highestLog2FC_genes_edger)
  genenames_deseq2_dn = rownames(lowestLog2FC_genes_deseq2)
  genenames_edger_dn = rownames(lowestLog2FC_genes_edger)
  
  # only proceed if all four tables have rows, otherwise program will crash
  if (length(genenames_deseq2_up) > 0 &&
      length(genenames_edger_up) > 0 &&
      length(genenames_deseq2_dn) > 0 &&
      length(genenames_edger_dn) > 0) {
    
    print ("== + generating venn diagram...")
    
    v1 <- list("edgeR up"=genenames_edger_up, 
               "edgeR dn"=genenames_edger_dn,
               "DESeq2 up"=genenames_deseq2_up,
               "DESeq2 dn"=genenames_deseq2_dn)
    
    pdf(paste(myOutputDir,
              "/DEG_venndiagram_",
              myExperimentalCondition,
              "_deseq2VSedger.pdf",
              sep=""),
        width=10, 
        height=8
    )
    p <- plot(euler(v1),quantities = TRUE, adjust_labels=TRUE)
    print(p)
    dev.off()
    
    
  } else {
    print ("== + NOT ABLE TO create venn diagram because there is at least one table with zero up/down-regulated DEGs...")
    
  }
  
  
  
  ####################
  ## now to find   --->    common genes of deseq2 vs edgeR up/down-regulated genes
  ####################
  print ("== + check if can determine the common genes between deseq2 and edgerR output...")
  
  # only proceed if all four tables have rows, otherwise program will crash
  if (length(genenames_deseq2_up) > 0 &&
      length(genenames_edger_up) > 0 &&
      length(genenames_deseq2_dn) > 0 &&
      length(genenames_edger_dn) > 0) {
    
    print ("== + determining common genes ...")
    
    
    print ("[DEBUG] subset the current data so the columns match exactly")
    deseq2_up = subset(highestLog2FC_genes_deseq2, select=
                         c("log2FoldChange", "pvalue", "padj"))
    deseq2_dn = subset(lowestLog2FC_genes_deseq2, select=
                         c("log2FoldChange", "pvalue", "padj"))
    edger_up = subset(highestLog2FC_genes_edger, select=
                        c("log2FoldChange", "pvalue", "padj"))
    edger_dn = subset(lowestLog2FC_genes_edger, select=
                        c("log2FoldChange", "pvalue", "padj"))
    
    
    print ("[DEBUG] make sure the data is dataframes ")
    deseq2_up = as.data.frame(deseq2_up)
    deseq2_dn = as.data.frame(deseq2_dn)
    edger_up = as.data.frame(edger_up)
    edger_dn = as.data.frame(edger_dn)
    
    
    print ("[DEBUG] find the intersection of genes among the two dataframes")
    #  (notice) not a perfect merge, since for each set the values of log2fc, pval, and padj slightly different
    intersection_up = merge(deseq2_up, edger_up, by="row.names") 
    row.names(intersection_up) = intersection_up$Row.names   # need to reset rownames after merge and get rid of new column made
    intersection_up = within(intersection_up, rm(Row.names))
    colnames(intersection_up) = c("log2FoldChange.deseq2", "pvalue.deseq2", "padj.deseq2", "log2FoldChange.edger", "pvalue.edger", "padj.edger")
    intersection_up = intersection_up[ order( intersection_up$log2FoldChange.deseq2 ,decreasing=TRUE), ] 

    intersection_dn = merge(deseq2_dn, edger_dn, by="row.names")   
    row.names(intersection_dn) = intersection_dn$Row.names   # need to reset rownames after merge and get rid of new column made
    intersection_dn = within(intersection_dn, rm(Row.names))
    colnames(intersection_dn) = c("log2FoldChange.deseq2", "pvalue.deseq2", "padj.deseq2", "log2FoldChange.edger", "pvalue.edger", "padj.edger")
    intersection_dn = intersection_dn[ order( intersection_dn$log2FoldChange.deseq2 ), ]
  
    
    print ("[DEBUG] using BIOMART library to get gene symbol and descriptions")
    intersection_up$ensembl <- sapply( strsplit( rownames(intersection_up), split="\\+" ), "[", 1)
    intersection_dn$ensembl <- sapply( strsplit( rownames(intersection_dn), split="\\+" ), "[", 1)
    genemap <- getBM(  attributes = c("ensembl_gene_id","external_gene_name","interpro_description", "description", "refseq_mrna"),
                       filters = "ensembl_gene_id",
                       values = intersection_up$ensembl,
                       mart = ENSEMBL_BIOMART)
    idx <- match( intersection_up$ensembl, genemap$ensembl_gene_id )
    intersection_up$gene_name <- genemap$external_gene_name[ idx ]
    intersection_up$gene_description_interpro <- genemap$interpro_description[ idx ]
    intersection_up$gene_description_hgnc <- genemap$description[ idx ]
    intersection_up$refseq <- genemap$refseq_mrna[ idx ]
    # intersection_up$pw_mgi <- genemap$mgi_id[ idx ]
    # intersection_up$pw_kegg <- genemap$kegg_enzyme[ idx ]
    # intersection_up$pw_entrez <- genemap$entrezgene_id[ idx ]
    head(intersection_up, 4)
                                                                                      # (notice) getting more for downstream GO/GSEA/KEGG pathway analyses
    genemap <- getBM(  attributes = c("ensembl_gene_id","external_gene_name","interpro_description", "description", "refseq_mrna"),
                       filters = "ensembl_gene_id",
                       values = intersection_dn$ensembl,
                       mart = ENSEMBL_BIOMART)
    idx <- match( intersection_dn$ensembl, genemap$ensembl_gene_id )
    intersection_dn$gene_name <- genemap$external_gene_name[ idx ]
    intersection_dn$gene_description_interpro <- genemap$interpro_description[ idx ]
    intersection_dn$gene_description_hgnc <- genemap$description[ idx ]
    intersection_dn$refseq <- genemap$refseq_mrna[ idx ]
    # intersection_dn$pw_mgi <- genemap$mgi_id[ idx ]
    # intersection_dn$pw_kegg <- genemap$kegg_enzyme[ idx ]
    # intersection_dn$pw_entrez <- genemap$entrezgene_id[ idx ]
    head(intersection_dn, 4)
    
    
    print ("[DEBUG] filter out the genes where gene symbol is NA; could still have an ENSEMBL id but be pseudogene or ncrna")
    row.has.na <- apply(intersection_up, 1, function(x){any(is.na(x))})
    intersection_up = intersection_up[!row.has.na,] 
    row.has.na <- apply(intersection_dn, 1, function(x){any(is.na(x))})
    intersection_dn = intersection_dn[!row.has.na,] 

    
    print ("[DEBUG] print the significant genes and the list of non-NA common genes sorted by log2FC")
    
    
    print (paste("    (deseq2) signif. genes --- ", length(row.names(resSig_deseq2)) ))
    print (paste("    (deseq2)   log2FC > 0  --- ", length(row.names(subset(resSig_deseq2, log2FoldChange > 0)))))
    print (paste("    (deseq2)   log2FC < 0  --- ", length(row.names(subset(resSig_deseq2, log2FoldChange < 0)))))
    print (paste("    (edger)  signif. genes --- ", length(row.names(resSig_edger)) ))
    print (paste("    (deseq2) upreg DEGs    --- ", length(row.names(deseq2_up)) ))
    print (paste("    (edger_) upreg DEGs    --- ", length(row.names(edger_up)) ))
    print (paste("    (common) upreg DEGs    --- ", length(row.names(intersection_up)) ))
    
    print (paste("    (deseq2) dnreg DEGs    --- ", length(row.names(deseq2_dn)) ))
    print (paste("    (edger_) dnreg DEGs    --- ", length(row.names(edger_dn)) ))
    print (paste("    (common) dnreg DEGs    --- ", length(row.names(intersection_dn)) ))
    
    
    intersection_all = rbind(intersection_up, intersection_dn)
    intersection_all = intersection_all[ order( intersection_all$log2FoldChange.deseq2 ,decreasing=TRUE), ] 
    print (paste("    (common) all reg DEGs --- ", length(row.names(intersection_all)) ))
    
    write.table(intersection_up[ , c("gene_name", "ensembl","refseq", "log2FoldChange.deseq2", 
                                     "pvalue.deseq2", "padj.deseq2", "log2FoldChange.edger", "pvalue.edger", "padj.edger", 
                                     "gene_description_interpro", "gene_description_hgnc")],
                file=paste(myOutputDir,
                           "/",
                           "DEG_",
                           "commonUpRegulatedGenes",
                           "_",
                           myExperimentalCondition,
                           ".txt",
                           sep=""),
                sep="\t",
                quote=F,
                row.names=FALSE)
    write.table(intersection_dn[ , c("gene_name", "ensembl","refseq", "log2FoldChange.deseq2", 
                                     "pvalue.deseq2", "padj.deseq2", "log2FoldChange.edger", "pvalue.edger", "padj.edger", 
                                     "gene_description_interpro", "gene_description_hgnc")],
                file=paste(myOutputDir,
                           "/",
                           "DEG_",
                           "commonDownRegulatedGenes",
                           "_",
                           myExperimentalCondition,
                           ".txt",
                           sep=""),
                sep="\t",
                quote=F,
                row.names=FALSE)
    write.table(intersection_all[ , c("gene_name", "ensembl","refseq", "log2FoldChange.deseq2", 
                                     "pvalue.deseq2", "padj.deseq2", "log2FoldChange.edger", "pvalue.edger", "padj.edger", 
                                     "gene_description_interpro", "gene_description_hgnc")],
                file=paste(myOutputDir,
                           "/",
                           "DEG_",
                           "commonAllDiffentiallyRegulatedGenes",
                           "_",
                           myExperimentalCondition,
                           ".txt",
                           sep=""),
                sep="\t",
                quote=F,
                row.names=FALSE)
    
    # (extra info) for the up/down-regulated tables, here is sorted by padj; most significant genes first
    upreg_Padj_genes = intersection_up[ order( intersection_up$padj.deseq2 ,decreasing=FALSE), ] 
    write.table(upreg_Padj_genes[ , c("gene_name", "ensembl","refseq", "log2FoldChange.deseq2", 
                                      "pvalue.deseq2", "padj.deseq2", "log2FoldChange.edger", "pvalue.edger", "padj.edger",
                                      "gene_description_interpro", "gene_description_hgnc")],
                file=paste(myOutputDir,
                           "/",
                           "DEG_",
                           "commonUpRegulatedGenes",
                           "_",
                           myExperimentalCondition,
                           ".sortedByPadj.txt",
                           sep=""),
                sep="\t",
                quote=F,
                row.names=FALSE)
    downreg_Padj_genes = intersection_dn[ order( intersection_dn$padj.deseq2 ,decreasing=FALSE), ] 
    write.table(downreg_Padj_genes[ , c("gene_name", "ensembl","refseq", "log2FoldChange.deseq2", 
                                        "pvalue.deseq2", "padj.deseq2", "log2FoldChange.edger", "pvalue.edger", "padj.edger", 
                                        "gene_description_interpro", "gene_description_hgnc")],
                file=paste(myOutputDir,
                           "/",
                           "DEG_",
                           "commonDownRegulatedGenes",
                           "_",
                           myExperimentalCondition,
                           ".sortedByPadj.txt",
                           sep=""),
                sep="\t",
                quote=F,
                row.names=FALSE)
    
    # (extra info) histogram to show where the genes lie with regards to log2fc
    pdf(paste(myOutputDir,
              "/DEG_histogram_",
              myExperimentalCondition,
              ".pdf",
              sep=""),
        width=10, 
        height=6
    )
    forHistogram = subset(resSig_deseq2, log2FoldChange > -7 & log2FoldChange < 7)
    hist( forHistogram$log2FoldChange, 
          breaks=c(-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7), 
          col="Gray",
          axes=FALSE,
          ylim = c(0,1500),
          xlim = c(-7,7),
          xlab = "Log2 (Fold Change)",
          ylab = "Number of Significant Genes",
          main = "Histogram of Log2FC for Significant Genes")
    axis(1, at=c(-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7), 
         labels=c("-7","-6","-5","-4","-3","-2","-1","0","1","2","3","4","5","6","7"))
    axis(2, at=c(0, 500, 1000, 1500), 
         labels=c("0", "500", "1000", "1500"))
    abline(h=0)
    
    dev.off()

    
    # pdf(paste(myOutputDir,
    #           "/DEG_histogram_",
    #           myExperimentalCondition,
    #           ".pdf",
    #           sep=""),
    #     width=10, 
    #     height=6
    # )
    # hist( resSig_deseq2$log2FoldChange, 
    #       breaks=c(-12, -10,-8,-6,-4,-2,0,2,4,6,8,10,12), 
    #       col="Gray",
    #       axes=FALSE,
    #       ylim = c(0,2000),
    #       xlim = c(-12,12),
    #       xlab = "Log2 (Fold Change)",
    #       ylab = "Number of Significant Genes",
    #       main = "Histogram of Log2FC for Significant Genes")
    # axis(1, at=c(-12, -10,-8,-6,-4,-2,0,2,4,6,8,10,12), 
    #      labels=c("-12", "-10","-8","-6","-4","-2","0","2","4","6","8","10","12"))
    # axis(2, at=c(0, 500, 1000, 1500,2000), 
    #      labels=c("0", "500", "1000", "1500","2000"))
    # abline(h=0)
    # 
    # dev.off()
    
    
    
    # (extra analysis) cross-ref with the HOMER-outputted genes that have ERRE
    print ("[DEBUG] match up the HOMER-based genes that were found to have erre sequence")
    erre_no_mistmatches <- read.table("outputfile_genesWithMotifFor_erre__zeroMismatches.txt", sep="\t", header=TRUE)
    row.names(erre_no_mistmatches) = erre_no_mistmatches$Ensembl
    
    intersection_all_erre = merge(intersection_all, erre_no_mistmatches, by="row.names") 
    write.table(intersection_all_erre[ , c("gene_name","ensembl","refseq", "Motif.Name","Sequence", "Strand", "MotifScore", 
                                        "log2FoldChange.deseq2", "pvalue.deseq2", "padj.deseq2",
                                        "gene_description_interpro", "gene_description_hgnc")],
                file=paste(myOutputDir,
                           "/",
                           "HOMER_",
                           "commonDEGsWithERREmotif",
                           "_",
                           myExperimentalCondition,
                           ".txt",
                           sep=""),
                sep="\t",
                quote=F,
                row.names=FALSE)
    
  } else {
    print ("== + NOT ABLE TO determine common genes because there is at least one table with zero up/down-regulated DEGs...")
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ####################
  ## now to normalize data ---> using rlog for non-DEG analysis; "normalize within"
  ####################
  
  # some PCA here only looking at 6 samples
    # (notice) here we can work with transformed/normalized counts
    # (notice) for DEG analysis above, must use raw counts ONLY
  dds_rlog = assay(rlog(dds_currentExpCond,blind=FALSE))   # gets the counts of rlog ---> to pca
  numericData = dds_rlog
  nonzeroData = numericData[rowSums(numericData==0, na.rm=TRUE)<ncol(numericData), ]
  dim(numericData)
  nonzeroDataTranspose = t(nonzeroData)
  myPCA = prcomp(nonzeroDataTranspose, center=T, scale=T)
  myPCA_scale_df = as.data.frame(myPCA["scale"])
  myPCA_rotation_df = as.data.frame(myPCA["rotation"])
  head(myPCA_rotation_df,1) 
  
  
  ### save this data to a table; 
  #   (NOTE - this has modification to allow geneid in topleft corner for dt analysis)
  write.table(
    data.frame("geneid"=rownames(nonzeroData),format(nonzeroData,digits=4,scientific = FALSE)), 
    file=paste(myOutputDir,
               "/PCA_",
               myExperimentalCondition,
               "_RLOGDATA_usedForPlots.txt",
               sep=""),
    sep="\t", 
    quote=F, 
    row.names=FALSE)
  
  
  
  ####################
  ## now to create   --->    bar plot
  ####################
  
  
  # What percentage of the total variance is included in the first two principal components? 
  # Make a barplot showing the variance in each component (aka "Scree Plot").
  print ("== + generating bar plot (scree plot)...")
  pdf(paste(myOutputDir,
                 "/PCA_",
                  myExperimentalCondition,
                 "_SCREEPLOT_barplotOfEachPCVariance.pdf",
                  sep=""),
      width=10, 
      height=6
      )
  myBarplot = barplot(summary(myPCA)[[6]][c(2,5,8,11,14,17,20,23,26,29)],
                      ylim=c(0,1.0),
                      ylab="Variance due to Principal Component",
                      xlab=c("Principle Component (Total Variance)"),
                      names.arg = c(gsub(" ", "", paste("PC1 \n (", summary(myPCA)[[6]][3],")")),
                                    gsub(" ", "", paste("PC2 \n (", summary(myPCA)[[6]][6],")")),
                                    gsub(" ", "", paste("PC3 \n (", summary(myPCA)[[6]][9],")")), 
                                    gsub(" ", "", paste("PC4 \n (", summary(myPCA)[[6]][12],")")),
                                    gsub(" ", "", paste("PC5 \n (", summary(myPCA)[[6]][15],")")),
                                    gsub(" ", "", paste("PC6 \n (", summary(myPCA)[[6]][18],")")),
                                    gsub(" ", "", paste("PC7 \n (", summary(myPCA)[[6]][21],")")),
                                    gsub(" ", "", paste("PC8 \n (", summary(myPCA)[[6]][24],")")),
                                    gsub(" ", "", paste("PC9 \n (",  summary(myPCA)[[6]][27],")")),
                                    gsub(" ", "", paste("PC10 \n (", summary(myPCA)[[6]][30],")"))
                      ), 
                      col="red",
                      main="Barplot of Proportion Variance due to each Principle Component"
  )
  abline(h=0)
  
  ## Add text at top of bars
  individualProportions = c(summary(myPCA)[[6]][2],
                            summary(myPCA)[[6]][5],
                            summary(myPCA)[[6]][8],
                            summary(myPCA)[[6]][11],
                            summary(myPCA)[[6]][14],
                            summary(myPCA)[[6]][17],
                            summary(myPCA)[[6]][20],
                            summary(myPCA)[[6]][23],
                            summary(myPCA)[[6]][26],
                            summary(myPCA)[[6]][29]
  )
  text(x = myBarplot, y = individualProportions, label = individualProportions, pos = 3, cex = 0.8, col = "red")
 
  dev.off()
  
  
  ####################
  ## now to create   --->    MA plot 
  ####################

  print ("== + generating MA plot...")
  
  # some additional diagonostic plots using DESeq2 library to quickly make
  # MA plot gives you idea of the average expression between the compared groups genome-wide
        # (defn of x-axis) It indicates the average expression between the groups you
        #   compared, which is simply the mean of the normalized counts of the samples
        #   used in the comparison group. The "meaning" is that genes with higher counts
        #   are more trustworthy (=have higher statistical power) than low count genes.
        #   Biologically this is because either the gene is lowly expressed and/or short
        #   and therefore intrinsically gets lower counts than a longer one at equal
        #   expression level.
  
  pdf(paste(myOutputDir,
            "/MAPLOT_",
            myExperimentalCondition,
            "_diagnosticPlot.pdf",
            sep=""),)
  p <- plotMA (res,
                ylim =c(-10,10 ),
                
                main=paste("MA plot of",
                           dds_currentExpCond$condition[1],
                           "at",
                           dds_currentExpCond$age[1],
                           "samples",
                           sep=" ")
  )

  dev.off()
  
  
  # ####################
  # ## now to create   --->    Volcano plot 
  # ####################
  # 
  # print ("== + generating Volcano plot...")
  # 
  # # some additional diagonostic plots using DESeq2 library to quickly make
  # # Volcano plot ___TODO___what_does_it_tell_us__
  # 
  # pdf(paste(myOutputDir,
  #           "/VOLCANOPLOT_",
  #           myExperimentalCondition,
  #           "_diagnosticPlot.pdf",
  #           sep=""),
  #     width=10, 
  #     height=6)
  # p <- plotMA (res,
  #              ylim =c(-10,10 ),
  #              
  #              main=paste("Volcano plot of",
  #                         dds_currentExpCond$condition[1],
  #                         "at",
  #                         dds_currentExpCond$age[1],
  #                         "samples",
  #                         sep=" ")
  # )
  # 
  # dev.off()
  # 
  
  ####################
  ## now to create   --->    PCA plots
  ####################

  
  print ("== + generating PCA plot...")
  pc <- as.data.frame(myPCA$x)
  
  # for multifactorial experimental design   (BY-COLOR)
  #     can use different sample types based on mutation and age
  pc$experimental_conditions <- ifelse(  grepl("retina_R1[012]", row.names(pc)), "wt_3mo",
                                         ifelse(  grepl("retina_R[456]", row.names(pc)), "tb_3mo",
                                                  ifelse(  grepl("retina_R[789]", row.names(pc)), "wt_24d",
                                                           ifelse(  grepl("retina_R[123]", row.names(pc)), "tb_24d",
                                                                    
                                                                    
                                                                    ifelse(  grepl("hypothalamus_R1[012]", row.names(pc)), "wt_3mo",
                                                                             ifelse(  grepl("hypothalamus_R[456]", row.names(pc)), "tb_3mo",
                                                                                      ifelse(  grepl("hypothalamus_R[789]", row.names(pc)), "wt_24d",
                                                                                               ifelse(  grepl("hypothalamus_R[123]", row.names(pc)), "tb_24d",
                                                                                                        "_TODO1_"))))))))
  # for multifactorial experimental design   (BY-SHAPE)
  #     can use hypothalamus vs retina here
  pc$tissue <- ifelse( grepl("hypothalamus", row.names(pc), fixed=TRUE), 
                       "hypothalamus",
                       "retina")
  
  ##   PCA Plot for PC1 vs PC2
  library(ggplot2)
  pc1_percent = summary(myPCA)[[6]][2]*100
  pc2_percent = summary(myPCA)[[6]][5]*100
  p <- ggplot(  pc,
                aes(x=PC1,y=PC2, color = experimental_conditions, shape = tissue)) +   #notice adding onto object
                geom_point(size = 4) +
    # scale_shape_manual(values=c(17)) +                                 ## (RETINA only as triangle)
    # guides(colour = guide_legend(override.aes = list(shape = 17))) +   ## (RETINA only as triangle)
    
    
    # plot labels
    xlab(paste0("PC1 (",round(pc1_percent,1),"%", ")" )) +
    ylab(paste0("PC2 (",round(pc2_percent,1),"%", ")" )) +
    labs(title = "PCA (PC1 vs PC2) of rlog-transformed data") +
    
    # custom colors
    scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#56B4E9"))
  
  print(p)
  
  ggsave(filename=paste(myOutputDir,
                    "/PCA_",
                    myExperimentalCondition,
                    "_PC1vsPC2.pdf",
                    sep=""),
         plot=p,
         width = 5, height = 4, dpi = 300, units = "in", device='pdf')
  
  
  
  ##   PCA Plot for PC3 vs PC4
  pc3_percent = summary(myPCA)[[6]][8]*100
  pc4_percent = summary(myPCA)[[6]][11]*100
  p <- ggplot(  pc,
                aes(x=PC3,y=PC4, color = experimental_conditions, shape = tissue)) +   #notice adding onto object
    geom_point(size = 4) +
    # scale_shape_manual(values=c(17)) +                                 ## (RETINA only as triangle)
    # guides(colour = guide_legend(override.aes = list(shape = 17))) +   ## (RETINA only as triangle)
    
    
    # plot labels
    xlab(paste0("PC3 (",round(pc3_percent,1),"%", ")" )) +
    ylab(paste0("PC4 (",round(pc4_percent,1),"%", ")" )) +
    labs(title = "PCA (PC3 vs PC4) of rlog-transformed data") +
    
    # custom colors
    scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#56B4E9"))
  
  print(p)
  
  ggsave(filename=paste(myOutputDir,
                        "/PCA_",
                        myExperimentalCondition,
                        "_PC3vsPC4.pdf",
                        sep=""),
         plot=p,
         width = 5, height = 4, dpi = 300, units = "in", device='pdf')
  
  ##   PCA Plot for PC3 vs PC5
  pc3_percent = summary(myPCA)[[6]][8]*100
  pc5_percent = summary(myPCA)[[6]][14]*100
  p <- ggplot(  pc,
                aes(x=PC3,y=PC5, color = experimental_conditions, shape = tissue)) +   #notice adding onto object
    geom_point(size = 4) +
    # scale_shape_manual(values=c(17)) +                                 ## (RETINA only as triangle)
    # guides(colour = guide_legend(override.aes = list(shape = 17))) +   ## (RETINA only as triangle)
    
    
    # plot labels
    xlab(paste0("PC3 (",round(pc3_percent,1),"%", ")" )) +
    ylab(paste0("PC5 (",round(pc5_percent,1),"%", ")" )) +
    labs(title = "PCA (PC3 vs PC5) of rlog-transformed data") +
    
    # custom colors
    scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#56B4E9"))
  
  print(p)
  
  ggsave(filename=paste(myOutputDir,
                        "/PCA_",
                        myExperimentalCondition,
                        "_PC3vsPC5.pdf",
                        sep=""),
         plot=p,
         width = 5, height = 4, dpi = 300, units = "in", device='pdf')
         
  ## ------ these are mainly for all_samples, retina only, hypothalamus only, not the exp conditions
  
  ##   PCA Plot for PC5 vs PC6
  pc5_percent = summary(myPCA)[[6]][14]*100
  pc6_percent = summary(myPCA)[[6]][17]*100
  
  if (pc6_percent > 0)
  {
    
    p <- ggplot(  pc,
                  aes(x=PC5,y=PC6, color = experimental_conditions, shape = tissue)) +   #notice adding onto object
      geom_point(size = 4) +
      # scale_shape_manual(values=c(17)) +                                 ## (RETINA only as triangle)
      # guides(colour = guide_legend(override.aes = list(shape = 17))) +   ## (RETINA only as triangle)
      
      
      # plot labels
      xlab(paste0("PC5 (",round(pc5_percent,1),"%", ")" )) +
      ylab(paste0("PC6 (",round(pc6_percent,1),"%", ")" )) +
      labs(title = "PCA (PC5 vs PC6) of rlog-transformed data") +
      
      # custom colors
      scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#56B4E9"))
    
    print(p)
    
    ggsave(filename=paste(myOutputDir,
                          "/PCA_",
                          myExperimentalCondition,
                          "_PC5vsPC6.pdf",
                          sep=""),
           plot=p,
           width = 5, height = 4, dpi = 300, units = "in", device='pdf')
    
    ##   PCA Plot for PC7 vs PC8
    pc7_percent = summary(myPCA)[[6]][20]*100
    pc8_percent = summary(myPCA)[[6]][23]*100
    p <- ggplot(  pc,
                  aes(x=PC7,y=PC8, color = experimental_conditions, shape = tissue)) +   #notice adding onto object
      geom_point(size = 4) +
      # scale_shape_manual(values=c(17)) +                                 ## (RETINA only as triangle)
      # guides(colour = guide_legend(override.aes = list(shape = 17))) +   ## (RETINA only as triangle)
      
      
      # plot labels
      xlab(paste0("PC7 (",round(pc7_percent,1),"%", ")" )) +
      ylab(paste0("PC8 (",round(pc8_percent,1),"%", ")" )) +
      labs(title = "PCA (PC7 vs PC8) of rlog-transformed data") +
      
      # custom colors
      scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#56B4E9"))
    
    print(p)
    
    ggsave(filename=paste(myOutputDir,
                          "/PCA_",
                          "PC7vsPC8.pdf",
                          sep=""),
           plot=p,
           width = 5, height = 4, dpi = 300, units = "in", device='pdf')
    
    ##   PCA Plot for PC9 vs PC10
    pc9_percent = summary(myPCA)[[6]][26]*100
    pc10_percent = summary(myPCA)[[6]][29]*100
    p <- ggplot(  pc,
                  aes(x=PC9,y=PC10, color = experimental_conditions, shape = tissue)) +   #notice adding onto object
      geom_point(size = 4) +
      # scale_shape_manual(values=c(17)) +                                 ## (RETINA only as triangle)
      # guides(colour = guide_legend(override.aes = list(shape = 17))) +   ## (RETINA only as triangle)
      
      
      # plot labels
      xlab(paste0("PC9 (",round(pc9_percent,1),"%", ")" )) +
      ylab(paste0("PC10 (",round(pc10_percent,1),"%", ")" )) +
      labs(title = "PCA (PC9 vs PC10) of rlog-transformed data") +
      
      # custom colors
      scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#56B4E9"))
    
    print(p)
    
    ggsave(filename=paste(myOutputDir,
                          "/PCA_",
                          myExperimentalCondition,
                          "_PC9vsPC10.pdf",
                          sep=""),
           plot=p,
           width = 5, height = 4, dpi = 300, units = "in", device='pdf')
    

  }
  
  
  
  
  
  
  
  
  print (paste("== LOOP ITERATION COMPLETE for  ----- ", dds_currentExpCond$condition[1],"_", dds_currentExpCond$age[1]))
  print ("")
}




















if (TRUE) {stop("End of script?")} #It should stop here; the stuff below not for sure working!!
print("Script did NOT end!") # but it doesn't, because this line is printed!
















#########################################v
#  PCA stuff on bigger picture like of retina vs hypo; (12 or 24 samples)
#

#   calculate PCA of rlog of dds
#   (spot to choose if is both tissues or retina-only or hypothalamus-only)
dds_rlog = assay(rlog(dds,blind=FALSE))   # gets the counts of rlog ---> to pca
# dds_rlog = assay(rlog(dds_retina,blind=FALSE))   # gets the counts of rlog ---> to pca
# dds_rlog = assay(rlog(dds_hypothalamus,blind=FALSE))   # gets the counts of rlog ---> to pca


dim(dds_rlog)
head(dds_rlog)[,1:4]         # function from DEseq2
#<NOTICE> the value is equal to raw counts divided by the size factor
                                                  # hypothalamus_R1 hypothalamus_R2 hypothalamus_R3 hypothalamus_R4
                              # ENSMUSG00000000001        8.431875        9.267349        8.244888        9.184564
                              # ENSMUSG00000000003        0.000000        0.000000        0.000000        0.000000
                              # ENSMUSG00000000028        4.698501        4.883920        4.795032        4.720438
                              # ENSMUSG00000000031        5.218583        5.275911        4.345575        4.793019
                              # ENSMUSG00000000037        4.268953        4.968680        4.258555        5.228940
                              # ENSMUSG00000000049        2.347812        2.200285        1.965741        2.180803





# Remove the zero counts 
numericData = dds_rlog
nonzeroData = numericData[rowSums(numericData==0, na.rm=TRUE)<ncol(numericData), ]
dim(numericData)
nonzeroDataTranspose = t(nonzeroData)
myPCA = prcomp(nonzeroDataTranspose, center=T, scale=T)
myPCA_scale_df = as.data.frame(myPCA["scale"])
myPCA_rotation_df = as.data.frame(myPCA["rotation"])
head(myPCA_rotation_df,1) 



### save this data to a table; 
#   (NOTE - this has modification to allow geneid in topleft corner for dt analysis)
write.table(
            data.frame("geneid"=rownames(nonzeroData),format(nonzeroData,digits=4,scientific = FALSE)), 
            file="data_rlog_counts.txt", 
            sep="\t", 
            quote=F, 
            row.names=FALSE)


# What percentage of the total variance is included in the first two principal components? 
# Make a barplot showing the variance in each component (aka "Scree Plot").
myBarplot = barplot(summary(myPCA)[[6]][c(2,5,8,11,14,17,20,23,26,29)],
        ylim=c(0,1.0),
        ylab="Variance due to Principal Component",
        xlab=c("Principle Component (Total Variance)"),
        names.arg = c(gsub(" ", "", paste("PC1 \n (", summary(myPCA)[[6]][3],")")),
                      gsub(" ", "", paste("PC2 \n (", summary(myPCA)[[6]][6],")")),
                      gsub(" ", "", paste("PC3 \n (", summary(myPCA)[[6]][9],")")), 
                      gsub(" ", "", paste("PC4 \n (", summary(myPCA)[[6]][12],")")),
                      gsub(" ", "", paste("PC5 \n (", summary(myPCA)[[6]][15],")")),
                      gsub(" ", "", paste("PC6 \n (", summary(myPCA)[[6]][18],")")),
                      gsub(" ", "", paste("PC7 \n (", summary(myPCA)[[6]][21],")")),
                      gsub(" ", "", paste("PC8 \n (", summary(myPCA)[[6]][24],")")),
                      gsub(" ", "", paste("PC9 \n (",  summary(myPCA)[[6]][27],")")),
                      gsub(" ", "", paste("PC10 \n (", summary(myPCA)[[6]][30],")"))
                      ), 
        col="red",
        main="Barplot of Proportion Variance due to each Principle Component"
)
abline(h=0)

## Add text at top of bars
individualProportions = c(summary(myPCA)[[6]][2],
                          summary(myPCA)[[6]][5],
                          summary(myPCA)[[6]][8],
                          summary(myPCA)[[6]][11],
                          summary(myPCA)[[6]][14],
                          summary(myPCA)[[6]][17],
                          summary(myPCA)[[6]][20],
                          summary(myPCA)[[6]][23],
                          summary(myPCA)[[6]][26],
                          summary(myPCA)[[6]][29]
                          )
text(x = myBarplot, y = individualProportions, label = individualProportions, pos = 3, cex = 0.8, col = "red")


myPCA$x[1:3,1:4]
                # PC1       PC2       PC3        PC4
                # hypothalamus_R1  -99.98774  48.70289  43.90888 -38.931630
                # hypothalamus_R2  -87.90940 -73.68699 -20.69452   3.459668
                # hypothalamus_R3 -102.12399  91.97317  89.15929  20.530111

pc <- as.data.frame(myPCA$x)

# for multifactorial experimental design   (BY-COLOR)
#     can use different sample types based on mutation and age
pc$experimental_conditions <- ifelse(  grepl("retina_R1[012]", row.names(pc)), "wt_3mo",
                                       ifelse(  grepl("retina_R[456]", row.names(pc)), "tb_3mo",
                                                ifelse(  grepl("retina_R[789]", row.names(pc)), "wt_24d",
                                                    ifelse(  grepl("retina_R[123]", row.names(pc)), "tb_24d",
                                                         
                                                         
                                                         ifelse(  grepl("hypothalamus_R1[012]", row.names(pc)), "wt_3mo",
                                                                  ifelse(  grepl("hypothalamus_R[456]", row.names(pc)), "tb_3mo",
                                                                           ifelse(  grepl("hypothalamus_R[789]", row.names(pc)), "wt_24d",
                                                                                    ifelse(  grepl("hypothalamus_R[123]", row.names(pc)), "tb_24d",
                                                                                             "_TODO1_"))))))))
# for multifactorial experimental design   (BY-SHAPE)
#     can use hypothalamus vs retina here
pc$tissue <- ifelse( grepl("hypothalamus", row.names(pc), fixed=TRUE), 
                     "hypothalamus",
                     "retina")




##   PCA Plot for PC1 vs PC2
library(ggplot2)
pc1_percent = summary(myPCA)[[6]][2]*100
pc2_percent = summary(myPCA)[[6]][5]*100
ggplot( pc,
        aes(x=PC1,y=PC2, color = experimental_conditions, shape = tissue)) +   #notice adding onto object
        geom_point(size = 4) +
        # scale_shape_manual(values=c(17)) +                                 ## (RETINA only as triangle)
        # guides(colour = guide_legend(override.aes = list(shape = 17))) +   ## (RETINA only as triangle)
  
      
        # plot labels
        xlab(paste0("PC1 (",round(pc1_percent,1),"%", ")" )) +
        ylab(paste0("PC2 (",round(pc2_percent,1),"%", ")" )) +
        labs(title = "PCA (PC1 vs PC2) of rlog-transformed data") +
  
        # custom colors
        scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#56B4E9"))

##   PCA Plot for PC3 vs PC4
pc3_percent = summary(myPCA)[[6]][8]*100
pc4_percent = summary(myPCA)[[6]][11]*100
ggplot( pc,
        aes(x=PC3,y=PC4, color = experimental_conditions, shape = tissue)) +   #notice adding onto object
        geom_point(size = 4) +
        # scale_shape_manual(values=c(17)) +                                 ## (RETINA only as triangle)
        # guides(colour = guide_legend(override.aes = list(shape = 17))) +   ## (RETINA only as triangle)
    
        # plot labels
        xlab(paste0("PC3 (",round(pc3_percent,1),"%", ")" )) +
        ylab(paste0("PC4 (",round(pc4_percent,1),"%", ")" )) +
        labs(title = "PCA (PC3 vs PC4) of rlog-transformed data on retina samples") +
  
        # custom colors
        scale_color_manual(values = c("#0072B2", "#E69F00", "#D55E00", "#56B4E9"))
        

# (to target candidate genes based on separation in retina samples of tub vs wt)
### <source> https://www.biostars.org/p/289196/
###    <ex code> data.frame(sort(abs(project.pca$rotation[,"PC4"]), decreasing=TRUE)[1:50])

top500gene = data.frame(sort(abs(myPCA$rotation[,"PC4"]), 
                            decreasing=TRUE)[1:500])
colnames(top500gene) = "pca_rotation"
updated_top500gene = top500gene
updated_top500gene["ensembl_gene_id"] = row.names(top500gene)
rownames(updated_top500gene) = 1:nrow(top500gene)
colnames(updated_top500gene) 
  

### get table of all gene names and their gene_id from ensemble
### <source> https://www.biostars.org/p/337072/#337073
library( "biomaRt" )
mart = useMart('ensembl')
listDatasets(mart)[107,]  
ensembl = useMart( "ensembl", dataset = "mmusculus_gene_ensembl" )  
listAttributes(ensembl)     # choose attributes of your interest
subsetOfAttributes <- getBM( attributes = c("ensembl_gene_id",
                              "external_gene_name",
                              "description",
                              "refseq_mrna"
                              
                              
                              
                              ),
               values = row.names(top500gene),
               mart = ensembl) 

# match the two previous files
require(data.table)
head(updated_top500gene,4)
head(subsetOfAttributes,4)
a <- data.table(updated_top500gene)
setkey(a,ensembl_gene_id)
b <- data.table(subsetOfAttributes)
setkey(b,ensembl_gene_id)
matched_df = b[a]
table_to_write = matched_df[order(matched_df$pca_rotation, decreasing = TRUE),]
write.table(table_to_write, file="data_top500genes_PC4retina", sep="\t", quote=F, col.names=NA)







## to plot in 3D; biometry way that can self rotate::
# (requires for mac) https://www.xquartz.org/
# library(rgl)

# #(tutorial code)
# iris.pca=prcomp(iris[,1:4], center=T, scale=T)
# color=c(rep("black", 50), rep("red", 50), rep("green", 50))
# plot3d(iris.pca$x[,1:3], col=color, size=10)

#(for my code)
# color=c(rep("green", 9), rep("blue", 9))
# plot3d(myPCA$x[,1:3], type="s", col=color, size=4)
# plot3d(myPCA$x[,2:4], type="s", col=color, size=4)
# plot3d(myPCA$x[,3:5], type="s", col=color, size=4)




### more PCA stuff
xx = assay(rlog(dds))   # gets the counts of rlog ---> to pca
yy = log(counts(dds,normalized=T) + 1)     # own log; not rlog
