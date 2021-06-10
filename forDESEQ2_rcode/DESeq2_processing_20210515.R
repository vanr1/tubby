# _TODO_ check if this below output is the same as the file i'm looking at??

library("DESeq2")
library("ggplot2")
## values coming out of deseq2.dds
# load("18samples_retinaVsHypothalamus.deseq2.dds.RData")
load("24samples_retinaVsHypothalamus.deseq2.dds.RData")
ls()

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




# (here) get gene symbol and description; add it to normalized_counts before outputted
all_genes <- sapply( strsplit( rownames(normalized_counts), split="\\+" ), "[", 1)
library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "mmusculus_gene_ensembl" )
# listAttributes(ensembl)     # choose attributes of your interest
genemap <- getBM(  attributes = c("ensembl_gene_id","external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = all_genes,
                   mart = ensembl)
idx <- match( all_genes, genemap$ensembl_gene_id )
gene_name <- genemap$external_gene_name[ idx ]
head(all_genes, 4)
head(gene_name, 4)

# add genename to first column, genedescription to last
normalized_counts = cbind(gene_name, normalized_counts) 

head(normalized_counts)[,1:4]         
write.table(
  data.frame("geneid"=rownames(normalized_counts),normalized_counts), 
  file="data_normalized_counts.txt", 
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



######
#  modify dds by adding other factors
#
dds$age = ifelse(  grepl(".*_R1[012]", colnames(dds)), "3mo",
                   ifelse(  grepl(".*_R[456]", colnames(dds)), "3mo",
                            ifelse(  grepl(".*_R[789]", colnames(dds)), "24d",
                                     ifelse(  grepl(".*_R[123]", colnames(dds)), "24d",
                                              "_TODO1_"))))
dds$age = as.factor(dds$age)
dds$genotype = ifelse(  grepl(".*_R1[012]", colnames(dds)), "wt",
                        ifelse(  grepl(".*_R[456]", colnames(dds)), "tb",
                                 ifelse(  grepl(".*_R[789]", colnames(dds)), "wt",
                                          ifelse(  grepl(".*_R[123]", colnames(dds)), "tb",
                                                   "_TODO1_"))))
dds$genotype = as.factor(dds$genotype)

colData(dds)[ ,c("genotype","age")]  # can also get tissue which is called "condition"



######
#  subsetting dds to different variables ; drop other variables to be safe right after
#

## (by tissue) two groups of 12 samples each
dds_retina = dds[,dds$condition == "retina"]
dim(dds_retina)
dds_hypothalamus = dds[,dds$condition == "hypothalamus"]
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

# (function) relevel  
# convenient way to make sure that "wt" is the first level in the treatment factor,
# so the default log2FC calculated as treatment over control and not other way around
dds_retina_24d$genotype = relevel(dds_retina_24d$genotype, "wt")
dds_retina_3mo$genotype = relevel(dds_retina_3mo$genotype, "wt")
dds_hypothalamus_24d$genotype = relevel(dds_hypothalamus_24d$genotype, "wt")
dds_hypothalamus_3mo$genotype = relevel(dds_hypothalamus_3mo$genotype, "wt")

as.data.frame( colData(dds_retina_24d))
as.data.frame( colData(dds_retina_3mo))
as.data.frame( colData(dds_hypothalamus_24d))
as.data.frame( colData(dds_hypothalamus_3mo))



#####
#  Running the DESEQ2 pipeline 
####

# #    change design to look at it based on genotype; rerun DESeq
design(dds_retina_24d) <- formula(~genotype)
dds_retina_24d <- DESeq(dds_retina_24d)
design(dds_retina_3mo) <- formula(~genotype)
dds_retina_3mo <- DESeq(dds_retina_3mo)
design(dds_hypothalamus_24d) <- formula(~genotype)
dds_hypothalamus_24d <- DESeq(dds_hypothalamus_24d)
design(dds_hypothalamus_3mo) <- formula(~genotype)
dds_hypothalamus_3mo <- DESeq(dds_hypothalamus_3mo)


listOfExpCond = list(dds_retina_24d, dds_retina_3mo, dds_hypothalamus_24d, dds_hypothalamus_3mo)
for(dds_currentExpCond in listOfExpCond)
{
  # print ("[DEBUG] performing DESeq DEG analysis pipeline to output two files for")
  # print (dds_currentExpCond$condition[1])
  # print (dds_currentExpCond$age[1])
  
  res = results(dds_currentExpCond)
  # res = results( dds_retina_24d )   ## CHOICE ONE; of the next 4 lines; how it was before turn into a loop
  # res = results( dds_retina_3mo )  
  # res = results( dds_hypothalamus_24d )  
  # res = results( dds_hypothalamus_3mo )  
  
  # debug stuff
  sum (res$pvalue < 0.01, na.rm=TRUE)
  sum (res$pvalue < 0.05, na.rm=TRUE)
  length (res [ which(res$padj <0.10 ), ] )
  length (res [ which(res$padj <0.05 ), ] )
  res
  table (is.na (res$pvalue))
  
  #    to add gene names
  res$ensembl <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1)
  
  library( "biomaRt" )
  ensembl = useMart( "ensembl", dataset = "mmusculus_gene_ensembl" )
  # listAttributes(ensembl)     # choose attributes of your interest
  genemap <- getBM(  attributes = c("ensembl_gene_id","external_gene_name","interpro_description"),
                     filters = "ensembl_gene_id",
                     values = res$ensembl,
                     mart = ensembl)
  idx <- match( res$ensembl, genemap$ensembl_gene_id )
  res$gene_name <- genemap$external_gene_name[ idx ]
  res$gene_description <- genemap$interpro_description[ idx ]
  head(res, 4)
  
  # subset to get top down/up-regulated genes
  resSig = res [ which(res$padj <0.1 ), ]             # fdr correction
  dim(resSig)
  head (resSig[ order( resSig$log2FoldChange ), ] )   # for downregulated genes  (double-check not flipped; how?)
  head (rev(resSig[ order( resSig$log2FoldChange ), ] ))   # for upregulated genes (double-check not flipped; how?)
  
  
  
  ### save this data to two tables; 
  highestLog2FC_genes = rev(resSig[ order( resSig$log2FoldChange ), ] )
  
  allLog2FC_genes = rev( res[ order( res$log2FoldChange ), ] )
  row.has.na <- apply(allLog2FC_genes, 1, function(x){any(is.na(x))})
  print (paste("    droppingNAs --- $sum(row.has.na) ", sum(row.has.na) ))
  allLog2FC_genes_filtered = allLog2FC_genes[!row.has.na,] 
  
  lowestLog2FC_genes = resSig[ order( resSig$log2FoldChange ), ]
  
  
  print ("[DEBUG] about to write two tables of DEG ordered upregulated then downregulated ")
  print (paste("        (expr conds) ", dds_currentExpCond$condition[1], dds_currentExpCond$age[1] ))
  print (paste("        (# of DEGs) ", dim(highestLog2FC_genes)[1] ))
  write.table(highestLog2FC_genes[ , c("gene_name", "ensembl", "log2FoldChange", "pvalue", "padj", "gene_description")],
              file=paste("output/DEGs",
                         dds_currentExpCond$condition[1],
                         dds_currentExpCond$age[1],
                         "sortedByHighestLog2FC.txt",
                         sep="_"),
              sep="\t",
              quote=F,
              row.names=FALSE)
  
  # (extra info) full list of genes and also without NA
  write.table(allLog2FC_genes[ , c("gene_name", "ensembl", "log2FoldChange", "pvalue", "padj", "gene_description")],
              file=paste("output/DEGs",
                         dds_currentExpCond$condition[1],
                         dds_currentExpCond$age[1],
                         "allGenes.txt",
                         sep="_"),
              sep="\t",
              quote=F,
              row.names=FALSE)
  write.table(allLog2FC_genes_filtered[ , c("gene_name", "ensembl", "log2FoldChange", "pvalue", "padj", "gene_description")],
              file=paste("output/DEGs",
                         dds_currentExpCond$condition[1],
                         dds_currentExpCond$age[1],
                         "allGenesFiltered.txt",
                         sep="_"),
              sep="\t",
              quote=F,
              row.names=FALSE)
  
  write.table(lowestLog2FC_genes[ , c("gene_name", "ensembl", "log2FoldChange", "pvalue", "padj", "gene_description")],
              file=paste("output/DEGs",
                         dds_currentExpCond$condition[1],
                         dds_currentExpCond$age[1],
                         "sortedByLowestLog2FC.txt",
                         sep="_"),
              sep="\t",
              quote=F,
              row.names=FALSE)
  
  
  # some quick PCA here only looking at 6 samples
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
    file=paste("RLOGDATA",
               dds_currentExpCond$condition[1],
               dds_currentExpCond$age[1],
               "sorted.txt",
               sep="_"),
    sep="\t", 
    quote=F, 
    row.names=FALSE)
  
  
  # What percentage of the total variance is included in the first two principal components? 
  # Make a barplot showing the variance in each component (aka "Scree Plot").
  print ("        + generating bar plot...")
  pdf(paste("output/SCREEPLOT",
            dds_currentExpCond$condition[1],
            dds_currentExpCond$age[1],
            "processed.pdf",
            sep="_"))
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
  
  
  ## now to create the PCA plots; 
  
  print ("        + generating PCA plot...")
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
  
  ggsave(filename=paste("output/PCAPLOT",
                        dds_currentExpCond$condition[1],
                        dds_currentExpCond$age[1],
                        "PC1vsPC2.pdf",
                        sep="_"),
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
  
  ggsave(filename=paste("output/PCAPLOT",
                        dds_currentExpCond$condition[1],
                        dds_currentExpCond$age[1],
                        "PC3vsPC4.pdf",
                        sep="_"),
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
  
  ggsave(filename=paste("output/PCAPLOT",
                        dds_currentExpCond$condition[1],
                        dds_currentExpCond$age[1],
                        "PC3vsPC5.pdf",
                        sep="_"),
         plot=p,
         width = 5, height = 4, dpi = 300, units = "in", device='pdf')
  
  
  
  # some additional diagonostic plots using DESeq2 library to quickly make
  print ("        + generating MA plot...")
  pdf(paste("output/MAPLOT",
            dds_currentExpCond$condition[1],
            dds_currentExpCond$age[1],
            "diagnosticPlot.pdf",
            sep="_"))
  plotMA (res, 
          ylim =c(-10,10 ),
          
          main=paste("MA plot of",
                     dds_currentExpCond$condition[1],
                     "at",
                     dds_currentExpCond$age[1],
                     "samples",
                     sep=" ")
  )
  
  # (defn of x-axis) It indicates the average expression between the groups you 
  #   compared, which is simply the mean of the normalized counts of the samples 
  #   used in the comparison group. The "meaning" is that genes with higher counts 
  #   are more trustworthy (=have higher statistical power) than low count genes. 
  #   Biologically this is because either the gene is lowly expressed and/or short 
  #   and therefore intrinsically gets lower counts than a longer one at equal 
  #   expression level.
  
  dev.off()
  
}







if (TRUE) {stop("End of script?")} #It should stop here
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
                                            "description"
                                            
                                            
                                            
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
