

cat("\n")
cat("     =======================================================\n")
cat("     [Starting R Script] -- INSIDE RSCRIPT myScript_deseq2_processing.r \n")
cat("\n")

# BiocManager::install(c("vsn"))
# BiocManager::install(c("pheatmap"))
# BiocManager::install(c("tximport"))
# BiocManager::install(c("tximportData"))
# BiocManager::install(c("GenomicFeatures"))
# BiocManager::install(c("DESeq2"))
# BiocManager::install("org.Mm.eg.db") # for convert ensembleID to gene name and symbols
# BiocManager::install(c("AnnotationDbi"))
# BiocManager::install("gskb")  # for gene ontology

# library("tximport")
# library("readr")
# library("tximportData")

# dir <- system.file("extdata", package="tximportData")

# samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
# samples$condition <- factor(rep(c("A","B"),each=3))
# rownames(samples) <- samples$run
# samples[,c("pop","center","run","condition")]

# dir = system.file("extdata", package = "tximportData")
# list.files(dir)
# library(tximportData)
# dir = "."



# making a variable for our main pipeline path
#TUBBY_PIPELINE = "/data2/han_lab/richardvan/tubby_pipeline"                            #server
#directory = paste(TUBBY_PIPELINE, "/data/deseq2_input_sets/wt_vs_tub_all", sep="")      #server

directory = "./"     #local


# samples is just a text file of how we perform DEG analysis against
# samples <- read.table(file.path(directory, "wt_vs_tub_all.txt"), header = F)
samples_wildtype  = grep("wt",list.files(directory),value=TRUE)
samples_mutant    = grep("tub",list.files(directory),value=TRUE)
cat("Reading in samples ... \n")
cat("Displaying the wildtype 'tub' samples:\n")
cat("  ", samples_wildtype, "\n")
cat("Displaying the mutant 'tub' samples:\n")
cat("  ", samples_mutant, "\n")
cat("Displaying the wildtype 'untreated' samples:\n")
cat("  ", samples_wildtype, "\n")
cat("Displaying the mutant 'treated' samples:\n")
cat("  ", samples_mutant, "\n")
cat("\n")

cat("Building sample tables ... \n")
sampleFiles <- grep("treated",list.files(directory),value=TRUE)
sampleCondition <- sub(".*_(.*treated).*","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)


cat("Loading the library DESEQ2 ...\n")
library("DESeq2")
cat("Finished loading the library DESEQ2 ...\n")
cat("\n")

cat("Building the DESeqDataSet ...")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
cat("Finished building the DESeqDataSe ...\n")
cat("Displaying 'ddsHTSeq' ...\n")
ddsHTSeq
cat("\n")

cat("Running the command for differential expression analysis (alpha = 0.05) ...\n")
dds <- DESeq(ddsHTSeq)
res <- results(dds)
cat("\n")
cat("Displaying 'res' ...\n")
res
cat("\n")
cat("Displaying 'summary(res)' ...\n")
summary(res)
cat("\n")


# cat("Plotting PCA using DESEQ2 and saving to file ...\n")
# plotPCA(vsd, intgroup=c("condition"))
# pdf(file="pca.pdf")

#< HARD QUIT FOR DEBUGING ------------------------------------ >
cat("\n")
cat("[DEBUG] about to hard quit; update declaring variables --\n\n")
quit(status=1)

#<  -------- >



 # transcripts_quant_1Tub24-R1_S4
 # transcripts_quant_2Tub24-R2_S6
 # transcripts_quant_3Tub24-R3_S3
 # transcripts_quant_7WT24-R1_S5
 # transcripts_quant_8WT24-R2_S2
 # transcripts_quant_9WT24-R3_S1

# files <- file.path(dir, "salmon", samples$V1, "quant.sf")
# names(files) <- paste0("sample", 1:6)
# all(file.exists(files))

# library(GenomicFeatures)
# library(readr)
# txdb = makeTxDbFromGFF("gencode.vM25.annotation.gtf.gz", format=c("auto"))
# k <- keys(txdb, keytype = "TXNAME")
# tx2gene <- select(txdb, k, "GENEID", "TXNAME")
# head(tx2gene)

# library(tximport)
# txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
# names(txi)
# head(txi$counts)



library("DESeq2")
samples$condition <- factor(rep(c("A","B"),each=3))

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
dds <- DESeq(ddsTxi)

res05 <- results(dds, alpha=0.05)
summary(res05)
res10 <- results(dds, alpha=0.10)
summary(res10)


# figure out which genes are significant -- start with pvalues for all the genes
library(dplyr)
res05_df = data.frame(res05$baseMean, res05$log2FoldChange, res05$padj)
row.names(res05_df) = row.names(res05)
sorted_res05_df = arrange(res05_df, res05_df[[3]])

# select only the "positives" which in DESEQ2 is defined as those with adj_pvalue < alpha; should do q-value vs p-value??
positives = filter(sorted_res05_df, sorted_res05_df[[3]] < 0.05)
positives10 = filter(sorted_res05_df, sorted_res05_df[[3]] < 0.10)


#  - output ensembleIDs so can manually convert to ensembleIDs
#       http://uswest.ensembl.org/biomart/martview/d1d566b708991b9bc8179645ed4d8e66

#  - convert ensembleIDs to actual genenames
#    (adapted from) http://combine-australia.github.io/RNAseq-R/09-applying-rnaseq-solutions.html
library(org.Mm.eg.db)
library(AnnotationDbi)
library(data.table)
library(gskb)

gene_ids_version = row.names(positives)
# strip off the ensemble gene id version provided by DESEQ2
geneEnsembleIDs =nth(tstrsplit(gene_ids_version, split ="\\."),n=1)
geneEnsembleIDs

#  OUTPUT -  only the ensembleIDs in a table for manual input
write.table(geneEnsembleIDs, 
              file="tubby_significant_genes_byEnsembleIDs.txt", 
              row.names=F, 
              col.names=T,
              sep="\t",
              quote = FALSE
        )
# (note) after write table for this list, need to manually delete 'first row x'



cat("  There were ", length(geneNames), " candidate genes by end of DESeq2", "\n")


columns(org.Mm.eg.db)
AnnotationDbi::keytypes(org.Mm.eg.db)
myKeys = AnnotationDbi::keys(org.Mm.eg.db, "ENSEMBL")
ann <- AnnotationDbi::select(org.Mm.eg.db,
                             keys=geneEnsembleIDs,
                             columns=c("ENSEMBL","ENTREZID","SYMBOL","GENENAME"),
                             keytype="ENSEMBL"
)
head(ann)

cat("  There were ", nrow(na.omit(ann)), " candidate genes by found in org.Mm.eg.db annotation", "\n")
cat("       (notice) some ensembleIDs didn't map so had (", length(geneEnsembleIDs), "-", nrow(na.omit(ann)), ") NA values that are dropped", "\n")
cat("       So only ", nrow(na.omit(ann)), " candidate genes proceed to GO analysis", "\n")

significant_genes = data.frame(ann$SYMBOL, ann$GENENAME, positives)
row.names(significant_genes) = ann$ENSEMBL
colnames(significant_genes) = c("SYMBOL","GENENAME", "BASE_MEAN", "LOG2FC", "P_ADJ")


# OUTPUT -- entire table (these are with bonferroni correction; strict at p-value less than alpha*)
write.table(significant_genes, 
            file="tubby_significant_genes.txt", 
            row.names=T, 
            col.names=NA,
            sep="\t",
            quote = FALSE
           )





# remove the NAs and output the found symbols for Gene Ontology
write.table(na.omit(significant_genes$SYMBOL), 
            file="tubby_significant_genes_symbolsONLY.txt", 
            row.names=FALSE, 
            sep="\t",
            quote = FALSE
)


# this table is outputted so have a list of gene symbols can copy and paste into geneontology.org (for by-hand)
write.table(na.omit(significant_genes$SYMBOL), 
            file="tubby_significant_genes_symbolsONLY.txt", 
            row.names=FALSE, 
            sep="\t",
            quote = FALSE
)


# Gene Ontology programatically
data(mm_GO)
mm_GO[[1]][1:10]
GSEA.prog.loc<- "http://ge-lab.org/gskb/GSEA.1.0.R"

GSEA(
        # Input/Output Files :------------------------------------------------

        # Input gene expression Affy dataset file in RES or GCT format
        #   _TODO____changeThisLater
        input.ds = "./mouse_data.gct", 

        # Input class vector (phenotype) file in CLS format
        #   _TODO____changeThisLater
        input.cls = "./mouse.cls",

        # Gene set database in GMT format
        gs.db = mm_GO,

        # Directory where to store output and results (default: "")
        output.directory = getwd(),
        
        # Program parameters :-----------------------------------------------
        doc.string = "mouse",
        non.interactive.run = T,
        reshuffling.type = "sample.labels",
        nperm = 1000,
        weighted.score.type = 1,
        nom.p.val.threshold = -1,
        fwer.p.val.threshold = -1,
        fdr.q.val.threshold = 0.25,
        topgs = 10,
        adjust.FDR.q.val = F,
        gs.size.threshold.min = 15,
        gs.size.threshold.max = 500,
        reverse.sign = F,
        preproc.type = 0,
        random.seed = 3338,
        perm.type = 0,
        fraction = 1.0,
        replace = F,
        save.intermediate.results = F,
        OLD.GSEA = F,
        use.fast.enrichment.routine = T
        )

# visualizations!!



res <- results(dds)
res
plotMA(res, ylim=c(-5,5))

# extract transformed values
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)


#PCA
plotPCA(vsd, intgroup=c("condition", "V1"))

# effects of transformations on variance
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

#heatmap of the count matrix
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","V1")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)

#heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


cat("\n")
cat("     [Ending R Script]  \n")
cat("     =======================================================\n")
cat("\n")
