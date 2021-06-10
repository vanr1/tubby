
library(tidyverse)

temp_table = read.table('results_18samples_retinaVsHypothalamus__deseq2.pca.vals.txt', header=TRUE)
  ## need to remember to manually delete the top left value "sample" so it reads in correctly


percentVar = c(70,9)
## might need to write function that extracts the pc values;


data_ORIGINAL = load('corinne_pca_data.rds')
# y = load('deseq2.dds.RData')
pcs = prcomp(t(log(hervh_mat+1)), center = TRUE)

percentVar_ORIGINAL = round(((pcs$sdev) ^ 2 / sum((pcs$sdev) ^ 2)* 100), 2) 

pc_ORIGINAL <- as.data.frame(pcs$x)
pc = as.data.frame(temp_table)


#(BY-COLOR)
# (analogous to) pc$cell_type; 
# can use different sample types based on mutation and age
# (consider change to) mutation and age
# predefined for color blind friendly
pc$experimental_conditions <- "_TODO_"
pc$experimental_conditions <- ifelse(  grepl("retina_R[123]", row.names(pc)), "tub_24d",
                              ifelse(  grepl("retina_R[456]", row.names(pc)), "wt_24d",
                              ifelse(  grepl("retina_R[789]", row.names(pc)), "wt_3m",
                          
  
                              ifelse(  grepl("hypothalamus_R[12]", row.names(pc)), "wt_24d",
                              ifelse(  grepl("hypothalamus_R[345]", row.names(pc)), "wt_3m",
                              ifelse(  grepl("hypothalamus_R[678]", row.names(pc)), "tub_24d",
                              ifelse(  grepl("hypothalamus_R[9]", row.names(pc)), "tub_3m",
                                       "_TODO1_")))))))




#(BY-SHAPE)
# (analogous to) pc$study 
# can use hypothalamus vs retina here
# (consider change to) tissue type
pc$tissue <- ifelse( grepl("hypothalamus", row.names(pc), fixed=TRUE), 
                    "hypothalamus",
                    "retina")

png("pca_genes_small.png", height = 1200, width = 1450, res = 225)

## PCA Plot
ggplot( pc,
        aes(x=PC1,y=PC2, color = experimental_conditions, shape = tissue)) +   #notice adding onto object
  geom_point(size = 2) +
  
  
  # plot labels
  xlab(paste0("PC1 (",percentVar[1],"% variance", ")" )) +
  ylab(paste0("PC2 (",percentVar[2],"% variance", ")" )) +
  labs(title = "PCA of rlog-transformed data")

dev.off()