library(dplyr)
library(stats)

# inputFile = read.csv("phages.csv")
load("18samples_retinaVsHypothalamus.deseq2.dds.RData")
ls()
cat("What is the dimensionality of this dataset? ", 
    dim(dds)[1], " rows X ", dim(dds)[2], " cols" , "\n")
# r


# save some variables for later
# phage.DNApackagingStrategy = as.factor(inputFile$DNA.packaging.strategy)
# phage.Institution = as.factor(inputFile$Institution)
# phage.Family = as.factor(inputFile$Family)
# phage.Source = as.factor(inputFile$Source)

# Remove the non-numeric variables. What is the dimensionality of the dataset now?
numericData = inputFile %>% dplyr::select(where(is.numeric))

###
# Perform a principal components analysis on the numerical variables. 
#
plot(numericData)

# What are the eigenvalues and principal components?
phage.pca = prcomp(numericData, center=T, scale=T)  
  # provides std dev (=sqrt of the eigenvalues)

  # to get eigenvalues, square the std dev
eigenValues = phage.pca$sdev * phage.pca$sdev

# What percentage of the total variance is included in the first two principal components? 
# Make a barplot showing the variance in each component.
barplot(summary(phage.pca)[[6]][c(2,5,8,11)],
        ylim=c(0,1.0),
        ylab="% Variance due to Principal Components",
        xlab=c("Principle Component"),
        names.arg = c("PC1", "PC2", "PC3", "PC4"), 
        col="red",
        main="Barplot of Proportion Variance due to each Principle Component"
       )
abline(h=0)


# Make plots of the first two principal components against each other, 
# where the phages are colored according to DNA packaging strategy, 
# institution, family, and source (so 4 plots total).
# (in other words) check plots based on different possible clustering from inputData non-numeric columns

plot(phage.pca$x[,1:2], col=phage.DNApackagingStrategy, pch=20, # best looks like this one
     main="PCA Plot of Clustering based on DNA Packaging Strategy")   
legend(1,1,
       legend = levels(phage.DNApackagingStrategy),
       col=1:length(levels(phage.DNApackagingStrategy)),
       pch=20
      )
plot(phage.pca$x[,1:2], col=phage.Institution, pch=20, 
     main="PCA Plot of Clustering based on Institution")   
legend(1,1,
       legend = levels(phage.Institution),
       col=1:length(levels(phage.Institution)),
       pch=20
)

plot(phage.pca$x[,1:2], col=phage.Family, pch=20, 
     main="PCA Plot of Clustering based on Family")   
legend(1,1,
       legend = levels(phage.Family),
       col=1:length(levels(phage.Family)),
       pch=20
)

plot(phage.pca$x[,1:2], col=phage.Source, pch=20, 
     main="PCA Plot of Clustering based on Source")   
legend(1,1,
       legend = levels(phage.Source),
       col=1:length(levels(phage.Source)),
       pch=20
)


# What patterns do you see in the plots, if any? 
# What can you say about these phages? 
# How do they cluster based on the PCA?



