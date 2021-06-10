
#(usage) Rscript convertCountsToRPKM.r arg1 arg2 arg3

args<-commandArgs(TRUE)


FILE_INPUT_COUNTS = args[1] 	# (any counts file from HTSeq) 7WT24-R1_S5_counts.txt ; make sure it has removed the lines at the end starting with 
FILE_INPUT_GENE_LENGTH = args[2] 	# ./Mus_musculus.GRCm38.geneLength.txt
FILE_OUTPUT_RPKM = args[3]


file_count<- read.delim(FILE_INPUT_COUNTS,header=F,sep="\t")
file_len<- read.delim(FILE_INPUT_GENE_LENGTH,header=F,sep="\t")
colnames(file_len)<- c("GeneName","Len")
colnames(file_count)<- c("GeneName","Count")  
file_count<-file_count[ !grepl("__", file_count$GeneName) ,]
total_count<- sum(file_count$Count)
oneB<-10^9
finallist <- merge(file_len,file_count,by="GeneName")
finallist$RPKM<-0
finallist[,2:4] <- (sapply(finallist[,2:4], as.double))
finallist$RPKM<- (oneB*finallist$Count)/(total_count*finallist$Len)
#finallist<-finallist[finallist$RPKM>1,]

write.table(finallist,file=FILE_OUTPUT_RPKM,sep="\t", col.names = T, row.names = F, quote = F)