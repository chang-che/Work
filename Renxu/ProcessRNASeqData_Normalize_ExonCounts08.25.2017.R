##  normalize the ExonCounts533.RData  

library(DESeq2)

dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/'
load(file=paste(dataPathRNASeq,'ExonCounts533.RData',sep=''))
##  check missing values
sum(is.na(exoncounts))    # there are no missing values

##  make a DESeq data set
coldata = data.frame(batch = c(rep("old", 330), rep("new", 203)), 
					 row.names = colnames(exoncounts))
dds = DESeqDataSetFromMatrix(countData = exoncounts,
                              colData = coldata,
                              design = ~ batch)
##  estimate size factors and calculate normalized counts  
dds = estimateSizeFactors(dds)                             
exoncounts_normalized = counts(dds, normalized = TRUE)
# dim(genecounts_normalized)
dataPath <- '/export/home/xurren/WTCProject/Data/'
save(exoncounts_normalized, file=paste(dataPath, "ExonCountsNormalize533.RData", sep = ''))





