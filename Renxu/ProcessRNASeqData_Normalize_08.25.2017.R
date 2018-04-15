##  normalize the WTCProjectGeneCounts533.RData  

library(DESeq2)

dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/'
load(file=paste(dataPathRNASeq,'GeneCounts533.RData',sep=''))
##  check missing values
sum(is.na(genecounts))    # there are no missing values

##  make a DESeq data set
coldata = data.frame(batch = c(rep("old", 330), rep("new", 203)), 
					 row.names = colnames(genecounts))
dds = DESeqDataSetFromMatrix(countData = genecounts,
                              colData = coldata,
                              design = ~ batch)
##  estimate size factors and calculate normalized counts  
dds = estimateSizeFactors(dds)                             
genecounts_normalized = counts(dds, normalized = TRUE)
# dim(genecounts_normalized)
dataPath <- '/export/home/xurren/WTCProject/Data/'
save(genecounts_normalized, file=paste(dataPath, "GeneCountsNormalize533.RData", sep = ''))





