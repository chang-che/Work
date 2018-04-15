##  Using whole blood exoncounts, filter low count genes 

library(DESeq2)
library(edgeR)
dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_RemovedDuplicates/'
load(file=paste(dataPathRNASeq,'ExonCounts.RData',sep='')) ## whole blood (processed)

group = rep(0, 330)
y = DGEList(counts=exoncounts,group=group)
keep <- rowSums(cpm(y)>1) >= 5
save(keep, file="/export/home/xurren/WTCProject/Data/edgeRfilter_exoncounts_330_WB.RData")


