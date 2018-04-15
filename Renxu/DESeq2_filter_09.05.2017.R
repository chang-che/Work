

library("DESeq2")
library("edgeR")
load(file = '/export/home/xurren/WTCProject/Data/clinical533_reordered.RData')  ## load clinical data
dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/'
load(file=paste(dataPathRNASeq,'ExonCounts533.RData',sep=''))
load(file = '/export/home/xurren/WTCProject/Data/TrainID_and_TestID.RData')  ## load train and test ID. testID2 is 135 test data, testID1 is 203 samples from new batch

genecounts = exoncounts

batch = c(rep("old", 330), rep("new", 203))
y = DGEList(counts = exoncounts, group = batch)

threshold = seq(from = 0, to = 50, by = 5)
keep_num = rep(0, length(threshold))
for(i in 1:length(threshold)){
  keep = rowSums(cpm(y)>1) >= threshold[i]
  keep_num[i] = sum(keep)
  
}
pdf(file = "/export/home/xurren/WTCProject/Results/numofgenes_vs_threshold.pdf")
plot(threshold, keep_num, xlab = "threshold(number of samples cpm greater than 1)", ylab = "number of genes", main = "genes filter")
dev.off()

