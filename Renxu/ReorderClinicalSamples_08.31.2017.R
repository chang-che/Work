

dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/'
load(file=paste(dataPathRNASeq,'FPKM533.RData',sep=''))
load(file=paste(dataPathRNASeq,'GeneCounts533.RData',sep=''))

clinical = read.csv("/export/home/xurren/WTCProject/Data/clinical533_original.csv")
dim(clinical)

match_samples = match(colnames(fpkm), clinical[, 1])
clinical = clinical[match_samples, ]
sum(clinical[,1]!=colnames(fpkm))

save(clinical, file = "/export/home/xurren/WTCProject/Data/clinical533_reordered.RData")

