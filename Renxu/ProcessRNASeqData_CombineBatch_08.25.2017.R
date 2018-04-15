### old 330 samples
dataPathRNASeq_old <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_RemovedDuplicates/'
load(file=paste(dataPathRNASeq_old,'clinical_rnaseq_16Aug2017.RData',sep=''))
load(paste(dataPathRNASeq_old,'GeneCounts.RData',sep='')) ### load genecounts
load(paste(dataPathRNASeq_old,'FPKM.RData',sep='')) ## load fpkm


### new 203 samples
dat1 <- read.csv('/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/RQ-009443_Quantitation_data_hg19/Counts.csv',header=T)
dat2 <- dat1[,-1]

fpkm2 <- read.csv('/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/RQ-009443_Quantitation_data_hg19/RPKMedgeR.csv',header=T)
fpkm3 <- fpkm2[,-1]

#clinical1 <- cbind(c(1:533),data.frame(c(as.character(clinical.rnaseq$w_mrn),colnames(dat2))),rep(c(1,2),c(dim(clinical.rnaseq)[1],dim(dat2)[2])))
#colnames(clinical1) <- c('No','w_mrn','batch')

#write.csv(clinical1,'/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/Merged533Samples/clinicalVar533Samples.csv',row.names=F)

### checking if the rows of 330 samples matched with rows of 203 samples
# tmp <- cbind(rownames(genecounts),as.character(dat1[,1]))
# tmp1 <- cbind(rownames(fpkm),as.character(fpkm2[,1]))
sum(rownames(genecounts) %in% as.character(dat1[,1]))  
sum(rownames(genecounts)!=as.character(dat1[,1]))   ##  this value is not 0, which means the two gene lists are not in the same order
sum(rownames(fpkm) %in% as.character(fpkm2[,1]))
sum(rownames(fpkm)!=as.character(fpkm2[,1]))   ##  this value is also not 0


### reorder the genes and combine the samples
match_counts = match(rownames(genecounts), as.character(dat1[,1]))
genecounts = cbind(genecounts, dat2[match_counts, ])
# dim(genecounts)
match_fpkm = match(rownames(fpkm), as.character(fpkm2[,1]))
fpkm = cbind(fpkm, fpkm3[match_fpkm, ])
# dim(fpkm)

match_2 = match(rownames(genecounts), rownames(fpkm))
fpkm <- fpkm[match_2,]
sum(rownames(fpkm)!=rownames(genecounts))

dataPathRNASeq <- '/export/home/xurren/WTCProject/Data/'
save(genecounts,file=paste(dataPathRNASeq,'GeneCounts533.RData',sep=''))
save(fpkm,file=paste(dataPathRNASeq,'FPKM533.RData',sep=''))


#> dim(fpkm)
#[1] 25830   533
#> genecounts <- cbind(genecounts,dat2)
#> dim(genecounts)
#[1] 25830   533

### check if these samples are outliers: W15659, W29979, W26904

#id <- which(rownames(genecounts)=='FKBP5')
#genecounts[id,]

### note, use xcell to estimate the cell types on FPKM, store the estimated cell type matrix into dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/
