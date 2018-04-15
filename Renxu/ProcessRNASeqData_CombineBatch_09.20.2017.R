### old 330 samples
dataPathRNASeq_old <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_RemovedDuplicates/'
load(file=paste(dataPathRNASeq_old,'clinical_rnaseq_16Aug2017.RData',sep=''))
load(paste(dataPathRNASeq_old,'GeneCounts.RData',sep='')) ### load genecounts
load(paste(dataPathRNASeq_old,'FPKM.RData',sep='')) ## load fpkm


### new 203 samples
dat1 <- read.csv('/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/RQ-009443_Quantitation_data_hg19/GeneBodyCounts_Wholeblood_2017.csv',header=T) #****#
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

sum(rownames(genecounts)!=dat1[match_counts,1])
sum(rownames(fpkm)!=as.character(fpkm2[match_fpkm,1]))
match_2 = match(rownames(genecounts), rownames(fpkm))
fpkm <- fpkm[match_2,]
sum(rownames(fpkm)!=rownames(genecounts))

dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/'
save(genecounts,file=paste(dataPathRNASeq,'GeneBodyCounts533.RData',sep=''))
save(fpkm,file=paste(dataPathRNASeq,'FPKM533.RData',sep=''))


#> dim(fpkm)
#[1] 25830   533
#> genecounts <- cbind(genecounts,dat2)
#> dim(genecounts)
#[1] 25830   533

### check if these samples are outliers: W15659, W29979, W26904

#id <- which(rownames(genecounts)=="FKBP5")
#genecounts[id,]

### note, use xcell to estimate the cell types on FPKM, store the estimated cell t
#ype matrix into dataPathRNASeq <- /export/home/pfkuan/WTCproject/Epigenetics/Data
#/RNASeq/ProcessedData_533Samples/

### process clinical data ###
dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/'
load(file=paste(dataPathRNASeq,'GeneCounts533.RData',sep=''))
load(file=paste(dataPathRNASeq,'FPKM533.RData',sep=''))

clin.unsort <- read.delim('/export/home/pfkuan/WTCproject/Epigenetics/Data/Clinical/clinical533_unsorted.txt',header=T,sep='\t')

sum(colnames(genecounts)!=as.character(clin.unsort$w_mrn))

match_clin <- match(colnames(genecounts), as.character(clin.unsort$w_mrn))
clin <- clin.unsort[match_clin,]
sum(colnames(genecounts)!=as.character(clin$w_mrn))


save(clin,file=paste(dataPathRNASeq,'clinical533_083117.RData',sep=''))

### cross-check with older data ###
dataPathRNASeq_old <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_RemovedDuplicates/'
load(file=paste(dataPathRNASeq_old,'clinical_rnaseq_16Aug2017.RData',sep=''))

clin.tmp <- clin[match(as.character(clinical.rnaseq$w_mrn),as.character(clin$w_mrn)),]

xtabs(~clin.tmp$PTSD_SCID_3grp+clinical.rnaseq[,57])

#> xtabs(~clin.tmp$PTSD_SCID_3grp+clinical.rnaseq[,57])
#                       clinical.rnaseq[, 57]
#clin.tmp$PTSD_SCID_3grp   1   2   3
#                      1  72   2   0
#                      2  11  40   1
#                      3   0   0 202

xtabs(~clin.tmp$PTSD_SCID_3grp[clin$batch==3]+clinical.rnaseq[clin$batch==3,57])
