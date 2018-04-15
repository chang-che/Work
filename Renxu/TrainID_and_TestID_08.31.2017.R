## specify the training ID and testing ID

dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_RemovedDuplicates/'
load(file=paste(dataPathRNASeq,'clinical_rnaseq_30Sept2016_Train.RData',sep='')) ## load clinical.rnaseq 
load(file=paste(dataPathRNASeq,'GeneCounts_Train.RData',sep=''))  ## load genecounts
dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/'
load(file=paste(dataPathRNASeq,'FPKM533.RData',sep=''))  ##  load all 533 samples

clinical.rnaseq[clinical.rnaseq<0] <- NA
clinical.rnaseq[clinical.rnaseq==5555] <- NA

f <- 208
pheno <- clinical.rnaseq[,f] ###this is the old PTSD, where 0:control/never PTSD, 1:current PTSD
phenoName <- colnames(clinical.rnaseq)[f]

idinc <- which(is.na(pheno)==FALSE&clinical.rnaseq$sex==0&clinical.rnaseq$batch==3)
sub.clinical <- clinical.rnaseq[idinc,] 

testID1 = colnames(fpkm)[331:533]  ##  test set 1 is 203 samples from new batch
testID2 = setdiff(colnames(fpkm)[1:330], sub.clinical$w_mrn) ##  test set 2 is non-training samples from old batch
trainID = setdiff(colnames(fpkm)[1:330], testID2)
## length(trainID)          ## training set: 195 samples
## length(testID1)          ## testing set 1: 203 samples
## length(testID2)          ## testing set 2: 135 samples

dataPath = '/export/home/xurren/WTCProject/Data/'
save(trainID, testID1, testID2, file = paste(dataPath, 'TrainID_and_TestID.RData', sep = ''))







