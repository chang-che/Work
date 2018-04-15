

library(glmnet)
library(ggplot2)
library(grid)
library(gridExtra)

load(file='/export/home/pfkuan/WTCproject/Epigenetics/Results/Results_Feb2017_RNASeq/cv.fit_DemoCellAdj_17Feb2017.RData')
load(file = '/export/home/xurren/WTCProject/Data/TrainID_and_TestID.RData')  ## load train and test ID
load(file = '/export/home/xurren/WTCProject/Data/GeneCountsNormalize533.RData')  ##  load normalized counts
load(file = '/export/home/xurren/WTCProject/Data/clinical533_reordered.RData')  ## load clinical data

##  center each sample by median, scale each column by mad
##  genecounts_normalized = log2(genecounts_normalized+1)
counts_madnorm = scale(genecounts_normalized, center = apply(genecounts_normalized, 2, median), scale = apply(genecounts_normalized, 2, mad))

ind = which(coef(cv.fit)!=0)  ## get 30 selected genes by glmnet
genes30 = rownames(coef(cv.fit))[ind]  
genes30 = genes30[-1]  ## delete intercept

################################################################################################################################################
##  this part is to load the PTSD status of training data
dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_RemovedDuplicates/'
load(file=paste(dataPathRNASeq,'clinical_rnaseq_30Sept2016_Train.RData',sep='')) ## load clinical.rnaseq 
load(file=paste(dataPathRNASeq,'GeneCounts_Train.RData',sep=''))  ## load genecounts

clinical.rnaseq[clinical.rnaseq<0] <- NA
clinical.rnaseq[clinical.rnaseq==5555] <- NA

f <- 208
pheno <- clinical.rnaseq[,f] ###this is the old PTSD, where 0:control/never PTSD, 1:current PTSD
phenoName <- colnames(clinical.rnaseq)[f]

idinc <- which(is.na(pheno)==FALSE&clinical.rnaseq$sex==0&clinical.rnaseq$batch==3)                                   
################################################################################################################################################

y = pheno[idinc]
trainID = as.character(trainID)
train_dat = data.frame(t(counts_madnorm[genes30, trainID]), y)
##  using 195 training samples, reestiamte the coefficients of the 30 selected genes
fit.train = glm(y~., data = train_dat, family = 'binomial')

pred_train = predict(fit.train, type='response')
match_train = match(names(pred_train), clinical$w_mrn)
PTSD_train = clinical[match_train, ]$PTSD_SCID_3grp
PTSD_train = as.factor(PTSD_train)
levels(PTSD_train) = c("current", "past", "control")
result_train = data.frame(pred_train, PTSD_train)
colnames(result_train) = c("riskscore", "PTSD")
# pdf(file = '/export/home/xurren/WTCProject/Results/boxplot log transform mad normalized.pdf')
pdf(file = '/export/home/xurren/WTCProject/Results/boxplot mad normalized.pdf')
##  make a box plot of score for train samples
ggplot(result_train, aes(x = PTSD, y = riskscore)) + geom_boxplot(aes(color = PTSD_train)) + ggtitle("Boxplot of score of train samples")


##  testID2 is 135 test data
testID2_male = clinical$w_mrn[ (clinical$w_mrn %in% testID2) & (clinical$Gender==0) ]
testID2_male = as.character(testID2_male)
test_dat_135 = data.frame(t(counts_madnorm[genes30, testID2_male]))
pred_135 = predict(fit.train, test_dat_135, type = 'response')
PTSD_135 = clinical[clinical$w_mrn %in% testID2_male, ]$PTSD_SCID_3grp
PTSD_135 = as.factor(PTSD_135)
levels(PTSD_135) = c("current", "past", "control")

result_135 = data.frame(pred_135, PTSD_135)
colnames(result_135) = c("riskscore", "PTSD")
##  make a box plot of score for 135 test samples
ggplot(result_135, aes(x = PTSD, y = riskscore)) + geom_boxplot(aes(color = PTSD_135)) + ggtitle("Boxplot of score 135 test samples")


##  testID1 is 203 test data
testID1_male = clinical$w_mrn[ (clinical$w_mrn %in% testID1) & (clinical$Gender==0) ]
testID1_male = as.character(testID1_male)
test_dat_203 = data.frame(t(counts_madnorm[genes30, testID1_male]))
pred_203 = predict(fit.train, test_dat_203, type = 'response')
PTSD_203 = clinical[clinical$w_mrn %in% testID1_male, ]$PTSD_SCID_3grp
PTSD_203 = as.factor(PTSD_203)
levels(PTSD_203) = c("current", "past", "control")

result_203 = data.frame(pred_203, PTSD_203)
colnames(result_203) = c("riskscore", "PTSD")
##  make a box plot of score for 203 test samples
ggplot(result_203, aes(x = PTSD, y = riskscore)) + geom_boxplot(aes(color = PTSD_203)) + ggtitle("Boxplot of score 203 test samples")
dev.off()








