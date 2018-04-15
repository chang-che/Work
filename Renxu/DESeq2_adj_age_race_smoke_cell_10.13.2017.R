##  run DESeq2 adjusting for age_RNA, race, and smoke_status on ONLY male subjects

library("DESeq2")
load(file = '/export/home/xurren/WTCProject/Data/clinical533_reordered.RData')  ## load clinical data
dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/'
load(file=paste(dataPathRNASeq,'ExonCounts533.RData',sep=''))
load(file = '/export/home/xurren/WTCProject/Data/TrainID_and_TestID.RData')  ## load train and test ID. testID2 is 135 test data, testID1 is 203 samples from new batch
load('/export/home/xurren/WTCProject/Results/CellProportion.RData')

genecounts = exoncounts

race1 = rep(NA, dim(clinical)[1])
race1[clinical$race8==7] = 1
race1[clinical$race8!=7] = 0

PTSD1 = rep(NA, dim(clinical)[1])
PTSD1[clinical$PTSD_hist_3grp==1] = 1
PTSD1[clinical$PTSD_hist_3grp==3] = 0

PTSD1 = factor(PTSD1)
race1 = factor(race1)
smoke_status = factor(clinical$smoker)
age = clinical$age_RNA

##  load cell proportion:
CD4T = t(proportion[5, ])
CD8T = t(proportion[10, ])
BCell = t(proportion[3, ])
Mono = t(proportion[39, ])
NK = t(proportion[41, ])


coldata = data.frame(age = age, PTSD1 = PTSD1, race1 = race1, smoke_status = smoke_status, CD4T = CD4T, CD8T = CD8T, BCell = BCell, Mono = Mono, NK = NK)
rownames(coldata) = clinical$w_mrn
colnames(coldata) = c("age", "PTSD1", "race1", "smoke_status", "CD4T", "CD8T", "BCell", "Mono", "NK")

###########################################################################################################################################
# a) on subset of 195 training data from the 330 samples
train195 = setdiff(clinical$w_mrn[1:330], testID2)  ##  train195 is 195 training data from the 330 samples 
# train195_male = which( (clinical$w_mrn %in% train195) & clinical$Gender==0 & is.na(PTSD1)==FALSE )
# train195_male = which( (clinical$w_mrn %in% train195) & clinical$Gender==0 & is.na(PTSD1)==FALSE & is.na(race1)==FALSE & is.na(smoke_status)==FALSE)
train195_male = which( (clinical$w_mrn %in% train195) & clinical$Gender==0 & is.na(PTSD1)==FALSE & is.na(race1)==FALSE)
count195_male = genecounts[, train195_male]
coldata195_male = coldata[train195_male, ]

## please note: NK for 195 samples is all 0, so I did not adjust NK
# dds195 = DESeqDataSetFromMatrix(countData = count195_male, colData = coldata195_male, design = ~ PTSD1+age+race1+smoke_status+CD4T+CD8T+BCell+Mono)
dds195 = DESeqDataSetFromMatrix(countData = count195_male, colData = coldata195_male, design = ~ PTSD1+age+race1+CD4T+CD8T+BCell+Mono)
# dds195 = DESeqDataSetFromMatrix(countData = count195_male, colData = coldata195_male, design = ~ PTSD1+age)
dds195 = DESeq(dds195)
res195 = results(dds195, contrast = c("PTSD1", "1", "0"))
lfc195 = res195$log2FoldChange


# b) on subset of 135 test data 
test135 = testID2  ##  test135 is 135 testing data from the 330 samples
# test135_male = which( (clinical$w_mrn %in% test135) & clinical$Gender==0 & is.na(PTSD1)==FALSE )
# test135_male = which( (clinical$w_mrn %in% test135) & clinical$Gender==0 & is.na(PTSD1)==FALSE & is.na(race1)==FALSE & is.na(smoke_status)==FALSE)
test135_male = which( (clinical$w_mrn %in% test135) & clinical$Gender==0 & is.na(PTSD1)==FALSE & is.na(race1)==FALSE)
count135_male = genecounts[, test135_male]
coldata135_male = coldata[test135_male, ]

## please note: all 135 test samples are non-smokers, and their NK are all 0
dds135 = DESeqDataSetFromMatrix(countData = count135_male, colData = coldata135_male, design = ~ PTSD1+age+race1+CD4T+CD8T+BCell+Mono)
# dds135 = DESeqDataSetFromMatrix(countData = count135_male, colData = coldata135_male, design = ~ PTSD1+age)
dds135 = DESeq(dds135)
res135 = results(dds135, contrast = c("PTSD1", "1", "0"))
lfc135 = res135$log2FoldChange


# c) on 330 samples 
count330_male = cbind(count195_male, count135_male)
coldata330_male = rbind(coldata195_male, coldata135_male)
## rownames(coldata330_male) == colnames(count330_male)
# please note: NK are all 0
dds330 = DESeqDataSetFromMatrix(countData = count330_male, colData = coldata330_male, design = ~ PTSD1+age+race1+CD4T+CD8T+BCell+Mono)
# dds330 = DESeqDataSetFromMatrix(countData = count330_male, colData = coldata330_male, design = ~ PTSD1+age+race1+smoke_status+CD4T+CD8T+BCell+Mono)
# dds330 = DESeqDataSetFromMatrix(countData = count330_male, colData = coldata330_male, design = ~ PTSD1+age+race1+CD4T+CD8T+BCell+Mono)
dds330 = DESeq(dds330)
res330 = results(dds330, contrast = c("PTSD1", "1", "0"))
lfc330 = res330$log2FoldChange  

#################################################################################
set.seed(100)
# randomly split the 203 samples into 102 train samples and 101 test samples
train_sample = clinical$w_mrn[sample(331:533, 102, replace = FALSE)]
train102 = clinical$w_mrn[clinical$w_mrn %in% train_sample]
test101 = setdiff(clinical$w_mrn[331:533], train102)
##################################################################################

# d) on training data from 203 samples 
# train102_male = which( (clinical$w_mrn %in% train102) & clinical$Gender==0 & is.na(PTSD1)==FALSE & is.na(race1)==FALSE & is.na(smoke_status)==FALSE)
train102_male = which( (clinical$w_mrn %in% train102) & clinical$Gender==0 & is.na(PTSD1)==FALSE & is.na(race1)==FALSE)
# train102_male = which( (clinical$w_mrn %in% train102) & clinical$Gender==0 & is.na(PTSD1)==FALSE )
count102_male = genecounts[, train102_male]
coldata102_male = coldata[train102_male, ]

dds102 = DESeqDataSetFromMatrix(countData = count102_male, colData = coldata102_male, design = ~ PTSD1+age+race1+CD4T+CD8T+BCell+Mono)
# dds102 = DESeqDataSetFromMatrix(countData = count102_male, colData = coldata102_male, design = ~ PTSD1+age+race1+smoke_status+CD4T+CD8T+BCell+Mono+NK)
# dds102 = DESeqDataSetFromMatrix(countData = count102_male, colData = coldata102_male, design = ~ PTSD1+age)
dds102 = DESeq(dds102)
res102 = results(dds102, contrast = c("PTSD1", "1", "0"))
lfc102 = res102$log2FoldChange


# e) on test data from 203 samples
test101_male = which( (clinical$w_mrn %in% test101) & clinical$Gender==0 & is.na(PTSD1)==FALSE & is.na(race1)==FALSE )
# test101_male = which( (clinical$w_mrn %in% test101) & clinical$Gender==0 & is.na(PTSD1)==FALSE & is.na(race1)==FALSE & is.na(smoke_status)==FALSE)
# test101_male = which( (clinical$w_mrn %in% test101) & clinical$Gender==0 & is.na(PTSD1)==FALSE )
count101_male = genecounts[, test101_male]
coldata101_male = coldata[test101_male, ]

dds101 = DESeqDataSetFromMatrix(countData = count101_male, colData = coldata101_male, design = ~ PTSD1+age+race1+CD4T+CD8T+BCell+Mono)
# dds101 = DESeqDataSetFromMatrix(countData = count101_male, colData = coldata101_male, design = ~ PTSD1+age+race1+smoke_status+CD4T+CD8T+BCell+Mono+NK)
# dds101 = DESeqDataSetFromMatrix(countData = count101_male, colData = coldata101_male, design = ~ PTSD1+age)
dds101 = DESeq(dds101)
res101 = results(dds101, contrast = c("PTSD1", "1", "0"))
lfc101 = res101$log2FoldChange


# f) on 203 samples
count203_male = cbind(count102_male, count101_male)
coldata203_male = rbind(coldata102_male, coldata101_male)
## rownames(coldata203_male) == colnames(count203_male)
dds203 = DESeqDataSetFromMatrix(countData = count203_male, colData = coldata203_male, design = ~ PTSD1+age+race1+CD4T+CD8T+BCell+Mono)
# dds203 = DESeqDataSetFromMatrix(countData = count203_male, colData = coldata203_male, design = ~ PTSD1+age+race1+smoke_status+CD4T+CD8T+BCell+Mono+NK)
# dds203 = DESeqDataSetFromMatrix(countData = count203_male, colData = coldata203_male, design = ~ PTSD1+age)
dds203 = DESeq(dds203)
res203 = results(dds203, contrast = c("PTSD1", "1", "0"))
lfc203 = res203$log2FoldChange 

save(lfc195, res195, lfc135, res135, lfc330, res330, 
     lfc102, res102, lfc101, res101, lfc203, res203, 
     file = '/export/home/xurren/WTCProject/Results/lfc.RData')


