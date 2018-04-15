

##  DE using 533 samples
library("DESeq2")
load(file = '/export/home/xurren/WTCProject/Data/clinical533_reordered.RData')  ## load clinical data
dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/'
load(file=paste(dataPathRNASeq,'ExonCounts533.RData',sep=''))
# load(file = '/export/home/xurren/WTCProject/Data/TrainID_and_TestID.RData')  ## load train and test ID. testID2 is 135 test data, testID1 is 203 samples from new batch
load('/export/home/xurren/WTCProject/Results/CellProportion.RData')

genecounts = exoncounts

batch = c(rep("old", 330), rep("new", 203))
batch = as.factor(batch)

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


coldata = data.frame(age = age, PTSD1 = PTSD1, race1 = race1, smoke_status = smoke_status, CD4T = CD4T, CD8T = CD8T, BCell = BCell, Mono = Mono, NK = NK, batch = batch)
rownames(coldata) = clinical$w_mrn
colnames(coldata) = c("age", "PTSD1", "race1", "smoke_status", "CD4T", "CD8T", "BCell", "Mono", "NK", "batch")

male_all = which( clinical$Gender==0 & is.na(PTSD1)==FALSE & is.na(race1)==FALSE )
count_male = genecounts[, male_all]
coldata_male = coldata[male_all, ]
dds_all = DESeqDataSetFromMatrix(countData = count_male, colData = coldata_male, design = ~ PTSD1+age+race1+CD4T+CD8T+BCell+Mono+batch)

dds_all = DESeq(dds_all)
res_all = results(dds_all, contrast = c("PTSD1", "1", "0"))

save(res_all, file = '/export/home/xurren/WTCProject/Results/DE_all.RData')





