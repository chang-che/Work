

library(DESeq2)
load(file='/export/home/xurren/WTCProject/Data/clinical533newsleep_reordered.RData')
load('/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/ExonCounts533.RData')

clinical = clinical[1:330, ] ## this file only deal with batch 3 samples
exoncounts = exoncounts[, 1:330]
exoncounts = exoncounts[, -which(clinical$smoker==1 | clinical$Gender==1)]
clinical = clinical[-which(clinical$smoker==1 | clinical$Gender==1), ]  ## remove female and smokers

race1 = rep(NA, dim(clinical)[1])
race1[clinical$race8==7] = 1
race1[clinical$race8!=7] = 0
race1 = factor(race1)

age = clinical$age_RNA

stop_breath = clinical$stop_breath
daytime_sleep = clinical$daytime_sleep
fall_stay_asleep = clinical$fall_stay_asleep
difficulty_fall_stay_asleep = clinical$difficulty_fall_stay_asleep

coldata = data.frame(age = age, race1 = race1, 
                     stop_breath = stop_breath, daytime_sleep = daytime_sleep, 
                     fall_stay_asleep = fall_stay_asleep, difficulty_fall_stay_asleep = difficulty_fall_stay_asleep)
rownames(coldata) = clinical$w_mrn

###########################################################################################################################################
## 1) DESeq on stop_breath
ind = which(stop_breath==0 | stop_breath==1)
coldata1 = coldata[ind, ]
coldata1$stop_breath = as.factor(coldata1$stop_breath)
exoncounts1 = exoncounts[, ind]

dds_stop_breath = DESeqDataSetFromMatrix(countData = exoncounts1, colData = coldata1, design = ~ stop_breath + age + race1)
dds_stop_breath = DESeq(dds_stop_breath)
res_stop_breath = results(dds_stop_breath, contrast = c("stop_breath", "1", "0"))


###########################################################################################################################################
## 2) DESeq on daytime_sleep
ind2 = which(!is.na(daytime_sleep))
coldata2 = coldata[ind2, ]
coldata2$daytime_sleep = as.factor(coldata2$daytime_sleep)
exoncounts2 = exoncounts[, ind2]

dds_daytime_sleep = DESeqDataSetFromMatrix(countData = exoncounts2, colData = coldata2, design = ~ daytime_sleep + age + race1)
dds_daytime_sleep = DESeq(dds_daytime_sleep)
res_daytime_sleep = results(dds_daytime_sleep, contrast = c("daytime_sleep", "1", "0"))


###########################################################################################################################################
## 3) DESeq on fall_stay_asleep
ind3 = which(!is.na(fall_stay_asleep))
coldata3 = coldata[ind3, ]
exoncounts3 = exoncounts[, ind3]

dds_fall_stay_asleep = DESeqDataSetFromMatrix(countData = exoncounts3, colData = coldata3, design = ~ fall_stay_asleep + age + race1)
dds_fall_stay_asleep = DESeq(dds_fall_stay_asleep)
res_fall_stay_asleep = results(dds_fall_stay_asleep, name = "fall_stay_asleep")

###########################################################################################################################################
## 4) DESeq on difficulty_fall_stay_asleep
ind4 = which(!is.na(difficulty_fall_stay_asleep))
coldata4 = coldata[ind4, ]
coldata4$difficulty_fall_stay_asleep = as.factor(coldata4$difficulty_fall_stay_asleep)
exoncounts4 = exoncounts[, ind4]

dds_difficulty_fall_stay_asleep = DESeqDataSetFromMatrix(countData = exoncounts4, colData = coldata4, design = ~ difficulty_fall_stay_asleep + age + race1)
dds_difficulty_fall_stay_asleep = DESeq(dds_difficulty_fall_stay_asleep)
res_difficulty_fall_stay_asleep = results(dds_difficulty_fall_stay_asleep, contrast = c("difficulty_fall_stay_asleep", "1", "0"))


save(res_stop_breath, res_daytime_sleep, res_fall_stay_asleep, res_difficulty_fall_stay_asleep, 
     file = '/export/home/xurren/WTCProject/Results/res_breath_sleep.RData')


