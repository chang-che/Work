## elastic net prediction, boxplot on scores 

library(glmnet)
library(ggplot2)
library(grid)
library(gridExtra)

load(file='/export/home/pfkuan/WTCproject/Epigenetics/Results/Results_Feb2017_RNASeq/plots_06Sept2017/cv.fit_DemoCellAdj_06Sept2017.RData')
load(file='/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/ExonCounts533.RData') ##  load exon counts
load(file = '/export/home/xurren/WTCProject/Data/TrainID_and_TestID.RData')  ## load train and test ID
load(file = '/export/home/xurren/WTCProject/Data/clinical533_reordered.RData')  ## load clinical data
# length(which(coef(cv.fit)!=0))

genes278 = rownames(coef(cv.fit))[-1]  ## get 278 differentially expressed genes, delete 1st element because it is intercept
# length(genes278)

exon_log = log2(exoncounts+1)   ## take log on exoncounts
exon_mad = scale(exon_log, center = apply(exon_log, 2, median), scale = apply(exon_log, 2, mad))  ##  center each column by median, scale each column by mad


## predict on 195 train data from batch 3
train195 =  t(exon_mad[genes278, trainID])
pred_train = predict(cv.fit, newx = train195, s = "lambda.1se", type = "response")

PTSD_train = clinical$PTSD_SCID_3grp[ clinical$w_mrn %in% trainID ]
##  convert PTSD to a factor
PTSD_train = as.factor(PTSD_train)
levels(PTSD_train) = c("current", "past", "control")

result_train = data.frame(pred_train, PTSD_train)
colnames(result_train) = c("riskscore", "PTSD")
pdf(file = "/export/home/xurren/WTCProject/Results/boxplot_for_riskscores_09.07.2017.pdf")
theme_update(plot.title = element_text(hjust = 0.5))
ggplot(result_train, aes(x = PTSD, y = riskscore)) + geom_boxplot(aes(color = PTSD_train)) + ggtitle("Risk score for 195 train samples in old batch")


## predict on 135 test data from batch 3
test135 = t(exon_mad[genes278, testID2]) 
pred_test135 = predict(cv.fit, newx = test135, s = "lambda.1se", type = "response")

PTSD_test135 = clinical$PTSD_SCID_3grp[ clinical$w_mrn %in% testID2 ]
##  convert PTSD to a factor
PTSD_test135 = as.factor(PTSD_test135)
levels(PTSD_test135) = c("current", "past", "control")

result_test135 = data.frame(pred_test135, PTSD_test135)
colnames(result_test135) = c("riskscore", "PTSD")
ggplot(result_test135, aes(x = PTSD, y = riskscore)) + geom_boxplot(aes(color = PTSD_test135)) + ggtitle("Risk score for 135 test samples in old batch")


## predict on 203 test data from batch 4
test203 = t(exon_mad[genes278, testID1]) 
pred_test203 = predict(cv.fit, newx = test203, s = "lambda.1se", type = "response")

PTSD_test203 = clinical$PTSD_SCID_3grp[ clinical$w_mrn %in% testID1 ]
##  convert PTSD to a factor
PTSD_test203 = as.factor(PTSD_test203)
levels(PTSD_test203) = c("current", "past", "control")

result_test203 = data.frame(pred_test203, PTSD_test203)
colnames(result_test203) = c("riskscore", "PTSD")
ggplot(result_test203, aes(x = PTSD, y = riskscore)) + geom_boxplot(aes(color = PTSD_test203)) + ggtitle("Risk score for 203 test samples in new batch")

dev.off()



