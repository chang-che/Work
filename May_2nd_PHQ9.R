library(ggplot2)
library(ggpubr)
library(remMap)
library(glmnet)
library(grid)
library(gridExtra)
library(gplots)
library(parallel)
library(DESeq2)
library("BiocParallel")
library(pROC)

library(doMC)
registerDoMC(cores=6)
#rm(list = ls())

getplot = function(main, xlab, ylab, pred, truelabel){
  setTitle = main
  df <- data.frame(pred=as.vector(pred),true=truelabel)
  p1 = ggplot(df, aes(x=pred, y=true)) +
    geom_point(shape=19,size=1) +   # Use hollow circles
    geom_smooth(method=lm)+ggtitle(setTitle)+ 
    theme(text = element_text(size=10),legend.text=element_text(size=10),legend.title=element_blank()) +
    xlab(xlab) +
    ylab(ylab) +
    theme(plot.title = element_text(hjust = 0.5))+stat_cor(method = "pearson")
  return(p1)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
load(file = '/export/home/chche/WTC_Project/Data/multiplot.Rdata')


load(file = '/export/home/chche/WTC_Project/Data/clinical533_02_17_2018.RData')  ## load clinical data
load(file="/export/home/chche/WTC_Project/Data/ExonCountsNormalize533.RData")   ## load normalized exoncounts
load(file = '/export/home/chche/WTC_Project/Data/exon330_genefilter.RData') 
load(file ='/export/home/chche/WTC_Project/Data/train_test_ID_12_14_2017.RData')

load(file ='/export/home/chche/WTC_Project/Data/MDD922_withoutMono.Rdata')
bio_cov <- read.table('/export/home/pfkuan/WTCproject/Epigenetics/Data/levinsonNIMH_MDD_EQTL_Covariates/covariates/Biological_and_hidden_factors.txt', header = T, sep ='\t')
MDDstat <- read.table('/export/home/pfkuan/WTCproject/Epigenetics/Data/levinsonNIMH_MDD_EQTL_Covariates/data_used_for_MDD_analysis/Dx_Case_status.txt', header = T, sep ='\t')
Clinical_v = read.table('/export/home/pfkuan/WTCproject/Epigenetics/Data/levinsonNIMH_MDD_EQTL_Covariates/data_used_for_MDD_analysis/Clinical_variables.txt', header = T, sep ='\t')
exon.norm.log = log2(exoncounts_normalized+1)  ## log transform normalized exoncounts

exon.norm.log = t(exon.norm.log)


dds <- estimateSizeFactors(MDD922)
dds <- counts(dds, normalized = T)
dds.log <- log2(dds+1)
dds.log <- dds.log[, bio_cov$Sex == 0] # exclude female samples




M.Xtest <- t(dds.log) #normalize the gene counts
M.Xtest <- scale(M.Xtest)
M.Ytest <- data.frame(as.factor(MDDstat$MDD.status-1), Clinical_v$Phq_Tot)
M.Ytest <- M.Ytest[bio_cov$Sex == 0,]
colnames(M.Ytest) <- c('MDD', 'PHQ9')


load(file = '/export/home/xurren/WTCProject/Data/XmYmPHQ9.RData')
Y.m3$MDD <- as.factor(Y.m3$MDD)
Y.m4$MDD <- as.factor(Y.m4$MDD)

genes.inter = intersect(colnames(X.m3), colnames(M.Xtest)) # has already intersected 

listofn <- c('MDD','PHQ9')


j = 2

# train on batch 3, test on batch 4
X.trainData <- X.m3
Y.trainData <- Y.m3

Y.testData <- Y.m4
ind4 <- which(clinical$batch==4 & clinical$Gender==0)
X.testData = exon.norm.log[ind4, colnames(X.trainData)]
X.testData <- scale(X.testData)
M.X = M.Xtest[,colnames(X.trainData)]
##########################################################
#Use the Elastic Net method
trainX_data = as.matrix(X.trainData)
trainY_data = as.matrix(Y.trainData[,j])
ind = !is.na(trainY_data)
trainX_data = trainX_data[ind,]
trainY_data = trainY_data[ind,]

alpha_cand = seq(0, 1, 0.01) ## candidate alpha values
cvms = rep(0, length(alpha_cand)) ## cvms stores the minimum cross validation rate for each alpha
# MDD
if (j == 1) { 
  cvElaNet = lapply(alpha_cand, function(a){ cv.glmnet(x = trainX_data, y = trainY_data, family = "binomial",
                                                       nfolds = 10, alpha = a, type.measure = "auc", parallel = TRUE) } )
}
# PHQ9
if (j == 2){ 
  cvElaNet = lapply(alpha_cand, function(a){ cv.glmnet(x = trainX_data, y = trainY_data, family = "gaussian",
                                                       nfolds = 10, alpha = a, type.measure = "mae", parallel = TRUE) } )
}

for(k in 1:length(alpha_cand)){
  
  cvms[k] = min(cvElaNet[[k]]$cvm)
  
}
opt_alpha = alpha_cand[which.min(cvms)]
opt_lambda = cvElaNet[[which.min(cvms)]]$lambda.1se
opt_model = cvElaNet[[which.min(cvms)]]

cvElaNet_pred_WTC = predict(opt_model, newx = X.testData, s = opt_lambda, type = "response")

cvElaNet_pred_Mostafari = predict(opt_model, newx = M.X, s = opt_lambda, type = "response")
###store plots in a list

listofplot <- list()
ela_auc <- c()
ela_auc[1] <- auc(roc(Y.testData[,1],as.vector(cvElaNet_pred_WTC), direction = "<"))
ela_auc[2] <- auc(roc(M.Ytest[,1],as.vector(cvElaNet_pred_Mostafari), direction = "<"))
listofplot[[1]] <- getplot(main = paste('PHQ9', 'Predicted score vs PHQ9 by ElasticNet'), 
                           xlab = 'PHQ9 Predicted score' ,ylab = paste0('PHQ9 in WTC--', 'True'), cvElaNet_pred_WTC, Y.testData[,2])
listofplot[[2]] <- getplot(main = paste('PHQ9', 'Predicted score vs PHQ9 by ElasticNet'), 
                           xlab = 'PHQ9 Predicted score' ,ylab = paste0('PHQ9 in Mostafari--', 'True'), cvElaNet_pred_Mostafari, M.Ytest[,2])
ela_auc
plot(listofplot[[1]])


#############################################################################################################
#train on batch 4, test on batch 3
X.trainData <- X.m4
Y.trainData <- Y.m4

Y.testData <- Y.m3
ind3 <- which(clinical$batch==3 & clinical$Gender==0)
X.testData = exon.norm.log[ind3, colnames(X.trainData)]
X.testData <- scale(X.testData)

##########################################################
#Use the Elastic Net method
trainX_data = as.matrix(X.trainData)
trainY_data = as.matrix(Y.trainData[,j])
ind = !is.na(trainY_data)
trainX_data = trainX_data[ind,]
trainY_data = trainY_data[ind,]
alpha_cand = seq(0, 1, 0.01) ## candidate alpha values
cvms = rep(0, length(alpha_cand)) ## cvms stores the minimum cross validation rate for each alpha
# MDD
if (j == 1) { 
  cvElaNet = lapply(alpha_cand, function(a){ cv.glmnet(x = trainX_data, y = trainY_data, family = "binomial",
                                                       nfolds = 10, alpha = a, type.measure = "auc", parallel = TRUE) } )
}
# PHQ9
if (j == 2){ 
  cvElaNet = lapply(alpha_cand, function(a){ cv.glmnet(x = trainX_data, y = trainY_data, family = "gaussian",
                                                       nfolds = 10, alpha = a, type.measure = "mae", parallel = TRUE) } )
}

for(k in 1:length(alpha_cand)){
  
  cvms[k] = min(cvElaNet[[k]]$cvm)
  
}
opt_alpha = alpha_cand[which.min(cvms)]
opt_lambda = cvElaNet[[which.min(cvms)]]$lambda.1se
opt_model = cvElaNet[[which.min(cvms)]]
cvElaNet_pred_WTC = predict(opt_model, newx = X.testData, s = opt_lambda, type = "response")
M.X = M.Xtest[,colnames(X.trainData)]
cvElaNet_pred_Mostafari = predict(opt_model, newx = M.X, s = opt_lambda, type = "response")
###store plots in a list

listofplot <- list()
ela_auc <- c()
ela_auc[1] <- auc(roc(Y.testData[,1],as.vector(cvElaNet_pred_WTC), direction = "<"))
ela_auc[2] <- auc(roc(M.Ytest[,1],as.vector(cvElaNet_pred_Mostafari), direction = "<"))
listofplot[[1]] <- getplot(main = paste('PHQ9', 'Predicted score vs PHQ9 by ElasticNet'), 
                           xlab = 'PHQ9 Predicted score' ,ylab = paste0('PHQ9 in WTC--', 'True'), cvElaNet_pred_WTC, Y.testData[,2])
listofplot[[2]] <- getplot(main = paste('PHQ9', 'Predicted score vs PHQ9 by ElasticNet'), 
                           xlab = 'PHQ9 Predicted score' ,ylab = paste0('PHQ9 in Mostafari--', 'True'), cvElaNet_pred_Mostafari, M.Ytest[,2])



