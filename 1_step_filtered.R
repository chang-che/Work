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

####################################################################################################
load(file = '/export/home/chche/WTC_Project/Data/topBatch3.Rdata')
load(file = '/export/home/chche/WTC_Project/Data/topBatch4.Rdata')
listofgenes <- list(list(topBatch3_MDD, topBatch3_PHQ9), list(topBatch4_MDD, topBatch4_PHQ9))
## log transform normalized exoncounts
exon.norm.log = log2(exoncounts_normalized+1)  # genes by sample

#scale over the samples
#exon.norm.log.mad <- apply(exon.norm.log, 2, function(x) (x-median(x))/mad(x)) 
exon.norm.log.mad <- scale(exon.norm.log)
exon.norm.log.mad = t(exon.norm.log.mad) # samples by genes

MDD533 <- clinical[,"MDD_hist_3grp"]

MDD533[MDD533<=2] = 1
MDD533[MDD533==3] = 0
M_P533 <- cbind(as.factor(MDD533),clinical$PHQ9)
colnames(M_P533) <- c('MDD','PHQ9')
##################################################
### count matrix from mostafavi
dds <- estimateSizeFactors(MDD922)
dds <- counts(dds, normalized = T)
dds.log <- log2(dds+1)
dds.log <- dds.log[, bio_cov$Sex == 0] # only 274 male samples 

# scale over the samples
#dds.log <- apply(dds.log, 2, function(x) (x-median(x))/mad(x)) 
M.Xtest <- scale(dds.log)
M.Xtest <- t(M.Xtest) #samples by genes
#No need to scale by genes
M.Ytest <- data.frame(as.factor(MDDstat$MDD.status-1), Clinical_v$Phq_Tot)
M.Ytest <- M.Ytest[bio_cov$Sex == 0,]
colnames(M.Ytest) <- c('MDD', 'PHQ9') # no NA values in the M.Ytest

# filtered gene list
#genelist <- candList <- c("GPR15","LONP2","HEXA","MARC1","SLC4A10","TMED8","RFWD3","KIF7","GPR162","VPRBP","DSN1","RETSAT","MAT2A", "ANO5","CCDC66","CDH2","CORIN","DBI","DESI2","ENKD1","EPS8","ERICH1-AS1","FANCL","FEM1C","FGD6","HDX","KANK1","MAPK14","MIR5047","MTOR","NADK2","PDXP","PFN4","PPT2-EGFL8","SEC24A","SLC35F1","SMAD1","ST14","ZNF100","ZNF493")
#genelist.inter <- intersect(genelist, colnames(M.Xtest))

respond <- c('MDD','PHQ9')
batch_num <- c(3, 4)
pdf(file = "/export/home/chche/WTC_Project/Results/med_mad.pdf")
textplot(paste('PHQ9 MDD prediction with log transformation \n','and scale by samples :(x-mean(x)/sd(x)'))

for(i in 1:2){
  
  batch_i <- batch_num[i]
  for(j in 1:2){
    batch_gene_filter <- listofgenes[[i]][[j]]
    repsonse <- respond[j]
    ind = which(clinical$batch==batch_i & clinical$Gender==0 & !is.na(M_P533[,respond[j]]))
    X.m = exon.norm.log.mad[ind, batch_gene_filter]
    X.m <- as.matrix(X.m)
    Y.m <- M_P533[ind, respond[j]]
    Y.m <- as.matrix(Y.m)
    #### Two Test data sets
    ind = which(clinical$batch==batch_num[-i] & clinical$Gender==0 & !is.na(M_P533[,respond[j]]))
    X.testData <-  exon.norm.log.mad[ind,batch_gene_filter]
    Y.testData <- M_P533[ind, ]
    Mostafavi.Xm <- M.Xtest[,batch_gene_filter]
    Mostafavi.Ym <- M.Ytest
    
    ##################################################
    ## Elastic net
    alpha_cand = seq(0, 1, 0.01) ## candidate alpha values
    cvms = rep(0, length(alpha_cand)) ## cvms stores the minimum cross validation rate for each alpha
    # MDD
    if (j == 1) { 
      cvElaNet = lapply(alpha_cand, function(a){ cv.glmnet(x = X.m, y = Y.m, family = "binomial",
                                                           nfolds = 10, alpha = a, type.measure = "auc", parallel = TRUE) } )
    }
    # PHQ9
    if (j == 2){ 
      cvElaNet = lapply(alpha_cand, function(a){ cv.glmnet(x = X.m, y = Y.m, family = "gaussian",
                                                           nfolds = 10, alpha = a, type.measure = "mae", parallel = TRUE) } )
    }
    for(k in 1:length(alpha_cand)){
      
      cvms[k] = min(cvElaNet[[k]]$cvm)
      
    }
    opt_alpha = alpha_cand[which.min(cvms)]
    opt_lambda = cvElaNet[[which.min(cvms)]]$lambda.1se
    opt_model = cvElaNet[[which.min(cvms)]]
    cvElaNet_pred_WTC = predict(opt_model, newx = X.testData, s = opt_lambda, type = "response")
    
    cvElaNet_pred_Mostafari = predict(opt_model, newx = Mostafavi.Xm, s = opt_lambda, type = "response")
    
    ###plot and auc
    a <- rep(NA, 2)
    listofplot <- list()
    a[1] <- auc(roc(Y.testData[,1],as.vector(cvElaNet_pred_WTC), direction = "<"))
    a[2] <- auc(roc(Mostafavi.Ym[,1],as.vector(cvElaNet_pred_Mostafari), direction = "<"))
    
    listofplot[[1]] <- getplot(main = paste(respond[j], 'Predicted score vs PHQ9 by ElasticNet'), 
                               xlab = paste(respond[j], 'Predicted score') ,ylab = paste0('PHQ9 in WTC--', 'True'), cvElaNet_pred_WTC, Y.testData[,2])
    listofplot[[2]] <- getplot(main = paste(respond[j], 'Predicted score vs PHQ9 by ElasticNet'), 
                               xlab = paste(respond[j], 'Predicted score') ,ylab = paste0('PHQ9 in Mostafari--', 'True'), cvElaNet_pred_Mostafari, Mostafavi.Ym[,2])
    par(mfrow = c(1,1))
    textplot(paste('Train on', respond[j], 'of batch', batch_i, '\n',
                   'Test on batch', batch_num[-i], 'and Mostafavi', '\n',
                   'Prediction Model on MDD of Batch',batch_num[-i], 'and Mostafavi:', '\n',
                   'AUC are', paste(a, collapse = ' '), 'respectively'), halign = 'left', valign = 'top', cex = 1.0)
    
    par(mfrow = c(1,1))
    multiplot(plotlist = listofplot, cols = 2)
  }
    
}
dev.off()




pdf(file = "/export/home/chche/WTC_Project/Results/med_mad.pdf")
#textplot(paste0('Plots and summary of 5 tests.', '\n',
#                'Data is randomly splitted into test set and train set,\n, with the ratio of 3:7.\n This is partI. \n From 1st to 3rd.'), halign = 'left', valign = 'top', cex = 0.7)

for (j in 1:length(listofn)){
  plots_remMap = list()
  beta_remMap = list() # list: store the name of non-zero variables for response j in the ith position
  plots_Elastic = list()
  beta_Elastic = list()
  cor_remMap = c()
  num_remMap = c()
  num_remMap_Elastic = c()
  #Plot the remMap plots and get a list of all remMap nonzero terms
  
  for (i in 1:nloop){
    #store the jth response variable's plot into a list for all 10 tests
    plots_remMap[[i]] = listofplot_remMap[[i]][[j]]
    if (!is.null(unlist(listofbeta_remMap[[i]][[j]]))) { #discard the only intercept model
      beta_remMap <- list(beta_remMap, listofbeta_remMap[[i]][[j]])
    } else {beta_remMap = beta_remMap}
    
  }
  num_remMap <- sapply(beta_remMap, length)
  multiplot(plotlist = plots_remMap, cols = 2)
  
  
  #summary text of remMap
  App_rem_10 = Reduce(intersect, beta_remMap) # the intesection of all 10 tests' nonzero terms
  textplot(paste0("Average number of nonzero coefficients in the remMap model for response ", listofn[j], ' is ',mean(num_remMap),'\n',
                  #'Average pearson correlation for 10 tests is', mean(cor_remMap), '\n',
                  "When using remMap model, ", length(App_rem_10), " of variables appear in all 10 tests", '\n',
                  "(Notice that only intercept models are deleted,\n because in that case the intersection of all 10 models will always be empty.)\n",
                  "They are ", paste(App_rem_10, collapse = ' ','\n')),
           halign = 'left', valign = 'top', cex = 0.7)
  
}
par(mfrow = c(1,1))
textplot(paste(listoflambdas), halign = 'left', valign = 'top', cex = 1.0)
textplot(paste("The range of lambda1 and lambda2 are",lamL1.string,'\n', lamL2.string,'\n'), halign = 'left', valign = 'top', cex = 0.8)
dev.off()

##################################################
j = 2
X.trainData <- X.m3
Y.trainData <- Y.m3

Y.testData <- Y.m4
ind4 <- which(clinical$batch==4 & clinical$Gender==0)
X.testData = exon.norm.log[ind4, colnames(X.trainData)]
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
listofplot[[1]] <- getplot(main = paste('MDD', 'Predicted score vs PHQ9 by ElasticNet'), 
                           xlab = 'MDD Predicted score' ,ylab = paste0('PHQ9 in WTC--', 'True'), cvElaNet_pred_WTC, Y.testData[,2])
listofplot[[2]] <- getplot(main = paste('MDD', 'Predicted score vs PHQ9 by ElasticNet'), 
                           xlab = 'MDD Predicted score' ,ylab = paste0('PHQ9 in Mostafari--', 'True'), cvElaNet_pred_Mostafari, M.Ytest[,2])
ela_auc
plot(listofplot[[1]])


#############################################################################################################
X.trainData <- X.m4
Y.trainData <- Y.m4

Y.testData <- Y.m3
ind3 <- which(clinical$batch==3 & clinical$Gender==0)
X.testData = exon.norm.log[ind3, colnames(X.trainData)]
X.testData <- scale(X.testData)
#X.testData = exon.norm.log[ind3, ]
##########################################################
#Use the Elastic Net method
trainX_data = as.matrix(X.trainData)
trainY_data = as.matrix(Y.trainData[,j])
ind = !is.na(trainY_data)
trainX_data = trainX_data[ind,]
trainY_data = trainY_data[ind,]
M.X = M.Xtest[,colnames(X.trainData)]
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
listofplot[[1]] <- getplot(main = paste('MDD', 'Predicted score vs PHQ9 by ElasticNet'), 
                           xlab = 'MDD Predicted score' ,ylab = paste0('PHQ9 in WTC--', 'True'), cvElaNet_pred_WTC, Y.testData[,2])
listofplot[[2]] <- getplot(main = paste('MDD', 'Predicted score vs PHQ9 by ElasticNet'), 
                           xlab = 'MDD Predicted score' ,ylab = paste0('PHQ9 in Mostafari--', 'True'), cvElaNet_pred_Mostafari, M.Ytest[,2])
ela_auc
plot(listofplot[[1]])  


