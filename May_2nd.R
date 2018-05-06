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
rm(list = ls())

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
MDDstat$MDD.status



  
exon.norm.log = log2(exoncounts_normalized+1)  ## log transform normalized exoncounts
exon.norm.log = exon.norm.log[keep,]
exon.norm.log = t(exon.norm.log)


listofn <- c('MDD','PHQ9') 

# delete the missing value and standardize both the response and predictor.
# !is.na(apply(clinical[listofn], 1, sum) is equivalent to !is.na(clinical$PCL) & !is.na(clinical$PHQ9) & !is.na(clinical$LRS_total)

ind3 <- which(clinical$batch==3 & clinical$Gender==0 & !is.na(apply(clinical[listofn], 1, sum)))# deleting the na's is not necessary
ind4 <- which(clinical$batch==4 & clinical$Gender==0)

exon.norm.log_1st = exon.norm.log[ind3,name.inter]

# feature engineering Mostafari data set
dds <- estimateSizeFactors(MDD922)
dds <- counts(dds, normalized = T)
dds.log <- log2(dds+1)
genes.inter = intersect(colnames(exon.norm.log), rownames(dds))
dds.log <- dds.log[genes.inter, bio_cov$Sex == 0] # exclude female samples
dds.log <- t(scale(dds.log))

exon.norm.log <- exon.norm.log[,genes.inter]
#Batch 3
X.m3 = scale(exon.norm.log_1st)
Y.m3 = clinical[ind3,listofn]




# Batch 4
exon.norm.log_1st = exon.norm.log[ind4,]
X.m4 = scale(exon.norm.log_1st)
Y.m4 = clinical[ind4,listofn]
colnames(Y.m4) = listofn


#Mostafari Dataset (used as test dataset)
M.Xtest <- 
M.Ytest <- 

nloop = 1
for (i in 1:nloop) 
{
  ## Batch 3 as train
  X.trainData <- X.m3
  Y.trainData <- Y.m3
  
  X.testData <- X.m4
  Y.testData <- Y.m4
  ##########################################################
  #Use the Elastic Net method
  trainX_data = as.matrix(X.trainData)
  
  for (j in 1:length(listofn)) {
    trainY_data = as.matrix(Y.trainData[,j])
    alpha_cand = seq(0, 1, 0.01) ## candidate alpha values
    cvms = rep(0, length(alpha_cand)) ## cvms stores the minimum cross validation rate for each alpha
    if (j == 1) { # MDD
      cvElaNet = lapply(alpha_cand, function(a){ cv.glmnet(x = trainX_data, y = trainY_data, family = "gaussian",
                                                           nfolds = 10, alpha = a, type.measure = "auc", parallel = TRUE) } )
    }
    if (j == 2){ # PHQ9
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
    cvElaNet_pred_Mostafari = predict(opt_model, newx = M.Xtest, s = opt_lambda, type = "response")
    ###store plots in a list
    
    listofplot <- list()
    ela_auc <- c()
    ela_auc[1] <- auc(cvElaNet_pred_WTC, Y.testData[1])
    ela_auc[2] <- auc(cvElaNet_pred_Mostafari, M.Ytest[1])
    listofplot[[1]] <- getplot(main = paste('MDD', 'fitted vs true by ElasticNet--',i), 
                               xlab = 'MDD Predicted score' ,ylab = paste0('PHQ9 in WTC--', 'True'), cvElaNet_pred_WTC, Y.testData[,2])
    listofplot[[2]] <- getplot(main = paste('MDD', 'fitted vs true by ElasticNet--',i), 
            xlab = 'MDD Predicted score' ,ylab = paste0('PHQ9 in Mostafari--', 'True'), cvElaNet_pred_Mostafari, M.Ytest[,2])
    
  
    }
  
  listofplot_Elastic[[i]] = listofplot
  listofbeta_Elastic[[i]] = listofbeta
  
}

###############################################
#the summary of the analysis for the comparison
pdf(file = "/export/home/chche/WTC_Project/Results/elastic_MDD_batch3_vs_4.pdf")
#textplot(paste0('Plots and summary of 5 tests.', '\n',
#                'Data is randomly splitted into test set and train set,\n, with the ratio of 3:7.'), halign = 'left', valign = 'top', cex = 0.7)

  
  
  par(mfrow = c(1,1))
  # summary text of Elastic Net
  App_net_10 = Reduce(intersect, beta_Elastic)
  textplot(paste0("Average number of nonzero coefficients in the Elastic net model for response ", listofn[j], ' is ',mean(sapply(beta_Elastic, length)),'\n',
                  #'Average pearson correlation for 10 tests is', mean(cor_Elastic), '\n',
                  "When using Elastic Net model, ", length(App_net_10), " of variables appear in all 10 tests", "\n",
                  "They are ", paste(App_net_10, collapse = ' ')), halign = 'left', valign = 'top', cex = 0.7)
  
}
dev.off()

