library(ggplot2)
library(ggpubr)
#library(remMap)
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
# load genes based on batch and response variable
load(file = '/export/home/xurren/WTCProject/Data/match_sign_genes_list.Rdata')
load(file = '/export/home/xurren/WTCProject/Data/match_sign_cand_genes.Rdata')
listofgenes <- list()
listofgenes[[1]] <- match_sign_genes_list
listofgenes[[2]] <- list(match_sign_cand_genes, match_sign_cand_genes)
## log transform normalized exoncounts
exon.norm.log = log2(exoncounts_normalized+1)  # genes by sample
MDD533 <- clinical[,"MDD_hist_3grp"]
MDD533[MDD533<=2] = 1
MDD533[MDD533==3] = 0


##################################################
### count matrix from mostafavi
dds <- estimateSizeFactors(MDD922)
dds <- counts(dds, normalized = T)





########################################
########################################
PHQ_phrase <- c('Log transformed PHQ9.', 'Original PHQ9.')
Male_phrase <- c('Only male are tested for Mostafavi.', 'Both gender are tested.')
glist_phrase <- c('Top 2000 genes by 2-step filtering.', '40 candidate genes by 2-step filtering.')

respond <- c('MDD', 'PHQ9')
batch_num <- c(3,4)
genes_retained_2 <- list(NA, NA)
for (PHQ_ind in 1:2){
  
  
  if (PHQ_ind == 1){
    M_P533 <- cbind(as.factor(MDD533),log(clinical$PHQ9+1))
    colnames(M_P533) <- c('MDD','PHQ9')
    M.Ytest <- data.frame(as.factor(MDDstat$MDD.status-1), log(Clinical_v$Phq_Tot+1))
  }else{
    M_P533 <- cbind(as.factor(MDD533),clinical$PHQ9)
    colnames(M_P533) <- c('MDD','PHQ9')
    M.Ytest <- data.frame(as.factor(MDDstat$MDD.status-1), Clinical_v$Phq_Tot)
  }
  for (male_ind in 1:1){
    if (male_ind == 1){
      dds.log <- log2(dds+1)
      dds.log <- dds.log[, bio_cov$Sex == 0]
      M.Ytest <- M.Ytest[bio_cov$Sex == 0,]
      colnames(M.Ytest) <- c('MDD', 'PHQ9')
    }else{
      dds.log <- log2(dds+1)
      M.Ytest <- data.frame(as.factor(MDDstat$MDD.status-1), Clinical_v$Phq_Tot)
    }
    for (glist_ind in 1:1){ # 1:Top 2000 genes by 2-step filtering   2: 40 candidate genes by 2-step filtering.
      auc_list <- list(list(NA, NA), list(NA, NA))
      plot_list <- list(list(NA, NA), list(NA, NA))
      
      for (batch_ind in 1:1){
        
        
        for (j in 2:2){ # MDD or PHQ9 as the response variable
          ######### innest loop
          ##################################################
          batch_gene_filter <- listofgenes[[glist_ind]][[batch_ind]][[j]]
          
          exon.norm.log.mad <- apply(exon.norm.log[batch_gene_filter,], 2, function(x) (x-median(x))/mad(x))
          exon.norm.log.mad = t(exon.norm.log.mad) # samples by genes
          repsonse <- respond[j]
          ind = which(clinical$batch==batch_num[batch_ind] & clinical$Gender==0 & !is.na(M_P533[,respond[j]]))
          X.m = exon.norm.log.mad[ind,]
          X.m <- as.matrix(X.m)
          Y.m <- M_P533[ind, respond[j]]
          Y.m <- as.matrix(Y.m)
          #### Two Test data sets
          ind = which(clinical$batch==batch_num[-batch_ind] & clinical$Gender==0 & !is.na(M_P533[,respond[j]]))
          X.testData <-  exon.norm.log.mad[ind,]
          Y.testData <- M_P533[ind, ]
          # scale over the samples
          
          M.Xtest <- apply(dds.log[batch_gene_filter, ], 2, function(x) (x-median(x))/mad(x)) 
          #M.Xtest <- scale(dds.log)
          M.Xtest <- t(M.Xtest) #samples by genes
          #No need to scale by gene
          
          # no NA values in the M.Ytest
          Mostafavi.Xm <- M.Xtest
          Mostafavi.Ym <- M.Ytest
          
          
          ## Elastic net
          alpha_cand = seq(0, 1, 0.01) ## candidate alpha values
          cvms = rep(0, length(alpha_cand)) ## cvms stores the minimum cross validation rate for each alpha
          if (j == 1) { 
            cvElaNet = lapply(alpha_cand, function(a){ cv.glmnet(x = X.m, y = Y.m, family = "binomial",
                                                                 nfolds = 10, alpha = a, type.measure = "auc", parallel = T) } )
          }
          # PHQ9
          if (j == 2){ 
            cvElaNet = lapply(alpha_cand, function(a){ cv.glmnet(x = X.m, y = Y.m, family = "gaussian",
                                                                 nfolds = 10, alpha = a, type.measure = "mae", parallel = T) } )
          }
          for(k in 1:length(alpha_cand)){
            
            cvms[k] = min(cvElaNet[[k]]$cvm)
            
          }
          opt_alpha = alpha_cand[which.min(cvms)]
          opt_lambda = cvElaNet[[which.min(cvms)]]$lambda.1se
          opt_model = cvElaNet[[which.min(cvms)]]
          
          genes_retained_2[[PHQ_ind]] <- list(rownames(coef(opt_model)[coef(opt_model)[,1]!=0,]))
          
          
          ############ end of innest loop
        }
      }
      ####
      
      
    }
  }
}
save(genes_retained_2, file = '/export/home/chche/WTC_Project/Data/genes_retained_2.Rdata')




