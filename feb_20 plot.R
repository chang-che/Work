library(ggplot2)
library(ggpubr)
library(remMap)

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
load(file = '/export/home/xurren/WTCProject/Data/clinical533_02_17_2018.RData')  ## load clinical data
load(file="/export/home/xurren/WTCProject/Data/ExonCountsNormalize533.RData")   ## load normalized exoncounts
load(file = '/export/home/xurren/WTCProject/Data/exon330_genefilter.RData') 

####################################################################################
# trying remMap on "PCL", "PHQ9", "LRS_total"

exon.norm.log = log2(exoncounts_normalized+1)  ## log transform normalized exoncounts
exon.norm.log = exon.norm.log[keep,]
exon.norm.log = t(exon.norm.log)

# delete the missing value and standardize both the response and predictor.
ind = which(clinical$batch==3 & clinical$Gender==0 & !is.na(clinical$PCL) & !is.na(clinical$PHQ9) & !is.na(clinical$LRS_total))
exon.norm.log_1st = exon.norm.log[ind,]
X.m1 = scale(exon.norm.log_1st, center = apply(exon.norm.log_1st, 2, mean), scale = apply(exon.norm.log_1st, 2, sd))

Y.m1 = cbind(clinical$PCL, clinical$PHQ9, clinical$LRS_total)
colnames(Y.m1) = c("PCL", "PHQ9", "LRS_total")
Y.m1 = Y.m1[ind,]
Y.m1 = scale(Y.m1, center = apply(Y.m1, 2, mean), scale = apply(Y.m1, 2, sd))

sam_ind = sample(nrow(X.m1))
X.m1 <- X.m1[sam_ind,]
Y.m1 <- Y.m1[sam_ind,]


#Create 5 equally size folds
folds <- cut(seq(1,nrow(Y.m1)),breaks=5,labels=FALSE)
#Perform 5 fold cross validation
for(i in 1:5){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  X.testData <- X.m1[testIndexes, ]
  Y.testData <- Y.m1[testIndexes, ]
  X.trainData <- X.m1[-testIndexes, ]
  Y.trainData <- Y.m1[-testIndexes, ]
  
  ## Tuning the parameters
  lamL1.v = exp(seq(log(10),log(20), length=3)) 
  lamL2.v = seq(0,5, length=3)
  list_4  = remMap.CV(X.trainData, Y.trainData, lamL1.v, lamL2.v, C.m=NULL, fold=10, seed=1)
  
  ############ use CV based on unshrinked estimator (ols.cv) 
  pick=which.min(as.vector(list_4$ols.cv))
  lamL1.pick=list_4$l.index[1,pick]  ##find the optimal (LamL1,LamL2) based on the cv score 
  lamL2.pick=list_4$l.index[2,pick]
  
  ##fit the remMap model under the optimal (LamL1,LamL2).
  result_1 = remMap(X.trainData, Y.trainData,lamL1=lamL1.pick, lamL2=lamL2.pick, phi0=NULL, C.m=NULL) 
  
  ## get the coefficients matrix
  phi.m = result_1$phi
  rownames(phi.m) = colnames(exon.norm.log)
  
  # get the names of genes with non-zero coefficients respectively
  rownames(phi.m) = colnames(exon.norm.log_1st)
  assign(paste0('beta_PCL.',i), names(phi.m[,1][phi.m[,1]!=0]))
  assign(paste0('beta_PHQ9.',i), names(phi.m[,2][phi.m[,2]!=0]))
  assign(paste0('beta_LRS_total.',i), names(phi.m[,3][phi.m[,3]!=0]))
  ##### examine the performance when using the test data.
  
  
  # get the fitted values for test data
  Y.fitted1 = as.data.frame(X.testData%*%phi.m)
  colnames(Y.fitted1) = c("PCL", "PHQ9", "LRS_total")
  Y.testData = as.data.frame(Y.testData)
  
  ## get the plot of "PCL", "PHQ9", "LRS_total"
  PCL = getplot(main = paste('PCL fitted vs true 5fold-cv by remMap----',i), 
                xlab = 'PCL PREDICTTED' ,ylab = 'PCL TRUE', Y.fitted1$PCL, Y.testData$PCL)
  PHQ9 = getplot(main = paste('PHQ9 fitted vs true 5fold-cv by remMap----', i),
                 xlab = 'PHQ9 PREDICTTED' ,ylab = 'PHQ9 TRUE',Y.fitted1$PHQ9, Y.testData$PHQ9)
  LRS_total = getplot(main = paste('LRS_total fitted vs true 5fold-cv by remMap----', i), xlab = 'LRS_total PREDICTTED', 
                      ylab = 'LRS_total TRUE', Y.fitted1$LRS_total, Y.testData$LRS_total)
  
  assign(paste0('PCL.',i), PCL)
  assign(paste0('PHQ9.',i), PHQ9)
  assign(paste0('LRS_total.',i), LRS_total)
  
}
for j in c("PCL", "PHQ9", "LRS_total") {
  for (i in 1:5) {
    plot(get(paste0('PCL.',i)))
  }
}

