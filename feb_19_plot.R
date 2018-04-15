library(ggplot2)
library(ggpubr)
library(remMap)
library(glmnet)

elaNetModel = function(train_data, test_data, train_response, test_response, response_type, penalty_fac){
  
  train_data = as.matrix(train_data)
  
  test_data = as.matrix(test_data)
  
  alpha_cand = seq(0, 1, 0.1)  ## candidate alpha values
  
  cvms = rep(0, length(alpha_cand))  ## cvms stores the minimum cross validation rate for each alpha
  
  if(response_type == "binary"){
    
    train_response = as.factor(train_response)  ## convert the response to a factor
    
    test_response = as.factor(test_response)  ## convert the response to a factor
    
    cvElaNet = lapply(alpha_cand, function(a){ cv.glmnet(x = train_data, y = train_response, family = "binomial", nfolds = nfolds, alpha = a, penalty.factor = penalty_fac, type.measure = "auc") } )
    
  }
  
  if(response_type == "continuous")
    
    cvElaNet = lapply(alpha_cand, function(a){ cv.glmnet(x = train_data, y = train_response, family = "gaussian", nfolds = nfolds, alpha = a, penalty.factor = penalty_fac, type.measure = "mae") } )
  
  for(i in 1:length(alpha_cand)){
    
    cvms[i] = min(cvElaNet[[i]]$cvm)
    
  }
  
  opt_alpha = alpha_cand[which.min(cvms)]
  
  opt_lambda = cvElaNet[[which.min(cvms)]]$lambda.1se
  
  message(paste("In elastic model, the optimal alpha is: ", opt_alpha))
  
  message(paste("In elastic model, the optimal lambda is: ", opt_lambda))
  
  opt_model = cvElaNet[[which.min(cvms)]]
  
  cvElaNet_trainpred = predict(opt_model, newx = train_data, s = opt_lambda, type = "response")
  
  cvElaNet_testpred = predict(opt_model, newx = test_data, s = opt_lambda, type = "response")
  
  res = list(train_score = cvElaNet_trainpred, test_score = cvElaNet_testpred)
  
  return(res)
  
}
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


## Tuning the parameters
lamL1.v = exp(seq(log(10),log(20), length=3)) 
lamL2.v = seq(0,5, length=3)
list_4  = remMap.CV(X.m1, Y.m1, lamL1.v, lamL2.v, C.m=NULL, fold=10, seed=1)

############ use CV based on unshrinked estimator (ols.cv) 
pick=which.min(as.vector(list_4$ols.cv))
lamL1.pick=list_4$l.index[1,pick]  ##find the optimal (LamL1,LamL2) based on the cv score 
lamL2.pick=list_4$l.index[2,pick]

##fit the remMap model under the optimal (LamL1,LamL2).
result_1 = remMap(X.m1, Y.m1,lamL1=lamL1.pick, lamL2=lamL2.pick, phi0=NULL, C.m=NULL) 

## get the coefficients matrix
phi.m = result_1$phi

# get the names of genes with non-zero coefficients
rownames(phi.m) = colnames(exon.norm.log_1st)
names(phi.m[1,][phi.m[1,]!=0])

# get the test responce true and fitted values
Y.fitted1 = as.data.frame(X.m1%*%phi.m)
colnames(Y.fitted1) = c("PCL", "PHQ9", "LRS_total")
Y.m1 = as.data.frame(Y.m1)

## get the plot of "PCL", "PHQ9", "LRS_total"
PCL.1 = getplot(main = 'PCL fitted vs true 5fold-cv by remMap', xlab = 'PCL PREDICTTED' ,ylab = 'PCL TRUE', Y.fitted1$PCL, Y.m1$PCL)
PCL.2
PCL.3

p1.2 = getplot(main = 'PHQ9 fitted vs true 5fold-cv by remMap', xlab = 'PHQ9 PREDICTTED' ,ylab = 'PHQ9 TRUE',Y.fitted1$PHQ9, Y.m1$PHQ9)
p1.3 = getplot(main = 'LRS_total fitted vs true 5fold-cv by remMap', xlab = 'LRS_total PREDICTTED' ,ylab = 'LRS_total TRUE', Y.fitted1$LRS_total, Y.m1$LRS_total)

pdf(file = "/export/home/chche/WTC_Project/Results/PCL&PHQ9&LRS_TOTAL.pdf")

plot(p1.1)
plot(p1.2)
plot(p1.3)

dev.off()


############################################################################################################
### trying remMap on "PCL_H", "PCL_A", "PCL_R","PCL_N"

ind = which(clinical$batch==3 & clinical$Gender==0 & !is.na(clinical$PCL_H) & !is.na(clinical$PCL_A) & !is.na(clinical$PCL_R) & !is.na(clinical$PCL_N))
exon.norm.log_2nd = exon.norm.log[ind,]
X.m2 = scale(exon.norm.log_2nd, center = apply(exon.norm.log_2nd, 2, mean), scale = apply(exon.norm.log_2nd, 2, sd))
Y.m2 = cbind(clinical$PCL_H, clinical$PCL_A, clinical$PCL_R, clinical$PCL_N)
colnames(Y.m2) = c("PCL_H", "PCL_A", "PCL_R","PCL_N")
Y.m2 = Y.m2[ind,]
Y.m2 = scale(Y.m2, center = apply(Y.m2, 2, mean), scale = apply(Y.m2, 2, sd))
Y.m2 = as.data.frame(Y.m2)


#Tuning parameters
lamL1.v = exp(seq(log(10),log(20), length=3)) 
lamL2.v = seq(0,5, length=3)
list_4_2nd  = remMap.CV(X.m2, Y.m2, lamL1.v, lamL2.v, C.m=NULL, fold=10, seed=1)

############ use CV based on unshrinked estimator (ols.cv) 

pick=which.min(as.vector(list_4_2nd$ols.cv))
lamL1.pick=list_4_2nd$l.index[1,pick]  ##find the optimal (LamL1,LamL2) based on the cv score 
lamL2.pick=list_4_2nd$l.index[2,pick]

##fit the remMap model under the optimal (LamL1,LamL2).

result_2 = remMap(X.m2, Y.m2,lamL1=lamL1.pick, lamL2=lamL2.pick, phi0=NULL, C.m=NULL) 
 

## get the coefficients matrix
phi.m = result_2$phi
rownames(phi.m)
phi.m[1,] 
Y.fitted2 = as.data.frame(X.m2%*%phi.m)
colnames(Y.fitted2) = c("PCL_H", "PCL_A", "PCL_R","PCL_N")
Y.m2 = as.data.frame(Y.m2)
## Plot Four fitted vs true value
pdf(file = "/export/home/chche/WTC_Project/Results/PCL_H&PCL_A&PCL_R&PCL_N.pdf")

p2.1 = getplot(main = 'PCL_H fitted vs true', xlab = 'PCL_H PREDICTTED', ylab = 'PCL_H TRUE', Y.fitted2$PCL_H, Y.m2$PCL_H)
p2.2 = getplot(main = 'PCL_A fitted vs true', xlab = 'PCL_A PREDICTTED', ylab = 'PCL_A TRUE', Y.fitted2$PCL_A, Y.m2$PCL_A)
p2.3 = getplot(main = 'PCL_R fitted vs true', xlab = 'PCL_R PREDICTTED', ylab = 'PCL_R TRUE', Y.fitted2$PCL_R, Y.m2$PCL_R)
p2.4 = getplot(main = 'PCL_N fitted vs true', xlab = 'PCL_N PREDICTTED', ylab = 'PCL_N TRUE', Y.fitted2$PCL_N, Y.m2$PCL_N)
plot(p2.1)
plot(p2.2)
plot(p2.3)
plot(p2.4)

dev.off()

