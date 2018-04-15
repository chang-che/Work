library(ggplot2)
library(ggpubr)
library(remMap)
library(glmnet)
library(grid)
library(gridExtra)
library(gplots)
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
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


load(file = '/export/home/xurren/WTCProject/Data/clinical533_02_17_2018.RData')  ## load clinical data
load(file="/export/home/xurren/WTCProject/Data/ExonCountsNormalize533.RData")   ## load normalized exoncounts
load(file = '/export/home/xurren/WTCProject/Data/exon330_genefilter.RData') 

####################################################################################
# trying remMap on "PCL", "PHQ9", "LRS_total"
listofn = c("PCL", "PHQ9", "LRS_total")
exon.norm.log = log2(exoncounts_normalized+1)  ## log transform normalized exoncounts
exon.norm.log = exon.norm.log[keep,]
exon.norm.log = t(exon.norm.log)

# delete the missing value and standardize both the response and predictor.
ind = which(clinical$batch==3 & clinical$Gender==0 & !is.na(clinical$PCL) & !is.na(clinical$PHQ9) & !is.na(clinical$LRS_total))
exon.norm.log_1st = exon.norm.log[ind,]
X.m1 = scale(exon.norm.log_1st, center = apply(exon.norm.log_1st, 2, mean), scale = apply(exon.norm.log_1st, 2, sd))

Y.m1 = cbind(clinical$PCL, clinical$PHQ9, clinical$LRS_total)
colnames(Y.m1) = listofn
Y.m1 = Y.m1[ind,]
Y.m1 = scale(Y.m1, center = apply(Y.m1, 2, mean), scale = apply(Y.m1, 2, sd))

sam_ind = sample(nrow(X.m1))
X.m1 <- X.m1[sam_ind,]
Y.m1 <- Y.m1[sam_ind,]


#Create 5 equally size folds
folds <- cut(seq(1,nrow(Y.m1)),breaks=5,labels=FALSE)
#Perform 5 fold cross validation
for (i in 1:5) {
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
  ##############################################################
  #############Use remMap method
  #### use CV based on unshrinked estimator (ols.cv) 
  pick=which.min(as.vector(list_4$ols.cv))
  lamL1.pick=list_4$l.index[1,pick]  ##find the optimal (LamL1,LamL2) based on the cv score 
  lamL2.pick=list_4$l.index[2,pick]
  
  ##fit the remMap model under the optimal (LamL1,LamL2).
  result_1 = remMap(X.trainData, Y.trainData,lamL1=lamL1.pick, lamL2=lamL2.pick, phi0=NULL, C.m=NULL) 
  
  ## get the coefficients matrix
  phi.m = result_1$phi
  
  ########### get the names of genes with non-zero coefficients respectively
  rownames(phi.m) = colnames(exon.norm.log_1st)
  for (j in 1:length(listofn)) {
    assign(paste0('remMap_beta_', listofn[j], '.',i), names(phi.m[,j][phi.m[,j]!=0]))
  }
  
  ##### examine the performance when using the test data.
  
  
  # get the fitted values for test data
  Y.fitted1 = as.data.frame(X.testData%*%phi.m)
  colnames(Y.fitted1) = listofn
  Y.testData = as.data.frame(Y.testData)
  
  ## get the plot of "PCL", "PHQ9", "LRS_total" and assign them to the named variables
  for (j in listofn) {
  assign(paste0('remMap_', j,'.',i), getplot(main = paste(j, 'fitted vs true 5fold-cv by remMap----',i), 
                    xlab = paste0(j, ' Predictted') ,ylab = paste0(j, ' True'), Y.fitted1[[j]], Y.testData[[j]]))
    }
  
  ##########################################################
  #Use the Elastic Net method
  trainX_data = as.matrix(X.trainData)
  
  for (j in 1:length(listofn)) {
    trainY_data = as.matrix(Y.trainData[,j])
    alpha_cand = seq(0, 1, 0.01) ## candidate alpha values
    cvms = rep(0, length(alpha_cand)) ## cvms stores the minimum cross validation rate for each alpha
    cvElaNet = lapply(alpha_cand, function(a){ cv.glmnet(x = trainX_data, y = trainY_data, family = "gaussian",
                                                         nfolds = 10, alpha = a, type.measure = "mae") } )
    for(k in 1:length(alpha_cand)){
      
      cvms[k] = min(cvElaNet[[k]]$cvm)
      
    }
    opt_alpha = alpha_cand[which.min(cvms)]
    opt_lambda = cvElaNet[[which.min(cvms)]]$lambda.1se
    opt_model = cvElaNet[[which.min(cvms)]]
    cvElaNet_pred = predict(opt_model, newx = X.testData, s = opt_lambda, type = "response")
    
    #get the plot of jth variables's elastic net when ith cross-validation is implemented
    assign(paste0('ElasticNet_', listofn[j],'.',i), getplot(main = paste(listofn[j], 'fitted vs true 5fold-cv by ElasticNet----',i), 
                                                            xlab = paste0(listofn[j], ' Predictted') ,ylab = paste0(listofn[j], ' True'), cvElaNet_pred, Y.testData[,j]))
    
    ############## get the names of genes with non-zero coefficients respectively for Elastic Net

    assign(paste0('Elastic_Net_beta_', listofn[j], '.',i), names(coef(opt_model)[coef(opt_model)[,1]!=0,]))
    
    
  }
}

###############################################
#the summary of the analysis for the comparison
pdf(file = "/export/home/chche/WTC_Project/Results/PCL&PHQ9&LRS_TOTAL_remMap_vs_ElasticNet.pdf")
for (j in listofn){
  
  all_remMap = list()
  all_Elastic = list()
  num_remMap = c()
  num_remMap_Elastic = c()
  #Plot the remMap plots and get a list of all remMap nonzero terms
  plots_remMap = list()
  for (i in 1:5){
    plots_remMap[[i]] = get(paste0('remMap_', j,'.',i))
    all_remMap = c(all_remMap, list(get(paste0('remMap_beta_', j, '.', i))))
  }
  multiplot(plotlist = plots_remMap, cols = 2)
  #Plot the ElasticNet plots and get a list of all Elastic net nonzero terms
  
  plots_Elastic = list()
  for (i in 1:5){
    plots_Elastic[[i]] = get(paste0('ElasticNet_', j,'.',i))
    all_Elastic = c(all_Elastic, list(get(paste0('Elastic_Net_beta_', j, '.', i))))
  }
  multiplot(plotlist = plots_Elastic, cols = 2)
  # summary text 
  num_remMap <- sapply(all_remMap, length)
  par(mfrow = c(1,1))
  for (k in 1:5){
    num_remMap_Elastic = c(num_remMap_Elastic, length(intersect(all_remMap[[k]], all_Elastic[[k]])))
  }
  textplot(paste0("Average number of nonzero coefficients in the remMap model for response ", j, ' is ',mean(num_remMap),'\n',
                  "Average number of nonzero coefficients in the Elastic Net model for response ", j, ' is ',mean(sapply(all_Elastic, length)),'\n',
                  "About ", 100*mean(num_remMap_Elastic/num_remMap), "% of nonzero terms in remMap model is in Elastic Net model" ), halign = 'left', valign = 'top', cex = 0.6)
  
}
dev.off()


