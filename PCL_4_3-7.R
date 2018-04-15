# Analysis of the variables in the response.
library(ggplot2)
library(ggpubr)
library(remMap)
library(glmnet)
library(grid)
library(gridExtra)
library(gplots)
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
load(file ='/export/home/xurren/WTCProject/Data/train_test_ID_12_14_2017.RData')
####################################################################################


#Create 10 sets of train and test data with size ratio 7:3
nloop = 10
listofplot_remMap = list()
listofplot_Elastic = list()
listofbeta_remMap = list()
listofbeta_Elastic = list()
test_list <- list()
for (l in 1:nloop){
  test_list[[l]] = sample(nrow(X.m1), floor(nrow(X.m1)*0.3))
}
for (i in 1:nloop) 
{
  testIndexes <- test_list[[i]]
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
  
  rownames(phi.m) = colnames(exon.norm.log_1st)
  # get the fitted values for test data
  Y.fitted1 = as.data.frame(X.testData%*%phi.m)
  colnames(Y.fitted1) = listofn
  Y.testData = as.data.frame(Y.testData)
  
  listofcor <-  c()
  listofplot <-  list()
  # Store the plots and beta's in two lists for ith test
  for (j in 1:length(listofn)) {
    
    listofbeta[[j]] = list(names(phi.m[,j][phi.m[,j]!=0]))
    listofplot[[j]] = getplot(main = paste(listofn[j], 'fitted vs true by remMap----',i), 
                              xlab = paste0(listofn[j], ' Predicted') ,ylab = paste0(j, ' True'), Y.fitted1[[j]], Y.testData[[j]]) 
  }
  
  listofplot_remMap[[i]] <-  listofplot
  listofbeta_remMap[[i]] <-  listofbeta
  ##### examine the performance when using the test data.
  
  
  ##########################################################
  #Use the Elastic Net method
  trainX_data = as.matrix(X.trainData)
  listofplot <-  list()
  listofbeta <-  list()
  
  for (j in 1:length(listofn)) {
    trainY_data = as.matrix(Y.trainData[,j])
    alpha_cand = seq(0, 1, 0.01) ## candidate alpha values
    cvms = rep(0, length(alpha_cand)) ## cvms stores the minimum cross validation rate for each alpha
    cvElaNet = lapply(alpha_cand, function(a){ cv.glmnet(x = trainX_data, y = trainY_data, family = "gaussian",
                                                         nfolds = 5, alpha = a, type.measure = "mae") } )
    for(k in 1:length(alpha_cand)){
      
      cvms[k] = min(cvElaNet[[k]]$cvm)
      
    }
    opt_alpha = alpha_cand[which.min(cvms)]
    opt_lambda = cvElaNet[[which.min(cvms)]]$lambda.1se
    opt_model = cvElaNet[[which.min(cvms)]]
    cvElaNet_pred = predict(opt_model, newx = X.testData, s = opt_lambda, type = "response")
    
    ###store plots in a list
    listofplot[[j]] <- getplot(main = paste(listofn[j], 'fitted vs true by ElasticNet--',i), 
                               xlab = paste0(listofn[j], ' Predicted') ,ylab = paste0(listofn[j], ' True'), cvElaNet_pred, Y.testData[,j])
    
    ### get the names of genes with non-zero coefficients respectively for Elastic Net
    listofbeta[[j]] = list(names(coef(opt_model)[coef(opt_model)[,1]!=0,]))
    
  }
  
  listofplot_Elastic[[i]] = listofplot
  listofbeta_Elastic[[i]] = listofbeta
}

###############################################
#the summary of the analysis for the comparison
pdf(file = "/export/home/chche/WTC_Project/Results/PCL_H&A&R&N_10*3:7.pdf")
textplot(paste0('Plots and summary of 10 tests.', '\n',
                'Data is randomly splitted into test set and train set,\n, with the ratio of 3:7.'), halign = 'left', valign = 'top', cex = 0.7)

for (j in 1:length(listofn)){
  plots_remMap = list()
  beta_remMap = list() # list: store the name of non-zero variables for response j in the ith position
  plots_Elastic = list()
  beta_Elastic = list()
  
  num_remMap = c()
  num_remMap_Elastic = c()
  #Plot the remMap plots and get a list of all remMap nonzero terms
  
  for (i in 1:nloop){
    #store the jth response variable's plot into a list for all 10 tests
    plots_remMap[[i]] = listofplot_remMap[[i]][[j]]
    if (!is.null(unlist(listofbeta_remMap[[i]][[j]]))) { #discard the only intercept model
      beta_remMap <- c(beta_remMap, listofbeta_remMap[[i]][[j]])
    } else {beta_remMap = beta_remMap}
    
  }
  num_remMap <- sapply(beta_remMap, length)
  multiplot(plotlist = plots_remMap[1:5], cols = 2)
  multiplot(plotlist = plots_remMap[6:10], cols = 2)
  
  #summary text of remMap
  App_rem_10 = Reduce(intersect, beta_remMap) # the intesection of all 10 tests' nonzero terms
  textplot(paste0("Average number of nonzero coefficients in the remMap model for response ", listofn[j], ' is ',mean(num_remMap),'\n',
                  #'Average pearson correlation for 10 tests is', mean(cor_remMap), '\n',
                  "When using remMap model, ", length(App_rem_10), " of variables appear in all 10 tests", '\n',
                  "They are ", paste(App_rem_10, collapse = ' ')), halign = 'left', valign = 'top', cex = 0.7)
  #####################################################################################
  #Plot the ElasticNet plots and get a list of all Elastic net nonzero terms
  
  
  for (i in 1:nloop) {
    #store the jth response variable's plot into a list for all 10 tests
    plots_Elastic[[i]] = listofplot_Elastic[[i]][[j]]
    if (!is.null(unlist(listofbeta_Elastic[[i]][[j]]))) { #discard the only intercept model
      beta_Elastic <- c(beta_Elastic, listofbeta_Elastic[[i]][[j]])
    } else {beta_Elastic = beta_Elastic}
    #cor_remMap = c(cor_remMap, get(paste0('Elastic_cor_', j,'.',i)))
  }
  multiplot(plotlist = plots_Elastic[1:5], cols = 2)
  multiplot(plotlist = plots_Elastic[6:10], cols = 2)
  
  par(mfrow = c(1,1))
  # summary text of Elastic Net
  App_net_10 = Reduce(intersect, beta_Elastic)
  textplot(paste0("Average number of nonzero coefficients in the Elastic net model for response ", listofn[j], ' is ',mean(sapply(beta_Elastic, length)),'\n',
                  #'Average pearson correlation for 10 tests is', mean(cor_Elastic), '\n',
                  "When using Elastic Net model, ", length(App_net_10), " of variables appear in all 10 tests", "\n",
                  "They are ", paste(App_net_10, collapse = ' ')), halign = 'left', valign = 'top', cex = 0.7)
  
  
  
}
dev.off()






