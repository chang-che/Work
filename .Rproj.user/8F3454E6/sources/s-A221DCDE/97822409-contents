# Analysis of the variables in the response.
library(ggplot2)
library(ggpubr)
library(remMap)
library(glmnet)
library(grid)
library(gridExtra)
library(gplots)
library(parallel)
library(igraph)

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

## PCLgenes is a vector contains non-zero gene names of PCL. The same with PHQ9genes and LRSgenes
getigraph = function(PCLgenes, PHQ9genes, LRSgenes){
  phenotypes = c(rep("PCL", length(PCLgenes)), rep("PHQ9", length(PHQ9genes)), rep("LRS", length(LRSgenes)))
  nonzerogenes = c(PCLgenes, PHQ9genes, LRSgenes)
  el = cbind(phenotypes, nonzerogenes)
  graph = graph_from_edgelist(el, directed = F)
  
  myvenn = venn(list(PCL = PCLgenes,PHQ9 = PHQ9genes,LRS = LRSgenes), show.plot=F)
  intersect.list = attr(myvenn, "intersections")
  intersect.names = names(attr(myvenn, "intersections"))
  mycolors = c("green", "red", "yellow", "skyblue", "pink", "brown", "cyan")
  
  ## this is to set vertex sizes
  vsizes = rep(4, length(vertex_attr(graph)$name))
  vsizes[vertex_attr(graph)$name %in% c("PCL", "PHQ9", "LRS")] = 10
  
  ## this is to set vertex colors
  vcolors = rep(NA, length(vertex_attr(graph)$name))
  for(i in 1:length(intersect.names)){
    vcolors[vertex_attr(graph)$name %in% intersect.list[[i]] ] = mycolors[i]
  }
  
  vcolors[vertex_attr(graph)$name == "PCL"] = mycolors[which(intersect.names == "PCL")]     
  vcolors[vertex_attr(graph)$name == "PHQ9"] = mycolors[which(intersect.names == "PHQ9")]   
  vcolors[vertex_attr(graph)$name == "LRS"] = mycolors[which(intersect.names == "LRS")]   
  
  V(graph)$color = vcolors
  plot(graph, vertex.size = vsizes, vertex.color = vcolors, vertex.label.cex = 0.5)
  legend("topleft", legend=intersect.names, col=mycolors, pch=21, pt.bg = mycolors)
  
}
load(file = '/export/home/xurren/WTCProject/Data/clinical533_02_17_2018.RData')  ## load clinical data
load(file="/export/home/xurren/WTCProject/Data/ExonCountsNormalize533.RData")   ## load normalized exoncounts
load(file = '/export/home/xurren/WTCProject/Data/exon330_genefilter.RData') 
load(file ='/export/home/xurren/WTCProject/Data/train_test_ID_12_14_2017.RData')
####################################################################################
# trying remMap on "PCL", "PHQ9", "LRS_total"
listofn = c("PCL", "PHQ9", "LRS_total")
exon.norm.log = log2(exoncounts_normalized+1)  ## log transform normalized exoncounts
exon.norm.log = exon.norm.log[keep,]
exon.norm.log = t(exon.norm.log)

# delete the missing value and standardize both the response and predictor.
ind = which(clinical$batch==3 & clinical$Gender==0 & !is.na(clinical$PCL) & !is.na(clinical$PHQ9) & !is.na(clinical$LRS_total))
exon.norm.log_1st = exon.norm.log[ind,]
X.m1 = scale(exon.norm.log_1st)

Y.m1 = cbind(clinical$PCL, clinical$PHQ9, clinical$LRS_total)
colnames(Y.m1) = listofn
Y.m1 = Y.m1[ind,]
Y.m1 = scale(Y.m1)

#Create 10 sets of train and test data with size ratio 7:3
nloop = 10
listofplot_remMap = list()
listofplot_Elastic = list()
listofbeta_remMap = list()
listofbeta_Elastic = list()
test_list <- list()
listoflambdas4 = list()

for (l in 1:nloop){
  test_list[[l]] = sample(nrow(X.m1), floor(nrow(X.m1)*0.3))
}

lamL1.v = exp(seq(log(10),log(20), length=30)) 
lamL2.v = seq(0,5, length=30)

remMapCV <- function(trainID){
  cv.list = remMap.CV(X.m1[trainID,], Y.m1[trainID,], lamL1.v, lamL2.v, C.m=NULL, fold=10, seed=1000)
  pick=which.min(as.vector(cv.list$ols.cv))
  lamL1.pick=cv.list$l.index[1,pick]  ##find the optimal (LamL1,LamL2) based on the cv score 
  lamL2.pick=cv.list$l.index[2,pick]
  rm(cv.list)
  return(list(lamL1.pick, lamL2.pick))
}

list_10_2lambda = mclapply(test_list, remMapCV, mc.cores = 12)

for (i in 1:nloop) 
{
  testIndexes <- test_list[[i]]
  X.testData <- X.m1[testIndexes, ]
  Y.testData <- Y.m1[testIndexes, ]
  X.trainData <- X.m1[-testIndexes, ]
  Y.trainData <- Y.m1[-testIndexes, ]
  
  ##############################################################
  #############Use remMap method
  #### use CV based on unshrinked estimator (ols.cv) 
  
  lamL1.pick= list_10_2lambda[[i]][[1]] ##find the optimal (LamL1,LamL2) based on the cv score 
  lamL2.pick= list_10_2lambda[[i]][[2]]
  
  listoflambdas4[[i]] = list(lamL1.pick, lamL2.pick)
  ##fit the remMap model under the optimal (LamL1,LamL2).
  result_1 = remMap(X.trainData, Y.trainData,lamL1=lamL1.pick, lamL2=lamL2.pick, phi0=NULL, C.m=NULL) 
  
  ## get the coefficients matrix
  phi.m = result_1$phi
  
  rownames(phi.m) = colnames(exon.norm.log_1st)
  # get the fitted values for test data
  Y.fitted1 = as.data.frame(X.testData%*%phi.m)
  colnames(Y.fitted1) = listofn
  Y.testData = as.data.frame(Y.testData)
  
  listofbeta <- list()
  listofplot <-  list()
  # Store the plots and beta's in two lists for ith test
  for (j in 1:length(listofn)) {
    
    listofbeta[[j]] = list(names(phi.m[,j][phi.m[,j]!=0]))
    listofplot[[j]] = getplot(main = paste(listofn[j], 'fitted vs true by remMap----',i), 
                              xlab = paste0(listofn[j], ' Predicted') ,ylab = paste0(j, ' True'), Y.fitted1[[j]], Y.testData[[j]]) 
  }
  
  listofplot_remMap[[i]] <-  listofplot
  listofbeta_remMap[[i]] <-  listofbeta
  
}

###############################################
#the summary of the analysis for the comparison
pdf(file = "/export/home/chche/WTC_Project/Results/test_10*3:7_v5.pdf")
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
                  "(Notice that only intercept models are deleted,\n because in that case the intersection of all 10 models will always be empty.)\n",
                  "They are ", paste(App_rem_10, collapse = ' ')), halign = 'left', valign = 'top', cex = 0.7)
  
}
dev.off()






