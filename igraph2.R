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

getigraph2 = function(PCLHgenes, PCLAgenes, PCLRgenes, PCLNgenes){
  phenotypes = c(rep("PCL_H", length(PCLHgenes)), rep("PCL_A", length(PCLAgenes)), rep("PCL_R", length(PCLRgenes)), rep("PCL_N", length(PCLNgenes)))
  nonzerogenes = c(PCLHgenes, PCLAgenes, PCLRgenes, PCLNgenes)
  el = cbind(phenotypes, nonzerogenes)
  graph = graph_from_edgelist(el, directed = F)
  
  myvenn = venn(list(PCL_H = PCLHgenes, PCL_A = PCLAgenes, PCL_R = PCLRgenes, PCL_N = PCLNgenes), show.plot=F)
  intersect.list = attr(myvenn, "intersections")
  intersect.names = names(attr(myvenn, "intersections"))
  mycolors = sample(colors(), length(intersect.names))
  
  ## this is to set vertex sizes
  vsizes = rep(4, length(vertex_attr(graph)$name))
  vsizes[vertex_attr(graph)$name %in% c("PCL_H", "PCL_A", "PCL_R", "PCL_N")] = 10
  
  ## this is to set vertex colors
  vcolors = rep(NA, length(vertex_attr(graph)$name))
  for(i in 1:length(intersect.names)){
    vcolors[vertex_attr(graph)$name %in% intersect.list[[i]] ] = mycolors[i]
  }
  
  vcolors[vertex_attr(graph)$name == "PCL_H"] = mycolors[which(intersect.names == "PCL_H")]     
  vcolors[vertex_attr(graph)$name == "PCL_A"] = mycolors[which(intersect.names == "PCL_A")]   
  vcolors[vertex_attr(graph)$name == "PCL_R"] = mycolors[which(intersect.names == "PCL_R")] 
  vcolors[vertex_attr(graph)$name == "PCL_N"] = mycolors[which(intersect.names == "PCL_N")] 
  
  V(graph)$color = vcolors
  plot(graph, vertex.size = vsizes, vertex.color = vcolors, vertex.label.cex = 0.5)
  legend("topleft", legend=intersect.names, col=mycolors, pch=21, pt.bg = mycolors)
  
}


load(file = '/export/home/xurren/WTCProject/Data/clinical533_02_17_2018.RData')  ## load clinical data
load(file="/export/home/xurren/WTCProject/Data/ExonCountsNormalize533.RData")   ## load normalized exoncounts
load(file = '/export/home/xurren/WTCProject/Data/exon330_genefilter.RData') 
load(file ='/export/home/xurren/WTCProject/Data/train_test_ID_12_14_2017.RData')

# trying remMap on "PCL_H", "PCL_A", "PCL_R","PCL_N"
listofn = c("PCL_H", "PCL_A", "PCL_R","PCL_N")
exon.norm.log = log2(exoncounts_normalized+1)  ## log transform normalized exoncounts
exon.norm.log = exon.norm.log[keep,]
exon.norm.log = t(exon.norm.log)

# delete the missing value and standardize both the response and predictor.
ind = which(clinical$batch==3 & clinical$Gender==0 & !is.na(clinical$PCL_H) & !is.na(clinical$PCL_A) & !is.na(clinical$PCL_R) & !is.na(clinical$PCL_N))
exon.norm.log_1st = exon.norm.log[ind,]
X.m1 = scale(exon.norm.log_1st)

Y.m1 = cbind(clinical$PCL_H, clinical$PCL_A, clinical$PCL_R, clinical$PCL_N)
colnames(Y.m1) = listofn
Y.m1 = Y.m1[ind,]
Y.m1 = scale(Y.m1)

i = 1


listofplot_remMap = list()
listofplot_Elastic = list()
listofbeta_remMap = list()
listofbeta_Elastic = list()
test_list <- list()
listoflambdas6 = list()

lamL1.v = exp(seq(log(10),log(20), length=3)) 
lamL2.v = seq(0,5, length=3)

list_4  = remMap.CV(X.m1, Y.m1, lamL1.v, lamL2.v, C.m=NULL, fold=10, seed=1)
pick=which.min(as.vector(list_4$ols.cv))
lamL1.pick=list_4$l.index[1,pick]  ##find the optimal (LamL1,LamL2) based on the cv score 
lamL2.pick=list_4$l.index[2,pick]

listoflambdas6[[i]] = list(lamL1.pick, lamL2.pick)
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

#############


trainX_data = as.matrix(X.m1)
Y.trainData = as.matrix(Y.m1)
listofplot <-  list()
listofbeta <-  list()

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
  
  ###store plots in a list
  listofplot[[j]] <- getplot(main = paste(listofn[j], 'fitted vs true by ElasticNet--',i), 
                             xlab = paste0(listofn[j], ' Predicted') ,ylab = paste0(listofn[j], ' True'), cvElaNet_pred, Y.testData[,j])
  
  ### get the names of genes with non-zero coefficients respectively for Elastic Net
  listofbeta[[j]] = list(names(coef(opt_model)[coef(opt_model)[,1]!=0,]))
  
}

listofplot_Elastic[[i]] = listofplot
listofbeta_Elastic[[i]] = listofbeta

getigraph2(unlist(listofbeta_Elastic[[1]][[1]]), unlist(listofbeta_Elastic[[1]][[2]]), unlist(listofbeta_Elastic[[1]][[3]]), unlist(listofbeta_Elastic[[1]][[4]]))

pdf(file = "/export/home/chche/WTC_Project/Results/igraph2.pdf")
textplot(paste0('Response variables are PCL_H, PCL_A, PCL_R,PCL_N', '\n',
                'Full Data is used\n',
                'First plot is remMap, second plot is Elastic Net'), halign = 'left', valign = 'top', cex = 0.9)
getigraph(unlist(listofbeta_remMap[[1]][[1]]), unlist(listofbeta_remMap[[1]][[2]]), unlist(listofbeta_remMap[[1]][[3]]))
getigraph(unlist(listofbeta_Elastic[[1]][[1]]), unlist(listofbeta_Elastic[[1]][[2]]), unlist(listofbeta_Elastic[[1]][[3]]))
dev.off()