# Analysis of the variables in the response.
library(ggplot2)
library(ggpubr)
library(remMap)
library(glmnet)
library(grid)
library(gridExtra)
library(gplots)
library(parallel)

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
load(file = '/export/home/chche/WTC_Project/Data/exon330_genefilter.RData') ## exclude low counts 
load(file ='/export/home/chche/WTC_Project/Data/train_test_ID_12_14_2017.RData')
####################################################################################
# trying remMap on "PCL", "PHQ9", "LRS_total"
#####################################################################################
######### the vector of names of response variables we want to investigate
listofn = c("PCL_H", "PCL_A", "PCL_R","PCL_N")

exon.norm.log = log2(exoncounts_normalized+1)  ## log transform normalized exoncounts
exon.norm.log = exon.norm.log[keep,]
exon.norm.log = t(exon.norm.log)

# delete the missing value and standardize both the response and predictor.
# !is.na(apply(clinical[listofn], 1, sum) is equivalent to !is.na(clinical$PCL) & !is.na(clinical$PHQ9) & !is.na(clinical$LRS_total)

ind = which(clinical$batch==3 & clinical$Gender==0 & !is.na(apply(clinical[listofn], 1, sum)))
exon.norm.log_1st = exon.norm.log[ind,]
X.m1 = scale(exon.norm.log_1st)
Y.m1 = clinical[,listofn]
colnames(Y.m1) = listofn
Y.m1 = Y.m1[ind,]
Y.m1 = scale(Y.m1)
rm(list = c('exoncounts_normalized','exon.norm.log','clinical'))
#Create 10 sets of train and test data with size ratio 7:3
nloop = 5
load(file = '/export/home/chche/WTC_Project/Data/mar_19_test_list.Rdata')
test_list <- test_list[1:nloop]
listofplot_remMap = list()
listofplot_Elastic = list()
listofbeta_remMap = list()
listofbeta_Elastic = list()

listoflambdas = list()


##################################################################################################
################################# define the range of lambda1 and lambda2!!!!!!!
lamL1.string = 'exp(seq(log(0.1),log(40), length=40))'
lamL2.string = 'seq(0,40, length=41)'
lamL1.v = eval(parse(text = lamL1.string)) 
lamL2.v = eval(parse(text = lamL2.string))
trainID <- test_list[[1]]
cv.list = remMap.CV(X.m1[trainID,], Y.m1[trainID,], lamL1.v, lamL2.v, C.m=NULL, fold=10, seed=1000)
pick=which.min(as.vector(cv.list$ols.cv))


lamL1.pick=cv.list$l.index[1,pick]  ##find the optimal (LamL1,LamL2) based on the cv score 
lamL2.pick=cv.list$l.index[2,pick]
save(cv.list, file = '/export/home/chche/WTC_Project/Data/remMapCV.Rdata')
?contour
plot(contour(x = exp(seq(log(0.1),log(40), length=40)), y =seq(0,40, length=41), z = cv.index.matrix,
             xlim = c(0.1, 40), ylim = c(0,40),zlim = range(cv.index.matrix), levels = pretty(range(cv.index.matrix), 1000)))
# Transform data to long form
ols.melt <- melt(cv.list$ols.cv, id.vars = c('l1', 'l2'), measure.vars = 'ols')
names(mtrx.melt) <- c('l1', 'l2', 'ols')
# Return data to numeric form
mtrx.melt$wt <- as.numeric(str_sub(mtrx.melt$wt, str_locate(mtrx.melt$wt, '=')[1,1] + 1))
mtrx.melt$hp <- as.numeric(str_sub(mtrx.melt$hp, str_locate(mtrx.melt$hp, '=')[1,1] + 1))

head(mtrx.melt)