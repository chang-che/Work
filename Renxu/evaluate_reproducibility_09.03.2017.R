## compare the estimated log2FoldChange from DESeq2 output

library(glmnet)
load(file = '/export/home/xurren/WTCProject/Results/lfc.RData')
load(file='/export/home/pfkuan/WTCproject/Epigenetics/Results/Results_Feb2017_RNASeq/cv.fit_DemoCellAdj_17Feb2017.RData')
dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/'
load(file=paste(dataPathRNASeq,'ExonCounts533.RData',sep=''))

genecounts = exoncounts

genes448 = rownames(coef(cv.fit))[-1]  ## get 448 genes used as input for glmnet
ind = which(coef(cv.fit)!=0)   
genes30 = rownames(coef(cv.fit))[ind]  
genes30 = genes30[-1]    ## get 30 selected genes by glmnet
genes448_ind = ( rownames(genecounts) %in% genes448 ) 
genes30_ind = ( rownames(genecounts) %in% genes30 )

# we have the log fold change from: 
# a) on subset of 195 training data from the 330 samples
# b) on subset of 135 test data 
# c) on 330 samples 
# d) on training data from 203 samples 
# e) on test data from 203 samples
# f) on 203 samples

# (a) vs (b)
# 1) all the genes
cor(lfc195, lfc135, use = "pairwise.complete.obs")
# 2) 448 genes 
cor(lfc195[genes448_ind], lfc135[genes448_ind], use = "pairwise.complete.obs")
# 3) 30 genes
cor(lfc195[genes30_ind], lfc135[genes30_ind], use = "pairwise.complete.obs")

# (d) vs (e)
# 1) all the genes
cor(lfc102, lfc101, use = "pairwise.complete.obs")
# 2) 448 genes 
cor(lfc102[genes448_ind], lfc101[genes448_ind], use = "pairwise.complete.obs")
# 3) 30 genes
cor(lfc102[genes30_ind], lfc101[genes30_ind], use = "pairwise.complete.obs")

# (c) vs (f)
# 1) all the genes
cor(lfc330, lfc203, use = "pairwise.complete.obs")
# 2) 448 genes 
cor(lfc330[genes448_ind], lfc203[genes448_ind], use = "pairwise.complete.obs")
# 3) 30 genes
cor(lfc330[genes30_ind], lfc203[genes30_ind], use = "pairwise.complete.obs")

# log2FoldChange for FKBP5
lfc195[which(rownames(genecounts)=="FKBP5")]
lfc135[which(rownames(genecounts)=="FKBP5")]
lfc330[which(rownames(genecounts)=="FKBP5")]
lfc102[which(rownames(genecounts)=="FKBP5")]
lfc101[which(rownames(genecounts)=="FKBP5")]
lfc203[which(rownames(genecounts)=="FKBP5")]

##########################################################################################################################

# sign of lfc
# (a) vs (b)
table(train = (lfc195 > 0), test = (lfc135 > 0))
# (d) vs (e)
table(train = (lfc102 > 0), test = (lfc101 > 0))
# (c) vs (f)
table(train = (lfc330 > 0), test = (lfc203 > 0))

