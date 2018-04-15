##  make PCA plots color by PTSD 

library(ggplot2)
library(grid)
library(gridExtra)

dataPath = '/export/home/xurren/WTCProject/Data/'
load(file = paste(dataPath, 'TopDevianceGene.RData', sep = ''))
load(file = paste(dataPath, 'clinical533_reordered.RData', sep = ''))
load(file = paste(dataPath, 'TopDevianceGeneFPKM.RData', sep = ''))
load(file = paste(dataPath, 'FPKM533.RData',sep=''))
load(file = '/export/home/xurren/WTCProject/Data/GeneCountsNormalize533.RData')

PTSD = clinical$PTSD_SCID_3grp
##  convert PTSD to a factor
PTSD = as.factor(PTSD)
levels(PTSD) = c("current", "past", "control")

##  do PCA on topDevianceGene
prinCompAnalysis = prcomp(t(topDevianceGene), retx = TRUE, center = TRUE, scale = TRUE)
pcVariance = prinCompAnalysis$sdev^2
pcProportion = pcVariance/sum(pcVariance)
pcProportion
cumsum(pcProportion)
pc1 = prinCompAnalysis$x[, 1]
pc2 = prinCompAnalysis$x[, 2]
pc3 = prinCompAnalysis$x[, 3]

df1 = data.frame(PC1=pc1, PC2=pc2)
df2 = data.frame(PC1=pc1, PC3=pc3)
df3 = data.frame(PC2=pc2, PC3=pc3)

resultsPath = '/export/home/xurren/WTCProject/Results/'
pdf(file = paste(resultsPath, "PC_plots_colored_by_PTSD.pdf", sep = ''), width=10, height=10)
p1_PTSD = ggplot(df1, aes(x=PC1, y=PC2, color = PTSD, label = colnames(topDevianceGene))) + geom_point() + geom_text(size = 2, aes(label = colnames(topDevianceGene)), hjust = 0.5, vjust = -0.5)
p2_PTSD = ggplot(df2, aes(x=PC1, y=PC3, color = PTSD, label = colnames(topDevianceGene))) + geom_point() + geom_text(size = 2, aes(label = colnames(topDevianceGene)), hjust = 0.5, vjust = -0.5)
p3_PTSD = ggplot(df3, aes(x=PC2, y=PC3, color = PTSD, label = colnames(topDevianceGene))) + geom_point() + geom_text(size = 2, aes(label = colnames(topDevianceGene)), hjust = 0.5, vjust = -0.5)
grid.arrange(p1_PTSD, p2_PTSD, p3_PTSD, ncol=2, top = "PC plots for normalized counts colored by PTSD")


##  do PCA on topDevianceGene_fpkm
pcAnalysis_FPKM = prcomp(t(topDevianceGene_FPKM), retx = TRUE, center = TRUE, scale = TRUE)
pcVariance_FPKM = pcAnalysis_FPKM$sdev^2
pcProportion_FPKM = pcVariance_FPKM/sum(pcVariance_FPKM)
pcProportion_FPKM
cumsum(pcProportion_FPKM)
pc1_FPKM = pcAnalysis_FPKM$x[, 1]
pc2_FPKM = pcAnalysis_FPKM$x[, 2]
pc3_FPKM = pcAnalysis_FPKM$x[, 3]

df1_FPKM = data.frame(PC1=pc1_FPKM, PC2=pc2_FPKM)
df2_FPKM = data.frame(PC1=pc1_FPKM, PC3=pc3_FPKM)
df3_FPKM = data.frame(PC2=pc2_FPKM, PC3=pc3_FPKM)

p1_FPKM_PTSD = ggplot(df1_FPKM, aes(x=PC1, y=PC2, color = PTSD, label = colnames(topDevianceGene_FPKM))) + geom_point() + geom_text(size = 2, aes(label = colnames(topDevianceGene_FPKM)), hjust = 0.5, vjust = -0.5)
p2_FPKM_PTSD = ggplot(df2_FPKM, aes(x=PC1, y=PC3, color = PTSD, label = colnames(topDevianceGene_FPKM))) + geom_point() + geom_text(size = 2, aes(label = colnames(topDevianceGene_FPKM)), hjust = 0.5, vjust = -0.5)
p3_FPKM_PTSD = ggplot(df3_FPKM, aes(x=PC2, y=PC3, color = PTSD, label = colnames(topDevianceGene_FPKM))) + geom_point() + geom_text(size = 2, aes(label = colnames(topDevianceGene_FPKM)), hjust = 0.5, vjust = -0.5)
grid.arrange(p1_FPKM_PTSD, p2_FPKM_PTSD, p3_FPKM_PTSD, ncol=2, top = "PC plots for FPKM colored by PTSD")

####################################################################################################
##  do PCA on new samples only  

PTSD_new = PTSD[331:533]
##  get only new samples
genecounts_normalized_new = genecounts_normalized[, 331:533]
# nonZeroRows = ( apply(genecounts_normalized_new, 1, sum)!=0 )
# genecounts_normalized_new = genecounts_normalized_new [nonZeroRows, ]
# dim(genecounts_normalized_new)

##  select top deviance genes
medianAbsoluteDeviation = apply(genecounts_normalized_new, 1, mad)
##  length(medianAbsoluteDeviation)
##  order the medianAbsoluteDeviation
order_ind = order(medianAbsoluteDeviation, decreasing = TRUE)
##  select 200 top deviance genes
Num = 200
topDevianceGene = genecounts_normalized_new[order_ind, ]
topDevianceGene = topDevianceGene[1:Num, ]

prinCompAnalysis = prcomp(t(topDevianceGene), retx = TRUE, center = TRUE, scale = TRUE)
pcVariance = prinCompAnalysis$sdev^2
pcProportion = pcVariance/sum(pcVariance)
cumsum(pcProportion)
pc1 = prinCompAnalysis$x[, 1]
pc2 = prinCompAnalysis$x[, 2]
pc3 = prinCompAnalysis$x[, 3]

df1 = data.frame(PC1=pc1, PC2=pc2)
df2 = data.frame(PC1=pc1, PC3=pc3)
df3 = data.frame(PC2=pc2, PC3=pc3)

p1_PTSD = ggplot(df1, aes(x=PC1, y=PC2, color = PTSD_new)) + geom_point() + geom_text(size = 2, aes(label = colnames(genecounts_normalized_new)), hjust = 0.5, vjust = -0.5)
p2_PTSD = ggplot(df2, aes(x=PC1, y=PC3, color = PTSD_new)) + geom_point() + geom_text(size = 2, aes(label = colnames(genecounts_normalized_new)), hjust = 0.5, vjust = -0.5)
p3_PTSD = ggplot(df3, aes(x=PC2, y=PC3, color = PTSD_new)) + geom_point() + geom_text(size = 2, aes(label = colnames(genecounts_normalized_new)), hjust = 0.5, vjust = -0.5)
grid.arrange(p1_PTSD, p2_PTSD, p3_PTSD, ncol=2, top = "PC plots for normalized counts 203 new samples colored by PTSD")


FPKM_new = fpkm[, 331:533]
# nonZeroRows_FPKM = ( apply(FPKM_new, 1, sum)!=0 )
# FPKM_new = FPKM_new [nonZeroRows_FPKM, ]
# dim(FPKM_new)

mad_FPKM = apply(FPKM_new, 1, mad)
order_ind_FPKM = order(mad_FPKM, decreasing = TRUE)
##  select 200 top deviance genes
Num = 200
topDevianceGene_FPKM = FPKM_new[order_ind_FPKM, ]
topDevianceGene_FPKM = topDevianceGene_FPKM[1:Num, ]

pcAnalysis_FPKM = prcomp(t(topDevianceGene_FPKM), retx = TRUE, center = TRUE, scale = TRUE)
pcVariance_FPKM = pcAnalysis_FPKM$sdev^2
pcProportion_FPKM = pcVariance_FPKM/sum(pcVariance_FPKM)
cumsum(pcProportion_FPKM)
pc1_FPKM = pcAnalysis_FPKM$x[, 1]
pc2_FPKM = pcAnalysis_FPKM$x[, 2]
pc3_FPKM = pcAnalysis_FPKM$x[, 3]

df1_FPKM = data.frame(PC1=pc1_FPKM, PC2=pc2_FPKM)
df2_FPKM = data.frame(PC1=pc1_FPKM, PC3=pc3_FPKM)
df3_FPKM = data.frame(PC2=pc2_FPKM, PC3=pc3_FPKM)

p1_FPKM_PTSD = ggplot(df1_FPKM, aes(x=PC1, y=PC2, color = PTSD_new)) + geom_point() + geom_text(size = 2, aes(label = colnames(FPKM_new)), hjust = 0.5, vjust = -0.5)
p2_FPKM_PTSD = ggplot(df2_FPKM, aes(x=PC1, y=PC3, color = PTSD_new)) + geom_point() + geom_text(size = 2, aes(label = colnames(FPKM_new)), hjust = 0.5, vjust = -0.5)
p3_FPKM_PTSD = ggplot(df3_FPKM, aes(x=PC2, y=PC3, color = PTSD_new)) + geom_point() + geom_text(size = 2, aes(label = colnames(FPKM_new)), hjust = 0.5, vjust = -0.5)
grid.arrange(p1_FPKM_PTSD, p2_FPKM_PTSD, p3_FPKM_PTSD, ncol=2, top = "PC plots for FPKM 203 new samples colored by PTSD")
dev.off()



