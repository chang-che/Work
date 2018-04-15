##  principal component analysis on WTCProjectGeneCountsNormalize533.RData and WTCProjectFPKM533.RData

library(ggplot2)
library(grid)
library(gridExtra)

dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/'
load(file=paste(dataPathRNASeq,'FPKM533.RData',sep=''))
load(file = '/export/home/xurren/WTCProject/Data/ExonCountsNormalize533.RData')   ## load exon counts 
# load(file = '/export/home/xurren/WTCProject/Data/GeneCountsNormalize533.RData')
##  dim(genecounts_normalized)
##  dim(fpkm)
batch = c(rep("old", 330), rep("new", 203))

genecounts_normalized = exoncounts_normalized

medianAbsoluteDeviation = apply(genecounts_normalized, 1, mad)
##  length(medianAbsoluteDeviation)
##  order the medianAbsoluteDeviation
order_ind = order(medianAbsoluteDeviation, decreasing = TRUE)
##  select 200 top deviance genes
Num = 200
topDevianceGene = genecounts_normalized[order_ind, ]
topDevianceGene = topDevianceGene[1:Num, ]
##  apply(topDevianceGene, 1, mad)
##  dim(topDevianceGene)
dataPath = '/export/home/xurren/WTCProject/Data/'
save(topDevianceGene, file = paste(dataPath, 'TopDevianceGene.RData', sep = ''))

##  do PCA on topDevianceGene
prinCompAnalysis = prcomp(t(topDevianceGene), retx = TRUE, center = TRUE, scale = TRUE)
pcVariance = prinCompAnalysis$sdev^2
pcProportion = pcVariance/sum(pcVariance)
pcProportion
cumsum(pcProportion)
pc1 = prinCompAnalysis$x[, 1]
pc2 = prinCompAnalysis$x[, 2]
pc3 = prinCompAnalysis$x[, 3]

resultsPath = '/export/home/xurren/WTCProject/Results/'
# pdf(file = paste(resultsPath, "PC_plots_for_normalized_counts.pdf", sep = ''),  width=10, height=10)
pdf(file = paste(resultsPath, "PC_plots_for_normalized_exon_counts.pdf", sep = ''),  width=10, height=10)
df1 = data.frame(PC1=pc1, PC2=pc2)
p1 = ggplot(df1, aes(x=PC1, y=PC2, color = batch)) + geom_point()
df2 = data.frame(PC1=pc1, PC3=pc3)
p2 = ggplot(df2, aes(x=PC1, y=PC3, color = batch)) + geom_point()
df3 = data.frame(PC2=pc2, PC3=pc3)
p3 = ggplot(df3, aes(x=PC2, y=PC3, color = batch)) + geom_point()
# grid.arrange(p1, p2, p3, ncol=2, top = "PC plots for normalized counts")
grid.arrange(p1, p2, p3, ncol=2, top = "PC plots for normalized exon counts")
dev.off()

#######################################################################################################

##  PCA plot on FPKM
mad_FPKM = apply(fpkm, 1, mad)
order_ind_FPKM = order(mad_FPKM, decreasing = TRUE)
##  select 200 top deviance genes
Num = 200
topDevianceGene_FPKM = fpkm[order_ind_FPKM, ]
topDevianceGene_FPKM = topDevianceGene_FPKM[1:Num, ]
save(topDevianceGene_FPKM, file = paste(dataPath, 'TopDevianceGeneFPKM.RData', sep = ''))

##  do PCA on topDevianceGene_fpkm
pcAnalysis_FPKM = prcomp(t(topDevianceGene_FPKM), retx = TRUE, center = TRUE, scale = TRUE)
pcVariance_FPKM = pcAnalysis_FPKM$sdev^2
pcProportion_FPKM = pcVariance_FPKM/sum(pcVariance_FPKM)
pcProportion_FPKM
cumsum(pcProportion_FPKM)
pc1_FPKM = pcAnalysis_FPKM$x[, 1]
pc2_FPKM = pcAnalysis_FPKM$x[, 2]
pc3_FPKM = pcAnalysis_FPKM$x[, 3]

resultsPath = '/export/home/xurren/WTCProject/Results/'
pdf(file = paste(resultsPath, "PC_plots_for_FPKM.pdf", sep = ''), width=10, height=10)
df1_FPKM = data.frame(PC1=pc1_FPKM, PC2=pc2_FPKM)
p1_FPKM = ggplot(df1_FPKM, aes(x=PC1, y=PC2, color = batch)) + geom_point()
df2_FPKM = data.frame(PC1=pc1_FPKM, PC3=pc3_FPKM)
p2_FPKM = ggplot(df2_FPKM, aes(x=PC1, y=PC3, color = batch)) + geom_point()
df3_FPKM = data.frame(PC2=pc2_FPKM, PC3=pc3_FPKM)
p3_FPKM = ggplot(df3_FPKM, aes(x=PC2, y=PC3, color = batch)) + geom_point()
grid.arrange(p1_FPKM, p2_FPKM, p3_FPKM, ncol=2, top = "PC plots FPKM")
dev.off()

##############################################################################################################################

##  check outliers
resultsPath = '/export/home/xurren/WTCProject/Results/'
pdf(file = paste(resultsPath, "PC_plots_check_outliers.pdf", sep = ''), width=10, height=10)
p1_label = ggplot(df1, aes(x=PC1, y=PC2, color = batch, label = colnames(topDevianceGene))) + geom_point() + geom_text(size = 2, aes(label = colnames(topDevianceGene)), hjust = 0.5, vjust = -0.5)
p2_label = ggplot(df2, aes(x=PC1, y=PC3, color = batch, label = colnames(topDevianceGene))) + geom_point() + geom_text(size = 2, aes(label = colnames(topDevianceGene)), hjust = 0.5, vjust = -0.5)
p3_label = ggplot(df3, aes(x=PC2, y=PC3, color = batch, label = colnames(topDevianceGene))) + geom_point() + geom_text(size = 2, aes(label = colnames(topDevianceGene)), hjust = 0.5, vjust = -0.5)
# grid.arrange(p1_label, p2_label, p3_label, ncol=2, top = "PC plots for normalized counts with sample labels")
grid.arrange(p1_label, p2_label, p3_label, ncol=2, top = "PC plots for normalized exon counts with sample labels")
p1_FPKM = ggplot(df1_FPKM, aes(x=PC1, y=PC2, color = batch, label = colnames(topDevianceGene_FPKM))) + geom_point() + geom_text(size = 2, aes(label = colnames(topDevianceGene_FPKM)), hjust = 0.5, vjust = -0.5)
p2_FPKM = ggplot(df2_FPKM, aes(x=PC1, y=PC3, color = batch, label = colnames(topDevianceGene_FPKM))) + geom_point() + geom_text(size = 2, aes(label = colnames(topDevianceGene_FPKM)), hjust = 0.5, vjust = -0.5)
p3_FPKM = ggplot(df3_FPKM, aes(x=PC2, y=PC3, color = batch, label = colnames(topDevianceGene_FPKM))) + geom_point() + geom_text(size = 2, aes(label = colnames(topDevianceGene_FPKM)), hjust = 0.5, vjust = -0.5)
grid.arrange(p1_FPKM, p2_FPKM, p3_FPKM, ncol=2, top = "PC plots FPKM with sample labels")
dev.off()

