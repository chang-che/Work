##  check outliers using top deviance genes and 203 new samples 

library(ggplot2)
library(grid)
library(gridExtra)

dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_533Samples/'
load(file=paste(dataPathRNASeq,'FPKM533.RData',sep=''))
load(file = '/export/home/xurren/WTCProject/Data/ExonCountsNormalize533.RData')  ## load normalized exon counts
# load(file = '/export/home/xurren/WTCProject/Data/GeneCountsNormalize533.RData')
##  dim(genecounts_normalized)
##  dim(fpkm)

genecounts_normalized = exoncounts_normalized

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

###################################################################
##  PCA on FPKM using top deviance genes and 203 new samples 
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

##  PC plots to check outliers
resultsPath = '/export/home/xurren/WTCProject/Results/'
pdf(file = paste(resultsPath, "PC_plots_new_samples_only_top_genes.pdf", sep = ''), width=10, height=10)
p1_label = ggplot(df1, aes(x=PC1, y=PC2, label = colnames(genecounts_normalized_new))) + geom_point() + geom_text(size = 2, aes(label = colnames(genecounts_normalized_new)), hjust = 0.5, vjust = -0.5)
p2_label = ggplot(df2, aes(x=PC1, y=PC3, label = colnames(genecounts_normalized_new))) + geom_point() + geom_text(size = 2, aes(label = colnames(genecounts_normalized_new)), hjust = 0.5, vjust = -0.5)
p3_label = ggplot(df3, aes(x=PC2, y=PC3, label = colnames(genecounts_normalized_new))) + geom_point() + geom_text(size = 2, aes(label = colnames(genecounts_normalized_new)), hjust = 0.5, vjust = -0.5)
# grid.arrange(p1_label, p2_label, p3_label, ncol=2, top = "PC plots for normalized counts with 203 new samples only")
grid.arrange(p1_label, p2_label, p3_label, ncol=2, top = "PC plots for normalized exon counts with 203 new samples only")

p1_FPKM = ggplot(df1_FPKM, aes(x=PC1, y=PC2, label = colnames(FPKM_new))) + geom_point() + geom_text(size = 2, aes(label = colnames(FPKM_new)), hjust = 0.5, vjust = -0.5)
p2_FPKM = ggplot(df2_FPKM, aes(x=PC1, y=PC3, label = colnames(FPKM_new))) + geom_point() + geom_text(size = 2, aes(label = colnames(FPKM_new)), hjust = 0.5, vjust = -0.5)
p3_FPKM = ggplot(df3_FPKM, aes(x=PC2, y=PC3, label = colnames(FPKM_new))) + geom_point() + geom_text(size = 2, aes(label = colnames(FPKM_new)), hjust = 0.5, vjust = -0.5)
grid.arrange(p1_FPKM, p2_FPKM, p3_FPKM, ncol=2, top = "PC plots for FPKM with 203 new samples only")
dev.off()

