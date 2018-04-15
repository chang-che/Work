##  heatmaps on WTCProjectTopDevianceGene.RData and WTCProjectTopDevianceGeneFPKM.RData

library(RColorBrewer)
library(pheatmap)

dataPath <- '/export/home/xurren/WTCProject/Data/'
load(file = paste(dataPath, 'TopDevianceGene.RData', sep = ''))
load(file = paste(dataPath, 'TopDevianceGeneFPKM.RData', sep = ''))

##  heapmap for WTCProjectTopDevianceGene.RData  
sampleDists = dist(t(topDevianceGene))
sampleDistMatrix = as.matrix(sampleDists)
batch = c(rep("old", 330), rep("new", 203))
rownames(sampleDistMatrix) = paste(colnames(topDevianceGene), batch, sep = "-")
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette( rev(brewer.pal(9, "RdYlBu")) )(256)
# pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, fontsize_row = 3, main = 'heatmap of normalized counts')

##  heatmap for WTCProjectTopDevianceGeneFPKM.RData
sampleDists_FPKM = dist(t(topDevianceGene_FPKM))
sampleDistMatrix_FPKM = as.matrix(sampleDists_FPKM)
batch = c(rep("old", 330), rep("new", 203))
rownames(sampleDistMatrix_FPKM) = paste(colnames(topDevianceGene_FPKM), batch, sep = "-")
colnames(sampleDistMatrix_FPKM) = NULL
colors = colorRampPalette( rev(brewer.pal(9, "RdYlBu")) )(256)
# pheatmap(sampleDistMatrix_FPKM, clustering_distance_rows=sampleDists_FPKM, clustering_distance_cols=sampleDists_FPKM, col=colors, fontsize_row = 3, main = 'heatmap of FPKM')

resultsPath = '/export/home/xurren/WTCProject/Results/'
pdf(file = paste(resultsPath, "heatmaps.pdf", sep = ''),  width=10, height=10)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, fontsize_row = 3, main = 'heatmap of normalized counts')
pheatmap(sampleDistMatrix_FPKM, clustering_distance_rows=sampleDists_FPKM, clustering_distance_cols=sampleDists_FPKM, col=colors, fontsize_row = 3, main = 'heatmap of FPKM')
dev.off()





