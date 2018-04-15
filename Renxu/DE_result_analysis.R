
## DE result analysis
library(ggplot2)
load(file = "/export/home/xurren/WTCProject/Results/DE_WB_result.RData")
load(file = "/export/home/xurren/WTCProject/Results/DE_4cell_result.RData")
load(file = "/export/home/xurren/WTCProject/Results/DE_PBMC_result.RData")


padj_cutoff = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
CD4T_num = rep(0, length(padj_cutoff))
CD8T_num = rep(0, length(padj_cutoff))
Mono_num = rep(0, length(padj_cutoff))
Bcell_num = rep(0, length(padj_cutoff))
PBMC_num = rep(0, length(padj_cutoff))
WB30_num = rep(0, length(padj_cutoff))
for (i in 1:length(padj_cutoff)){
  CD4T_num[i] = sum(res_4cell_CD4T$padj < padj_cutoff[i], na.rm = TRUE)
  CD8T_num[i] = sum(res_4cell_CD8T$padj < padj_cutoff[i], na.rm = TRUE)
  Mono_num[i] = sum(res_4cell_Mono$padj < padj_cutoff[i], na.rm = TRUE)
  Bcell_num[i] = sum(res_4cell_BCell$padj < padj_cutoff[i], na.rm = TRUE)
  PBMC_num[i] = sum(res_PBMC$padj < padj_cutoff[i], na.rm = TRUE)
  WB30_num[i] = sum(res30_WB$padj < padj_cutoff[i], na.rm = TRUE)
}
CD4T_num
CD8T_num
Mono_num
Bcell_num
PBMC_num
WB30_num

df = data.frame(celltype = rep(c("CD4T", "CD8T", "Mono", "Bcell", "PBMC", "WholeBlood"), each = length(padj_cutoff)),
                padjcutoff = rep(padj_cutoff, 6),
                num_of_DE_genes = c(CD4T_num, CD8T_num, Mono_num, Bcell_num, PBMC_num, WB30_num)
                )
ggplot(data = df, aes(x = padjcutoff, y = num_of_DE_genes, group = celltype)) + geom_line(aes(color = celltype)) + geom_point(aes(color = celltype)) + ggtitle("Number of DE genes versus padj cut off") + geom_text(data = subset(df, padjcutoff == 0.4), aes(label = celltype, color = celltype, x = Inf, y = num_of_DE_genes), hjust = 1.0) + theme(plot.title = element_text(hjust = 0.5)) 

library(gplots)
## top 200 genes
num = 200
CD4T_top = rownames(res_4cell_CD4T[order(res_4cell_CD4T$padj),])[1:num]
CD8T_top = rownames(res_4cell_CD8T[order(res_4cell_CD8T$padj),])[1:num]
Mono_top = rownames(res_4cell_Mono[order(res_4cell_Mono$padj),])[1:num]
Bcell_top = rownames(res_4cell_BCell[order(res_4cell_BCell$padj),])[1:num]
PBMC_top = rownames(res_PBMC[order(res_PBMC$padj),])[1:num]
WB30_top = rownames(res30_WB[order(res30_WB$padj),])[1:num]
WB330_top = rownames(res330_WB[order(res330_WB$padj),])[1:num]


#CD4T_top
#CD8T_top
#Mono_top
#Bcell_top
#PBMC_top
#WB30_top
#WB330_top

venn(list(CD4T = CD4T_top, CD8T = CD8T_top, Mono = Mono_top, Bcell = Bcell_top, PBMC = PBMC_top))

## compare each cell type with whole blood
par(mfrow = c(3, 2))
venn(list(CD4T = CD4T_top, Whole_Blood330 = WB330_top, Whole_Blood30 = WB30_top))
venn(list(CD8T = CD8T_top, Whole_Blood330 = WB330_top, Whole_Blood30 = WB30_top))
venn(list(Mono = Mono_top, Whole_Blood330 = WB330_top, Whole_Blood30 = WB30_top))
venn(list(Bcell = Bcell_top, Whole_Blood330 = WB330_top, Whole_Blood30 = WB30_top))
venn(list(PBMC = PBMC_top, Whole_Blood330 = WB330_top, Whole_Blood30 = WB30_top))


## one-to-one compare
par(mfrow = c(4, 3))
venn(list(CD4T = CD4T_top, CD8T = CD8T_top))
venn(list(CD4T = CD4T_top, Mono = Mono_top))
venn(list(CD4T = CD4T_top, Bcell = Bcell_top))
venn(list(CD4T = CD4T_top, PBMC = PBMC_top))
venn(list(CD8T = CD8T_top, Mono = Mono_top))
venn(list(CD8T = CD8T_top, Bcell = Bcell_top))
venn(list(CD8T = CD8T_top, PBMC = PBMC_top))
venn(list(Mono = Mono_top, Bcell = Bcell_top))
venn(list(Mono = Mono_top, PBMC = PBMC_top))
venn(list(Bcell = Bcell_top, PBMC = PBMC_top))



##  heatmap
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
count_total = cbind(counts(dds_4cell, normalize = TRUE), counts(dds_PBMC, normalize = TRUE))
count_log = log2(count_total + 1)
sampleDists = dist(t(count_log))
sampleDistMatrix = as.matrix(sampleDists)
celltype = c(rep("CD4T", 30), rep("CD8T", 30), rep("BCell", 30), rep("Mono", 30), rep("PBMC", 30))
rownames(sampleDistMatrix) = celltype
colnames(sampleDistMatrix) = celltype
colors = colorRampPalette( rev(brewer.pal(9, "RdYlBu")) )(256)

pheatmap(mat = sampleDistMatrix, clustering_distance_rows=sampleDists, 
         clustering_distance_cols=sampleDists, col=colors, fontsize_row = 3, fontsize_col = 3,
         main = 'heatmap for 5 cell types')

vsd.fast = vst(dds_4cell, blind=FALSE)
pcaData = plotPCA(vsd.fast, intgroup = c("PTSD", "cellID"), returnData = TRUE)
pdf(file = "/export/home/xurren/WTCProject/Results/PCPlot_4Celltype.pdf")
ggplot(pcaData, aes(PC1, PC2, color=PTSD, shape=cellID)) + geom_point(size=3) + coord_fixed() + ggtitle("Principal Component plot on 4 cell types") + theme(plot.title = element_text(hjust = 0.5))
dev.off()







