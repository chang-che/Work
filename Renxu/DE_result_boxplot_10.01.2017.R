

library(ggplot2)
library(grid)
library(gridExtra)

dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_RemovedDuplicates/'
load(file="/export/home/xurren/WTCProject/Data/edgeRfilter_exoncounts_330_WB.RData")   ## load filtered genes
load(file=paste(dataPathRNASeq,'4celltypeRNASeq.RData',sep=''))   ##  load cell type exoncounts data

PTSD_30 = as.factor(dat$clinical[, 208])
levels(PTSD_30) = c("Never", "Current")

makeBoxPlot = function(countMat, typeName, index){
  count = t(countMat[index, ])
  result = data.frame(count, PTSD_30)
  colnames(result) = c("counts", "PTSD_status")
  plot1 = ggplot(result, aes(x = PTSD_status, y = counts)) + geom_boxplot(aes(color = PTSD_30)) + ggtitle(typeName) + theme(axis.title.x=element_blank(), plot.title = element_text(hjust = 0.5))
  return(plot1)
}

##  FKBP5:
KFBP5_ind = which(dat$genelist == "FKBP5")
p1 = makeBoxPlot(countMat = dat$CD4Tmat, typeName = "CD4T", index = KFBP5_ind)
p2 = makeBoxPlot(countMat = dat$CD8Tmat, typeName = "CD8T", index = KFBP5_ind)
p3 = makeBoxPlot(countMat = dat$CD19Bmat, typeName = "Bcell", index = KFBP5_ind)
p4 = makeBoxPlot(countMat = dat$monomat, typeName = "Mono", index = KFBP5_ind)
p5 = makeBoxPlot(countMat = dat$unsortmat, typeName = "PBMC", index = KFBP5_ind)
p6 = makeBoxPlot(countMat = dat$old.genecounts, typeName = "Whole_Blood", index = KFBP5_ind)

pdf(file = '/export/home/xurren/WTCProject/Results/boxplots_for_celltype_FKBP5.pdf')
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2, top = "Boxplot of FKBP5 counts for each cell type")
dev.off()


##  NDUFA1:
NDUFA1_ind = which(dat$genelist == "NDUFA1")
p1 = makeBoxPlot(countMat = dat$CD4Tmat, typeName = "CD4T", index = NDUFA1_ind)
p2 = makeBoxPlot(countMat = dat$CD8Tmat, typeName = "CD8T", index = NDUFA1_ind)
p3 = makeBoxPlot(countMat = dat$CD19Bmat, typeName = "Bcell", index = NDUFA1_ind)
p4 = makeBoxPlot(countMat = dat$monomat, typeName = "Mono", index = NDUFA1_ind)
p5 = makeBoxPlot(countMat = dat$unsortmat, typeName = "PBMC", index = NDUFA1_ind)
p6 = makeBoxPlot(countMat = dat$old.genecounts, typeName = "Whole_Blood", index = NDUFA1_ind)

pdf(file = '/export/home/xurren/WTCProject/Results/boxplots_for_celltype_NDUFA1.pdf')
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2, top = "Boxplot of NDUFA1 counts for each cell type")
dev.off()


# CCDC85B:
CCDC85B_ind = which(dat$genelist == "CCDC85B")
p1 = makeBoxPlot(countMat = dat$CD4Tmat, typeName = "CD4T", index = CCDC85B_ind)
p2 = makeBoxPlot(countMat = dat$CD8Tmat, typeName = "CD8T", index = CCDC85B_ind)
p3 = makeBoxPlot(countMat = dat$CD19Bmat, typeName = "Bcell", index = CCDC85B_ind)
p4 = makeBoxPlot(countMat = dat$monomat, typeName = "Mono", index = CCDC85B_ind)
p5 = makeBoxPlot(countMat = dat$unsortmat, typeName = "PBMC", index = CCDC85B_ind)
p6 = makeBoxPlot(countMat = dat$old.genecounts, typeName = "Whole_Blood", index = CCDC85B_ind)

pdf(file = '/export/home/xurren/WTCProject/Results/boxplots_for_celltype_CCDC85B.pdf')
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=2, top = "Boxplot of CCDC85B counts for each cell type")
dev.off()











