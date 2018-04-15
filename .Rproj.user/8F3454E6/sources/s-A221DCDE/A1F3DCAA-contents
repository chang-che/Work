##  differential expression analysis on 4 cell types: CD4T, CD8T, Mono, Bcell

library(DESeq2)
library(edgeR)

dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_RemovedDuplicates/'
load(file="/export/home/xurren/WTCProject/Data/edgeRfilter_exoncounts_330_WB.RData")   ## load filtered genes
load(file=paste(dataPathRNASeq,'4celltypeRNASeq.RData',sep=''))   ##  load cell type exoncounts data

##  colData
clinical_4cell = rbind(dat$clinical, dat$clinical, dat$clinical, dat$clinical)
PTSD_col = 208
PTSD = as.factor(clinical_4cell[, PTSD_col])
age = clinical_4cell$AgeBloodDraw
race = as.factor(clinical_4cell$race_binary)
batch = as.factor(rep(dat$batchID, 4))
cellID = c(rep('CD4T', 30), rep('CD8T', 30), rep('BCell', 30), rep('Mono', 30))
PTSDCell = paste(PTSD, cellID, sep = '_') 
colData = data.frame(PTSD, age, race, batch, cellID, PTSDCell)

## exoncounts for 4 cell types
countData = cbind(dat$CD4Tmat, dat$CD8Tmat, dat$CD19Bmat, dat$monomat)
rownames(countData) = dat$genelist
## change the colnames of countData to 1:120 so that the colnames of countData match the rownames of colData, which is required by DESeq2
## R does not accept non-unique rownames, but our sample IDs of the cell types are the same, so convert the sample IDs to 1:120
colnames(countData) = 1:120   

##  a) adjust for PTSDCell, age, race, batch
dds_a = DESeqDataSetFromMatrix(countData = countData, colData = colData, 
                                   design = ~ PTSDCell + age + race + batch)
dds_a = dds_a[keep, ]
dds_a = DESeq(dds_a)
res_a_CD4T = results(dds_a, contrast=c('PTSDCell','1_CD4T','0_CD4T'))
res_a_CD8T = results(dds_a, contrast=c('PTSDCell','1_CD8T','0_CD8T'))
res_a_BCell = results(dds_a, contrast=c('PTSDCell','1_BCell','0_BCell'))
res_a_Mono = results(dds_a, contrast=c('PTSDCell','1_Mono','0_Mono'))

dds_4cell = dds_a
res_4cell_CD4T = res_a_CD4T
res_4cell_CD8T = res_a_CD8T
res_4cell_BCell = res_a_BCell
res_4cell_Mono = res_a_Mono

##  b) do not adjust for batch, adjust for PTSDCell, age, race
dds_b = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ PTSDCell + age + race )
dds_b = dds_b[keep, ]
dds_b = DESeq(dds_b)
res_b_CD4T = results(dds_b, contrast=c('PTSDCell','1_CD4T','0_CD4T'))
res_b_CD8T = results(dds_b, contrast=c('PTSDCell','1_CD8T','0_CD8T'))
res_b_BCell = results(dds_b, contrast=c('PTSDCell','1_BCell','0_BCell'))
res_b_Mono = results(dds_b, contrast=c('PTSDCell','1_Mono','0_Mono'))

save(dds_4cell, res_4cell_CD4T, res_4cell_CD8T, res_4cell_BCell, res_4cell_Mono,
     file = "/export/home/xurren/WTCProject/Results/DE_4cell_result.RData")




