
##  differential expression analysis on whole blood

library(DESeq2)
library(edgeR)

dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_RemovedDuplicates/'
load(file="/export/home/xurren/WTCProject/Data/edgeRfilter_exoncounts_330_WB.RData")   ## load filtered genes
load(file=paste(dataPathRNASeq,'ExonCounts.RData',sep='')) ## whole blood (processed)
load(file=paste(dataPathRNASeq,'clinical_rnaseq_07July2017.RData',sep='')) ## load clinical data
load(file=paste(dataPathRNASeq,"cell_est_rnaseq_11Feb2017_SumOne.Rdata",sep='')) ## load cell proportion

## DE analysis on 330 samples
noNA = which(is.na(clinical.rnaseq$AgeBloodDraw)==FALSE & is.na(clinical.rnaseq[, 208])==FALSE
             & is.na(clinical.rnaseq$race_binary)==FALSE)
clinical.rnaseq = clinical.rnaseq[noNA, ]
PTSD_col = 208
PTSD = as.factor(clinical.rnaseq[, PTSD_col])
age = clinical.rnaseq$AgeBloodDraw
race = as.factor(clinical.rnaseq$race_binary)
cell.est330 = cell.est.rnaseq[noNA, ]
colData = data.frame(PTSD, age, race, cell.est330)
rownames(colData) = clinical.rnaseq$w_mrn

countData = exoncounts[, noNA]  ## exoncounts for 4 whole blood data
dds330 = DESeqDataSetFromMatrix(countData = countData, colData = colData, 
                                design = ~ PTSD + age + race + CD8T + CD4T + NK + Bcell + Mono)
dds330 = dds330[keep, ]  
dds330 = DESeq(dds330)
res330 = results(dds330, contrast = c('PTSD', '1', '0'))

dds330_WB = dds330
res330_WB = res330


## DE analysis on 30 samples
load(file=paste(dataPathRNASeq,'4celltypeRNASeq.RData',sep='')) 
sampleID = colnames(dat$unsortmat)
rownames(colData) %in% sampleID 
colData30 = colData[rownames(colData) %in% sampleID, ]
countData30 = countData[, rownames(colData) %in% sampleID]
dds30 = DESeqDataSetFromMatrix(countData = countData30, colData = colData30, 
              design = ~ PTSD + age + race + CD8T + CD4T + NK + Bcell + Mono)
dds30 = dds30[keep, ]  
dds30 = DESeq(dds30)
res30 = results(dds30, contrast = c('PTSD', '1', '0'))

dds30_WB = dds30
res30_WB = res30

save(dds330_WB, res330_WB, dds30_WB, res30_WB, file = "/export/home/xurren/WTCProject/Results/DE_WB_result.RData")

























