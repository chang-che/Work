
##  differential expression analysis on PBMC

library(DESeq2)
library(edgeR)

dataPathRNASeq <- '/export/home/pfkuan/WTCproject/Epigenetics/Data/RNASeq/ProcessedData_RemovedDuplicates/'
load(file="/export/home/xurren/WTCProject/Data/edgeRfilter_exoncounts_330_WB.RData")   ## load filtered genes
load(file=paste(dataPathRNASeq,'4celltypeRNASeq.RData',sep=''))   ##  load cell type exoncounts data
load(file=paste(dataPathRNASeq,"cell_est_rnaseq_11Feb2017_SumOne.Rdata",sep=''))
# str(cell.est.rnaseq)

## exclude Gran so that the other 5 cell types sum up to 1: 
rownames(cell.est.rnaseq) = substr(rownames(cell.est.rnaseq), 1, 6)
cell.est.5types = t( apply(cell.est.rnaseq[, 1:5], 1, function(x) x/sum(x) ) )
# save(cell.est.5types, file = "/export/home/xurren/WTCProject/Data/cell_5types.R")

##  run DESeq on 30 samples
match30 = match(colnames(dat$unsortmat), rownames(cell.est.5types))
cell.est30 = cell.est.5types[match30, ]
PTSD_col = 208
PTSD = as.factor(dat$clinical[, PTSD_col])
age = dat$clinical$AgeBloodDraw
race = as.factor(dat$clinical$race_binary)
batch = as.factor(dat$batchID)
colData = data.frame(PTSD, age, race, batch, cell.est30)

## exoncounts for PBMC
countData = dat$unsortmat
rownames(countData) = dat$genelist
dds_PBMC = DESeqDataSetFromMatrix(countData = countData, colData = colData, 
                               design = ~ PTSD + age + race + batch + CD8T + CD4T + Bcell + Mono)
dds_PBMC = dds_PBMC[keep, ]  ## filter genes
dds_PBMC = DESeq(dds_PBMC)
dds_PBMC = nbinomWaldTest(dds_PBMC, betaPrior=FALSE, maxit = 5000)
res_PBMC = results(dds_PBMC, contrast=c('PTSD','1','0'))

save(dds_PBMC, res_PBMC, file = "/export/home/xurren/WTCProject/Results/DE_PBMC_result.RData")












