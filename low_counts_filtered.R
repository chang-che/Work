library(DESeq2)
library(edgeR)
library("BiocParallel")
register(MulticoreParam(multicoreWorkers()))
load(file = '/export/home/xurren/WTCProject/Data/clinical533_02_17_2018.RData') ## load clinical data
load(file="/export/home/xurren/WTCProject/Data/ExonCountsNormalize533.RData") ## load normalized exoncounts
load(file="/export/home/xurren/WTCProject/Data/ExonCounts533.RData")
load(file = '/export/home/xurren/WTCProject/Data/exon330_genefilter.RData')
# Mostafavi dataset
load(file ='/export/home/chche/WTC_Project/Data/MDD922_withoutMono.Rdata')
load(file = '/export/home/chche/WTC_Project/Data/MDD_int.Rdata')
MDDstat <- read.table('/export/home/pfkuan/WTCproject/Epigenetics/Data/levinsonNIMH_MDD_EQTL_Covariates/data_used_for_MDD_analysis/Dx_Case_status.txt', header = T, sep ='\t')
MDD.status <- as.factor(MDDstat$MDD.status)
#filter low counts genes in Mostafavi
countData <- t(MDD_inte)
coldata922 <- data.frame(MDD.status)
colnames(coldata922) <- c('MDD')
dds <- DESeqDataSetFromMatrix(countData = countData , colData = coldata922, design = ~ 1)
y <- DGEList(counts=countData, group = as.integer(MDD.status))
keep <- rowSums(cpm(y)>1) >= 10 ### note you can change 10 to other number
genes_mostafavi <- rownames(countData[keep,])

# filter low counts genes in WTC
MDD_WTC <- clinical[,"MDD_hist_3grp"]
MDD_WTC[MDD_WTC<=2] <- 0
MDD_WTC[MDD_WTC==3] <- 1
coldata535 <- data.frame(MDD = as.factor(MDD_WTC))
dds <- DESeqDataSetFromMatrix(countData = exoncounts , colData = coldata535, design = ~ 1)
y <- DGEList(counts=exoncounts, group = MDD_WTC)
keep <- rowSums(cpm(y)>1) >= 10 ### note you can change 10 to other number
genes_wtc <- rownames(exoncounts[keep,])




# intersection between filtered Mostafavi and filtered WTC
gene.inter <- intersect(genes_mostafavi, genes_wtc)

## train on 330 test on 203
clinical330 = clinical[1:330,]
clinical203 = clinical[331:533,]

ind = which(!is.na(clinical330$MDD_hist_3grp) & clinical330$Gender==0 & clinical330$smoker==0)
clinical330 = clinical330[ind, ]
ind = which(!is.na(clinical203$MDD_hist_3grp) & clinical203$Gender==0 & clinical203$smoker==0)
clinical203 = clinical203[ind, ]

MDD330 = clinical330$MDD_hist_3grp
MDD330[MDD330<=2] = 1
MDD330[MDD330==3] = 0

MDD203 = clinical203$MDD_hist_3grp
MDD203[MDD203<=2] = 1
MDD203[MDD203==3] = 0

coldata = data.frame(MDD = MDD330)
coldata$MDD = as.factor(coldata$MDD)
rownames(coldata) = as.character(clinical330$w_mrn)
countdata = exoncounts[gene.inter, as.character(clinical330$w_mrn)]



dds1 = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ MDD)


dds1 = DESeq(dds1, parallel = T)
res1 = results(dds1, contrast = c("MDD", "1", "0"))


res1 = res1[order(res$pvalue),]

topBatch3_MDD = rownames(res1)[1:2000]
res1_LFC <- res1$log2FoldChange[1:2000]

############################################################
coldata = data.frame(MDD = MDD203)
coldata$MDD = as.factor(coldata$MDD)
rownames(coldata) = as.character(clinical203$w_mrn)
countdata = exoncounts[topBatch3_MDD, as.character(clinical203$w_mrn)]
dds2 = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ MDD)


dds2 = DESeq(dds2, parallel = T)
res2 = results(dds2, contrast = c("MDD", "1", "0"))
res2_LFC <- res2$log2FoldChange
ind <-  ((res1_LFC>0)==(res2_LFC>0))
match_sign_genes_MDD_3 <- topBatch3_MDD[ind]

save(match_sign_genes_MDD_3, file = '/export/home/chche/WTC_Project/Data/match_sign_genes_MDD_2.Rdata')

##################################################
# PHQ9 top genes
clinical330 = clinical[1:330,]
clinical203 = clinical[331:533,]

ind = which(!is.na(clinical330$PHQ9) & clinical330$Gender==0 & clinical330$smoker==0)
clinical330 = clinical330[ind,]
ind = which(!is.na(clinical203$PHQ9) & clinical203$Gender==0 & clinical203$smoker==0)
clinical203 = clinical203[ind,]

PHQ9330 = clinical330$PHQ9
PHQ9203 = clinical203$PHQ9

coldata = data.frame(PHQ9 = PHQ9330)
rownames(coldata) = as.character(clinical330$w_mrn)
countdata = exoncounts[gene.inter,as.character(clinical330$w_mrn)]
dds = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ PHQ9)

dds = DESeq(dds, parallel = T)
res = results(dds, name = "PHQ9")


res = res[order(res$pvalue),]

topBatch3_PHQ9 = rownames(res)[1:2000]
res1_LFC <- res$log2FoldChange[1:2000]

##################################################
coldata = data.frame(PHQ9 = PHQ9203)
rownames(coldata) = as.character(clinical203$w_mrn)
countdata = exoncounts[topBatch3_PHQ9,as.character(clinical203$w_mrn)]
dds = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ PHQ9)

dds = DESeq(dds, parallel = T)
res = results(dds, name = "PHQ9")
res2_LFC <- res$log2FoldChange

match_sign_genes_PHQ9 <- topBatch3_PHQ9[(res1_LFC>0)==(res2_LFC>0)]
save(match_sign_genes_PHQ9, file = '/export/home/chche/WTC_Project/Data/match_sign_genes_PHQ9_2.Rdata')


