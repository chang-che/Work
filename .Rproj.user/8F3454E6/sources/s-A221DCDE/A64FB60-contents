library(DESeq2)
library("BiocParallel")
load(file = '/export/home/chche/WTC_Project/Data/MDD_int.Rdata')
register(MulticoreParam(multicoreWorkers()))#register the number of cores in that computer
load(file = '/export/home/chche/WTC_Project/Data/MDD_int.Rdata')
bio_cov <- read.table('/export/home/pfkuan/WTCproject/Epigenetics/Data/levinsonNIMH_MDD_EQTL_Covariates/covariates/Biological_and_hidden_factors.txt', header = T, sep ='\t')
xcell <- read.table('/export/home/chche/WTC_Project/Data/xCell_xcell_xCell_1852040918.txt', header = T, sep ='\t')
MDDstat <- read.table('/export/home/pfkuan/WTCproject/Epigenetics/Data/levinsonNIMH_MDD_EQTL_Covariates/data_used_for_MDD_analysis/Dx_Case_status.txt', header = T, sep ='\t')
MDD.status <- as.factor(MDDstat$MDD.status)
age <- bio_cov$PPAGE_at_interview
CD4T = t(xcell[5,-1])
Bcell = t(xcell[3,-1])
CD8T = t(xcell[10,-1])
Mono = t(xcell[39,-1])
NK = t(xcell[41,-1]) # for later use
coldata922 <- data.frame(age,CD4T,Bcell,CD8T, MDD.status)
colnames(coldata922) <- c('age', 'CD4T', 'Bcell','CD8T','MDD')

x = DESeqDataSetFromMatrix(countData = t(MDD_inte), 
                           colData = coldata922, design = ~ MDD+age+CD4T+CD8T+Bcell)
MDD922 = DESeq(x, parallel = T)
save(MDD922, file ='/export/home/chche/WTC_Project/Data/MDD922_withoutMono.Rdata')
res922MDD = results(MDD922, contrast = c("MDD", "1", "2"))
save(res922MDD, file = '/export/home/chche/WTC_Project/Data/res922M.Rdata')
