library(DESeq2)
library(BiocParallel)
register(MulticoreParam(multicoreWorkers()))#register the number of cores in that computer
load(file = '/export/home/chche/WTC_Project/Data/MDD_int.Rdata') # MDD_integer
tech_cov <- read.table('/export/home/pfkuan/WTCproject/Epigenetics/Data/levinsonNIMH_MDD_EQTL_Covariates/covariates/Technical_factors.txt', header = T, sep ='\t')
bio_cov <- read.table('/export/home/pfkuan/WTCproject/Epigenetics/Data/levinsonNIMH_MDD_EQTL_Covariates/covariates/Biological_and_hidden_factors.txt', header = T, sep ='\t')
MDDstat <- read.table('/export/home/pfkuan/WTCproject/Epigenetics/Data/levinsonNIMH_MDD_EQTL_Covariates/data_used_for_MDD_analysis/Dx_Case_status.txt', header = T, sep ='\t')
MDD <- as.factor(MDDstat$MDD.status)
all_cov = cbind(MDD,tech_cov[,-1],bio_cov[,-1])

x <- paste(colnames(all_cov) ,collapse = '+')
fomu = as.formula(paste('~',x))
x = DESeqDataSetFromMatrix(countData = t(MDD_inte), 
                           colData = all_cov, design = fomu)
MDD922 = DESeq(x, parallel = T)
save(MDD922, file ='/export/home/chche/WTC_Project/Data/MDD922_all_cov.Rdata')
res922MDD = results(MDD922, contrast = c("MDD", "1", "2"))
save(res922MDD, file = '/export/home/chche/WTC_Project/Data/res922M_all_cov.Rdata')
