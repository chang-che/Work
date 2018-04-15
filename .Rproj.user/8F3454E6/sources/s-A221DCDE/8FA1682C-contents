library(edgeR)
?rpkm
MDD <- read.table(file = '/export/home/pfkuan/WTCproject/Epigenetics/Data/levinsonNIMH_MDD_EQTL_Covariates/data_used_for_eqtl_study/raw_counts.txt', header = T, sep = '\t')
l.sort = sort(apply(MDD, 2, function(x) sum(as.numeric(is.na(x)))), decreasing = T)
MDD1 <- MDD[,-ncol(MDD)]
MDD_numeric = apply(MDD1[,2:ncol(MDD1)], 2, as.numeric)
MDD_numeric <- cbind(MDD1[,1], MDD_numeric)

xcell = read.table('/Users/changche/Downloads/xCell_xcell_xCell_1852040918_RAW.txt', header = T, sep = '\t')
