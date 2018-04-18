#check the correlation of covariates from the paper
tech_cov <- read.table('/export/home/pfkuan/WTCproject/Epigenetics/Data/levinsonNIMH_MDD_EQTL_Covariates/covariates/Technical_factors.txt', header = T, sep ='\t')
bio_cov <- read.table('/export/home/pfkuan/WTCproject/Epigenetics/Data/levinsonNIMH_MDD_EQTL_Covariates/covariates/Biological_and_hidden_factors.txt', header = T, sep ='\t')
all_cov = cbind(MDD,tech_cov[,-1],bio_cov[,-1])
all_cov = apply(all_cov, 2, as.numeric)
x = cor(all_cov)
xdf = data.frame(as.table(x))
# output the number of pairs of covariates according to three thresholds of correlation
# Notice that the output number is 2 times pairs of covariates
dim(subset(xdf, abs(Freq) >=0.9&Freq!=1))
dim(subset(xdf, abs(Freq) >=0.8&Freq!=1))
dim(subset(xdf, abs(Freq) >=0.7&Freq!=1))

###the result of WTC project from Xu
load(file = '/export/home/xurren/WTCProject/Code/DESeq_additional/MDD.RData')
res330_df = as.data.frame(res330MDD)
sum(is.na(subset(res330_df, padj<0.05)$padj))#check whether na has been removed, and YES
gene330_0.05 =rownames(subset(res330_df, padj<0.05))

### the result of paper
load(file = '/export/home/chche/WTC_Project/Data/res922M.Rdata')
res922_df = as.data.frame(res922MDD)

gene922_0.05 =rownames(subset(res922_df, padj<0.05))

gene_intersect <- intersect(gene330_0.05, gene922_0.05)
length(gene_intersect)
gene_union <- union(gene330_0.05, gene922_0.05)

#only choose the genes that has non-NA log2FoldChange
gene_union_corrected <- gene_union[!is.na(res330_df[gene_union,]$log2FoldChange)&!is.na(res922_df[gene_union,]$log2FoldChange)]
cor(res330_df[gene_union_corrected,]$log2FoldChange, res922_df[gene_union_corrected,]$log2FoldChange)
