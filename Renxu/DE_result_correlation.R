

load(file = "/export/home/xurren/WTCProject/Results/DE_WB_result.RData")
load(file = "/export/home/xurren/WTCProject/Results/DE_4cell_result.RData")
load(file = "/export/home/xurren/WTCProject/Results/DE_PBMC_result.RData")

CD4T_top = rownames(res_4cell_CD4T[order(res_4cell_CD4T$padj),])
CD8T_top = rownames(res_4cell_CD8T[order(res_4cell_CD8T$padj),])
Mono_top = rownames(res_4cell_Mono[order(res_4cell_Mono$padj),])
Bcell_top = rownames(res_4cell_BCell[order(res_4cell_BCell$padj),])
PBMC_top = rownames(res_PBMC[order(res_PBMC$padj),])

K = c(50, 100, 500, 1000, 2000, 5000, 10000, length(PBMC_top))

cat("Pair-wise correlation among log2 fold change")
for (i in 1:length(K)){
  unionGenes = unique(c(CD4T_top[1:K[i]], CD8T_top[1:K[i]], Mono_top[1:K[i]], Bcell_top[1:K[i]], PBMC_top[1:K[i]]))
  lfc_mat =cbind(res_4cell_CD4T[unionGenes, ]$log2FoldChange, 
                 res_4cell_CD8T[unionGenes, ]$log2FoldChange,
                 res_4cell_Mono[unionGenes, ]$log2FoldChange,
                 res_4cell_BCell[unionGenes, ]$log2FoldChange,
                 res_PBMC[unionGenes, ]$log2FoldChange)
  colnames(lfc_mat) = c("CD4T", "CD8T", "Mono", "BCell", "PBMC")
  cat(paste("Top", K[i], "genes: "))
  cat("\n")
  print(cor(lfc_mat, use = "pairwise.complete.obs"))
  cat("\n")
  
}

print(" ", quote = FALSE)
print(" ", quote = FALSE)
print(" ", quote = FALSE)
print(" ", quote = FALSE)
cat("Pair-wise correlation among estimated -log(p-value)")
for (i in 1:length(K)){
  unionGenes = unique(c(CD4T_top[1:K[i]], CD8T_top[1:K[i]], Mono_top[1:K[i]], Bcell_top[1:K[i]], PBMC_top[1:K[i]]))
  pval_mat =cbind(res_4cell_CD4T[unionGenes, ]$pvalue, 
                  res_4cell_CD8T[unionGenes, ]$pvalue,
                  res_4cell_Mono[unionGenes, ]$pvalue,
                  res_4cell_BCell[unionGenes, ]$pvalue,
                  res_PBMC[unionGenes, ]$pvalue)
  negLogPval = -log(pval_mat)
  colnames(negLogPval) = c("CD4T", "CD8T", "Mono", "BCell", "PBMC")
  cat(paste("Top", K[i], "genes: "))
  cat("\n")
  print(cor(negLogPval, use = "pairwise.complete.obs"))
  cat("\n")
  
}









