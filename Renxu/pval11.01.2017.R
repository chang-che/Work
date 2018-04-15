

library(RobustRankAggreg)
load('/export/home/xurren/WTCProject/Data/methy_res.RData')
head(methy_res)
methy_res = methy_res[which(methy_res$RefSeq!=""), ]

genes = unique(methy_res$RefSeq)  ## genes store the gene ref sequence
numOfPro = rep(0, length(genes))  ## numOfPro is a vector contains the number of probes for all genes
minPVal = rep(0, length(genes))
corrPVal = rep(0, length(genes))
rho = rep(0, length(genes))
for(i in 1:length(genes)){
  ind = which(methy_res$RefSeq == genes[i])
  numOfPro[i] = length(ind)
  minPVal[i] = min(methy_res$pval[ind])
  corrPVal[i] = min(betaScores(methy_res$pval[ind]))
  rho[i] = rhoScores(methy_res$pval[ind])
}

save(numOfPro, minPVal, corrPVal, rho, file = "/export/home/xurren/WTCProject/Data/robust.RData")

