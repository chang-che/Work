all_cov = cbind(MDD,tech_cov[,-1],bio_cov[,-1])
all_cov_factor <- all_cov
apply(all_cov[x], 2, function(x) length(unique(x)))
x = apply(all_cov, 2, function(x) length(unique(x))<=24)
factor.choice = x
factor.choice = colnames(all_cov)[x]
factor.choice = factor.choice[-c(1,2,3,4)]


for (i in factor.choice){
  all_cov_factor[,i] = factor(all_cov_factor[,i])
}

apply(all_cov_factor, 2, function(x) length(unique(x)<=24))
all_cov_factor <-  all_cov_factor[setdiff(names(all_cov_factor), c('Tc_act','PC', 'NK_act'))]
x <- paste(colnames(all_cov_factor) ,collapse = '+')
fomu = as.formula(paste('~',x))
x = DESeqDataSetFromMatrix(countData = t(MDD_inte), 
                           colData = all_cov_factor, design = fomu)

rownames(res922MDD)[which.min(res922MDD$padj)]

