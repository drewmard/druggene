library(data.table)
library(survival)

dataf.mg <- fread('/athena/elementolab/scratch/anm2868/druggene/output/pgs_drug_names.txt',data.table = F,stringsAsFactors = F)
drug_names <- colnames(dataf.mg)[27:ncol(dataf.mg)]
drug_names <- drug_names[(apply(dataf.mg[,drug_names],2,sum) > 1000)]

# marginal drug association:
res <- list()
for (i in 1:length(drug_names)) {
  drug <- drug_names[i]
  if (i%%10==0) {print(paste0('Drug: ',drug,' (',i,'/',length(drug_names),')'))}
  mod <- coxph(Surv(days, disease) ~ dataf.mg[,drug]+bmi+age+
                 menopause + number_live_birth +
                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = dataf.mg)
  res[[i]] <- data.frame(PGS='pgs',drug=drug,N=sum(dataf.mg[,drug]),hr=summary(mod)$coef["dataf.mg[, drug]",2], p=summary(mod)$coef["dataf.mg[, drug]",5])
}
results <- do.call(rbind,res)
# results <- results[1:60,]
results$fdr <- p.adjust(results$p,method = 'fdr')
results[order(results$p)[1:20],]
fwrite(results[,-1],'/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/marginal_drug_assoc.txt',sep = '\t',quote = F,na = "NA",row.names = F,col.names = T)

#
res <- list()
for (i in 1:length(drug_names)) {
  drug <- drug_names[i]
  if (i%%10==0) {print(paste0('Drug: ',drug,' (',i,'/',length(drug_names),')'))}
  mod <- coxph(Surv(days, disease) ~ pgs*dataf.mg[,drug]+bmi+age+
                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = dataf.mg)
  res[[i]] <- data.frame(PGS='pgs',drug=drug,hr=summary(mod)$coef["pgs:dataf.mg[, drug]",2], p=summary(mod)$coef["pgs:dataf.mg[, drug]",5])
}
results <- do.call(rbind,res)
# results <- results[1:60,]
results$fdr <- p.adjust(results$p,method = 'fdr')
results[order(results$p)[1:20],]
fwrite(results,'/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/interaction_results.txt',sep = '\t',quote = F,na = "NA",row.names = F,col.names = T)

sum(results$p<.05); mean(results$p<.05)
sum(results$fdr<.1); mean(results$fdr<.1)

marg <- fread('/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/marginal_drug_assoc.txt',data.table = F,stringsAsFactors = F)
gxdrug <- fread('/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/interaction_results.txt',data.table = F,stringsAsFactors = F)
marg[order(marg$p)[1:20],]
results <- subset(gxdrug,drug %in% subset(marg,fdr>0.01)$drug)
sum(results$p<.05); mean(results$p<.05)
sum(results$fdr<.1); mean(results$fdr<.1)
results[order(results$p)[1:sum(results$p<.05)],]




