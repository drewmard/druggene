library(data.table)
library(survival)
library(BEDMatrix)

dataf.mg <- fread('/athena/elementolab/scratch/anm2868/druggene/output/pgs_drug.txt',data.table = F,stringsAsFactors = F)

direc.vec <- c("/athena/elementolab/scratch/anm2868/druggene/output/ss/PGS000007",
               "/athena/elementolab/scratch/anm2868/druggene/output/ss/michailidou/clumped/P_0.00001.R2_0.1.KB_250")

iter=0
res <- list()
drug='S01BA'
# mod <- coxph(Surv(days, disease) ~ pgs*dataf.mg[,drug]+bmi+age+menopause+
#                number_live_birth+one_birth+
#                PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = dataf.mg)
for (ind.dir in 1:1) {
# for (ind.dir in 1:length(direc.vec)) {
  for (CHR in 1:22) {
    f.geno <- paste0(direc.vec[ind.dir],'/geno.',CHR)
    # f.geno <- paste0('/athena/elementolab/scratch/anm2868/DrugPGS/DrugPGS_Interactions/output/PGS/PGS_Catalog/PGS000007/geno.',CHR)
    geno <- BEDMatrix(f.geno)
    geno_names <- unlist(lapply(strsplit(rownames(geno),'_'),function(x) x[[2]]))
    ind <- match(dataf.mg$eid,geno_names)
    P=ncol(geno)
    print(paste0('Pt ',ind.dir,': Running chr',CHR,'...'))
    for (i in 1:P) {
      if (i %% 10 == 0) {print(paste0(i,'/',P))}
      iter=iter+1
      dataf.mg$SNP <- geno[ind,i]
      drug <- 'S01BA'
      mod <- coxph(Surv(days, disease) ~ SNP*dataf.mg[,drug]+bmi+age+menopause+
                     number_live_birth+one_birth+
                     PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = dataf.mg)
      res[[iter]] <- data.frame(
        ind=ind.dir,
        CHR=CHR,
        SNP=colnames(geno)[i],
        beta_snp=summary(mod)$coef['SNP','exp(coef)'],
        p_snp=summary(mod)$coef['SNP','Pr(>|z|)'],
        beta_int=summary(mod)$coef['SNP:dataf.mg[, drug]','exp(coef)'],
        p_int=summary(mod)$coef['SNP:dataf.mg[, drug]','Pr(>|z|)']
      )
    }
  }
}
res.save <- do.call(rbind,res)
res.save$ind[duplicated(res.save$SNP,fromLast = T)] <- 3
res.save <- res.save[!duplicated(res.save$SNP),]
res.save[order(res.save$p_int,decreasing = F)[1:10],]
f.out <- paste0('/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.both_pgs.txt')
fwrite(res.save,f.out,quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)


############