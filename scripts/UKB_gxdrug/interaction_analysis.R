library(data.table)
library(reshape2)
library(survival)

origdir <- '/athena/elementolab/scratch/anm2868/druggene'
iter <- 0
polyscore <- list()
for (PTHRES in c("0.05","0.005","0.0005","0.00001","5e-8")) {
  for (R2THRES in c("0.1")) {
    for (KBTHRES in c("250")) {
      iter <- iter + 1
      type <- paste0('P_',PTHRES,'.R2_',R2THRES,'.KB_',KBTHRES)
      print(paste0(iter,': ',type))
      dir <- paste0(origdir,'/output/ss/michailidou/clumped/',type)
      polyscore2 <- fread(paste0(dir,'/score_chrALL.profile'),data.table=F,stringsAsFactors = F)
      colnames(polyscore2)[2:23] <- paste0('score_',1:22)
      polyscore2[,2:ncol(polyscore2)] <- polyscore2[,2:ncol(polyscore2)]
      polyscore2$score_all <- apply(polyscore2[,-1],1,sum)
      polyscore2$score_all <- polyscore2$score_all#*-1
      polyscore2$score_all <- scale(polyscore2$score_all)[,1]
      colnames(polyscore2)[1] <- 'eid'; polyscore2[,type] <- polyscore2$score_all # this is necessary to work w/ previous script
      polyscore[[iter]] <- polyscore2[,c('eid',type)]
    }
  }
}
# for (PTHRES in c("0.0005","0.00001","5e-8")) {
#   for (R2THRES in c("0.1")) {
#     for (KBTHRES in c("250","5000")) {
#       iter <- iter + 1
#       type <- paste0('P_',PTHRES,'.R2_',R2THRES,'.KB_',KBTHRES)
#       print(paste0(iter,': ',type))
#       dir <- paste0(origdir,'/output/ss/michailidou/clumped/',type)
#       polyscore2 <- fread(paste0(dir,'/score_chrALL.profile'),data.table=F,stringsAsFactors = F)
#       colnames(polyscore2)[2:23] <- paste0('score_',1:22)
#       polyscore2[,2:ncol(polyscore2)] <- polyscore2[,2:ncol(polyscore2)]
#       polyscore2$score_all <- apply(polyscore2[,-1],1,sum)
#       polyscore2$score_all <- polyscore2$score_all#*-1
#       polyscore2$score_all <- scale(polyscore2$score_all)[,1]
#       colnames(polyscore2)[1] <- 'eid'; polyscore2[,type] <- polyscore2$score_all # this is necessary to work w/ previous script
#       polyscore[[iter]] <- polyscore2[,c('eid',type)]
#     }
#   }
# }


polyscore <- do.call(cbind,polyscore)
polyscore <- polyscore[,!duplicated(colnames(polyscore))]

iter <- 0
polyscore2.save <- list()
for (type in c('PGS000004','PGS000007')) {
  for (FLAG in c('.no_remove','')) {
  # for (FLAG in c('')) {
    iter <- iter + 1
    print(paste0(iter,': ',type,FLAG))
    f <- paste0(origdir,'/output/ss/',type,'/score_chrALL',FLAG,'.profile')
    polyscore2 <- fread(f,data.table=F,stringsAsFactors = F)
    colnames(polyscore2)[2:ncol(polyscore2)] <- paste0('score_',1:c(ncol(polyscore2)-1))
    polyscore2$score_all <- apply(polyscore2[,-1],1,sum)
    polyscore2$score_all <- polyscore2$score_all*-1
    polyscore2$score_all <- scale(polyscore2$score_all)[,1]
    colnames(polyscore2)[1] <- 'eid'; polyscore2[,paste0(type,FLAG)] <- polyscore2$score_all # this is necessary to work w/ previous script
    polyscore2.save[[iter]] <- polyscore2[,c('eid',paste0(type,FLAG))]
  }
}
polyscore2 <- do.call(cbind,polyscore2.save)
polyscore2 <- polyscore2[,!duplicated(colnames(polyscore2))]

polyscore.full <- merge(polyscore,polyscore2,by='eid')
cor.mat <- cor(polyscore.full[,-1])
fwrite(as.data.frame(cor.mat),'/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/cor_mat.txt',quote = F,na = 'NA',sep = '\t',row.names = T,col.names = T)
cor.df <- melt(cor.mat)
fwrite(cor.df,'/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/cor_df.txt',quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)

#######
# set up disease data too:

polyscore <- polyscore.full[,c(1:2)]
colnames(polyscore)[2] <- 'breast_cancer'

# load
level=4; drug_data <- fread(paste0('/athena/elementolab/scratch/anm2868/DrugPGS/UKB_drug_data/level',level,'.txt'),data.table = F,stringsAsFactors = F)
pheno2 <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb33822.csv.gz',data.table=F,stringsAsFactors = F)
drug_names_full <- colnames(drug_data)[which(colnames(drug_data) != 'eid')]

disease_lst <- c('breast_cancer');j = 1
use_percentile <- FALSE; remainder_popu=TRUE; percentile <- .2
returned <- prepare_data_for_analysis(j)
disease_name <- returned[[1]]
drug_names <- returned[[2]]
dataf.mg <- returned[[3]]


type.lst <- c()
for (PTHRES in c("0.05","0.005","0.0005","0.00001","5e-8")) {
  for (R2THRES in c("0.1")) {
    for (KBTHRES in c("250")) {
# for (PTHRES in c("0.0005","0.00001","5e-8")) {
#   for (R2THRES in c("0.1")) {
#     for (KBTHRES in c("250","5000")) {
      type <- paste0('P_',PTHRES,'.R2_',R2THRES,'.KB_',KBTHRES)
      type.lst <- c(type.lst,type)
    }
  }
}
for (type in c('PGS000004','PGS000007')) {
  for (FLAG in c('')) {
  # for (FLAG in c('.no_remove','')) {
    type2 <- paste0(type,FLAG)
    type.lst <- c(type.lst,type2)
  }
}

dataf.mg <- merge(dataf.mg,polyscore.full,by='eid')
disease_name <- disease_lst[j]
p.vec <- c()
for (i in 1:length(type.lst)) {
  mod <- coxph(Surv(days, disease) ~ dataf.mg[,type.lst[i]]+bmi+age+menopause+
                 number_live_birth+one_birth+
                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = dataf.mg)
  p <- summary(mod)$coef[1,2]
  p.vec <- c(p.vec,p)
  
}
cbind(type.lst,p.vec)
mod<-glm(disease ~ dataf.mg[,type.lst[i]]+bmi+age+menopause+
               number_live_birth+one_birth+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = dataf.mg,family = binomial(link='logit'))
summary(mod)$coef[1,]

summary(glm(rbinom(5,1,0.5)~rnorm(5),family = binomial(link='logit')))

print('Starting analysis...')
results <- list()
for (j in 1:length(type.lst)) {
  print(paste0(j,': ',type.lst[j]))
  res <- list()
  for (i in 1:length(drug_names)) {
  # for (i in 58:58) {
    drug <- drug_names[i]
    if (i%%10==0) {print(paste0('Drug: ',drug,' (',i,'/',length(drug_names),')'))}
    mod <- coxph(Surv(days, disease) ~ dataf.mg[,type.lst[j]]*dataf.mg[,drug]+bmi+age+menopause+
                   number_live_birth+one_birth+
                   PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = dataf.mg)
    res[[i]] <- data.frame(PGS=type.lst[j],drug=drug,hr=summary(mod)$coef["dataf.mg[, type.lst[j]]:dataf.mg[, drug]",2], p=summary(mod)$coef["dataf.mg[, type.lst[j]]:dataf.mg[, drug]",5])
  }
  results[[j]] <- do.call(rbind,res)
}
results_all <- do.call(rbind,results)
results_all$N <- apply(dataf.mg[,drug_names],2,sum)
results_all[order(results_all$p)[1:20],]


# results_all <- subset(results_all,N>=10000)
results_all[[k]]$fdr <- p.adjust(results_all[[k]]$p,method = 'fdr')
results_all[[k]]$k <- k
# results_all[[k]][order(results_all[[k]]$p)[1:20],]
# f.out <- paste0('/athena/elementolab/scratch/anm2868/DrugPGS/UKB_drug_results/level',level,'_screen.percentile_',percentile,'.remainder_',remainder_popu,'.txt')
# fwrite(results_all,f.out,quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)

}
results_all_save <- do.call(rbind,results_all)
results_all_save[order(results_all_save$p)[1:20],]




#####################
dataf.mg$breast_cancer_pgs <- apply(dataf.mg[,c('P_0.0005.R2_0.1.KB_5000','PGS000007')],1,mean)
dataf.mg$breast_cancer_pgs <- apply(dataf.mg[,c('P_0.0005.R2_0.1.KB_250','PGS000007')],1,mean)
# dataf.mg$breast_cancer_pgs <- apply(dataf.mg[,c('P_0.00001.R2_0.1.KB_5000','PGS000007')],1,mean)
for (i in 1:length(drug_names)) {
  drug <- drug_names[i]
  if (i%%10==0) {print(paste0('Drug: ',drug,' (',i,'/',length(drug_names),')'))}
  mod <- coxph(Surv(days, disease) ~ breast_cancer_pgs*dataf.mg[,drug]+bmi+age+menopause+
                 number_live_birth+one_birth+
                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = dataf.mg)
  res[[i]] <- data.frame(PGS='breast_cancer_pgs',drug=drug,hr=summary(mod)$coef["breast_cancer_pgs:dataf.mg[, drug]",2], p=summary(mod)$coef["breast_cancer_pgs:dataf.mg[, drug]",5])
}
results <- do.call(rbind,res)
results$fdr <- p.adjust(results$p,method = 'fdr')
results[order(results$p)[1:20],]

