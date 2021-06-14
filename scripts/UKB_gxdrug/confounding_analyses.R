library(data.table)
library(survival)
library(parallel)

df.mg <- fread("/athena/elementolab/scratch/anm2868/druggene/output/other/gene_drug_propenscore.txt",data.table = F,stringsAsFactors = F)
df.mg$CRP[is.na(df.mg$CRP)] <- median(df.mg$CRP,na.rm = T)
df.mg$bmi[is.na(df.mg$bmi)] <- median(df.mg$bmi,na.rm = T)
df.mg$townsend[is.na(df.mg$townsend)] <- median(df.mg$townsend,na.rm = T)

snp_data <- fread("/athena/elementolab/scratch/anm2868/druggene/output/pgs_drug_plusSNP.txt",data.table = F,stringsAsFactors = F)
# i <- match(snp_data$eid,df.mg$eid)
i <- match(df.mg$eid,snp_data$eid)
df.mg$pgs <- snp_data[i,'pgs']
df.mg$rs4784227_chr16 <- snp_data[i,'rs4784227_chr16']
df.mg$rs62119267_chr19 <- snp_data[i,'rs62119267_chr19']


# pheno1 <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)

tmp <- data.frame(eid=pheno1$eid,smoking=as.numeric(pheno1[,'20116-0.0']>0),employ=pheno1[,'26412-0.0'])

df.mg2 <- merge(df.mg,tmp,by='eid')

simvastatin <- fread("/athena/elementolab/scratch/anm2868/DrugPGS/UKB_drug_data/simvastatin",data.table = F,stringsAsFactors = F)
df.mg2 <- merge(df.mg2,simvastatin,by='eid')

mod <- coxph(Surv(days,disease) ~ rs62119267_chr19*S01BA+
               bmi*rs62119267_chr19+age*rs62119267_chr19+
               menopause*rs62119267_chr19+number_live_birth*rs62119267_chr19+one_birth*rs62119267_chr19+
               smoking*rs62119267_chr19+employ*rs62119267_chr19+townsend*rs62119267_chr19+
               X1140861958*rs62119267_chr19+
               bmi*S01BA+age*S01BA+
               menopause*S01BA+number_live_birth*S01BA+one_birth*S01BA+
               smoking*S01BA+employ*S01BA+townsend*S01BA+
               X1140861958*S01BA+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df.mg2)
summary(mod)$coef

mod <- coxph(Surv(days,disease) ~ pgs*S01BA+
               bmi*pgs+age*pgs+
               menopause*pgs+number_live_birth*pgs+one_birth*pgs+
               smoking*pgs+employ*pgs+townsend*pgs+
               X1140861958*pgs+
               bmi*S01BA+age*S01BA+
               menopause*S01BA+number_live_birth*S01BA+one_birth*S01BA+
               smoking*S01BA+employ*S01BA+townsend*S01BA+
               X1140861958*S01BA+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df.mg2)
summary(mod)$coef

summary(glm(X1140861958~rs62119267_chr19,df.mg2,family = binomial(link="logit")))


