library(data.table)
library(survival)
set.seed(031995)
df <- fread("/Users/andrewmarderstein/Documents/Research/druggene/output/pgs_drug_plusSNP.txt",data.table = F,stringsAsFactors = F)

df.dis0.drug0 <- subset(df,disease==0 & S01BA==0)
df.dis1.drug0 <- subset(df,disease==1 & S01BA==0)
df.dis0.drug1 <- subset(df,disease==0 & S01BA==1)
df.dis1.drug1 <- subset(df,disease==1 & S01BA==1)

res.chr16.lst <- list()
res.chr19.lst <- list()
res.pgs.lst <- list()
for (i in 1:100) {
  print(i)
  m=Inf
  n=Inf
  iter=0
  # while (m>36 & n>1) {
  while (m>36) {
    iter=iter+1;if (iter %% 50 == 0) {print(iter)}
    df1 <- 
      rbind(
        df.dis0.drug0[sample(1:nrow(df.dis0.drug0),1093,replace = F),],
        df.dis1.drug0[sample(1:nrow(df.dis1.drug0),1380,replace = F),],
        df.dis0.drug1[sample(1:nrow(df.dis0.drug1),33,replace = F),],
        df.dis1.drug1[sample(1:nrow(df.dis1.drug1),42,replace = F),]
      )
    m=nrow(subset(df1,S01BA>0 & rs4784227_chr16>0))
    # n=nrow(subset(df1,S01BA>0 & rs62119267_chr19>0))
    table(df1[,c('S01BA','rs4784227_chr16')])
    table(df1[,'S01BA']*df1[,'rs4784227_chr16'],df1[,'disease'])
    
  }
  iter=0; m=Inf
  while (m>58) {
    iter=iter+1;if (iter %% 50 == 0) {print(iter)}
    df2 <- 
      rbind(
        df.dis0.drug0[sample(1:nrow(df.dis0.drug0),3639,replace = F),],
        df.dis1.drug0[sample(1:nrow(df.dis1.drug0),3126,replace = F),],
        df.dis0.drug1[sample(1:nrow(df.dis0.drug1),63,replace = F),],
        df.dis1.drug1[sample(1:nrow(df.dis1.drug1),61,replace = F),]
      )
    m=nrow(subset(df2,S01BA>0 & rs4784227_chr16>0))
  }
  
  # mod1 <- coxph(Surv(days, disease) ~ rs62119267_chr19*S01BA+bmi+age+menopause+
  #                 number_live_birth+one_birth+
  #                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df1)
  # mod2 <- coxph(Surv(days, disease) ~ rs62119267_chr19*S01BA+bmi+age+menopause+
  #                 number_live_birth+one_birth+
  #                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df2)
  # res1 <- summary(mod1)$coef["rs62119267_chr19:S01BA",]
  # res2 <- summary(mod2)$coef["rs62119267_chr19:S01BA",]
  # names(res1) <- paste0(names(res1),".1")
  # names(res2) <- paste0(names(res2),".2")
  # res.mg <- cbind(data.frame(t(res1)),data.frame(t(res2)))
  # res.chr19.lst[[i]] <- res.mg

  mod1 <- coxph(Surv(days, disease) ~ rs4784227_chr16*S01BA+bmi+age+menopause+
                  number_live_birth+one_birth+
                  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df1)
  mod2 <- coxph(Surv(days, disease) ~ rs4784227_chr16*S01BA+bmi+age+menopause+
                  number_live_birth+one_birth+
                  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df2)
  res1 <- summary(mod1)$coef["rs4784227_chr16:S01BA",]
  res2 <- summary(mod2)$coef["rs4784227_chr16:S01BA",]
  names(res1) <- paste0(names(res1),".1")
  names(res2) <- paste0(names(res2),".2")
  res.mg <- cbind(data.frame(t(res1)),data.frame(t(res2)))
  res.chr16.lst[[i]] <- res.mg
  
  mod1 <- coxph(Surv(days, disease) ~ pgs*S01BA+bmi+age+menopause+
                  number_live_birth+one_birth+
                  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df1)
  mod2 <- coxph(Surv(days, disease) ~ pgs*S01BA+bmi+age+menopause+
                  number_live_birth+one_birth+
                  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df2)
  res1 <- summary(mod1)$coef["pgs:S01BA",]
  res2 <- summary(mod2)$coef["pgs:S01BA",]
  names(res1) <- paste0(names(res1),".1")
  names(res2) <- paste0(names(res2),".2")
  res.mg <- cbind(data.frame(t(res1)),data.frame(t(res2)))
  res.pgs.lst[[i]] <- res.mg
  
}

res.chr16 <- as.data.frame(do.call(rbind,res.chr16.lst))
res.pgs <- as.data.frame(do.call(rbind,res.pgs.lst))
# res.chr19 <- as.data.frame(do.call(rbind,res.chr19.lst))
# mean(res.chr19[,"Pr...z...1"] < 0.05)
# mean(res.chr19[,"Pr...z...2"] < 0.05)
mean(res.chr16[,"Pr...z...1"] < 0.05)
mean(res.chr16[,"Pr...z...2"] < 0.05)
mean(res.pgs[,"Pr...z...1"] < 0.05)
mean(res.pgs[,"Pr...z...2"] < 0.05)

summary(coxph(Surv(days, disease) ~ rs4784227_chr16+bmi+age+menopause+
                number_live_birth+one_birth+
                PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df2))

summary(coxph(Surv(days, disease) ~ rs62119267_chr19*S01BA+bmi+age+menopause+
                number_live_birth+one_birth+
                PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df))$coef["rs62119267_chr19:S01BA",]

summary(coxph(Surv(days, disease) ~ rs4784227_chr16*S01BA+bmi+age+menopause+
                number_live_birth+one_birth+
                PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df))$coef["rs4784227_chr16:S01BA",]

