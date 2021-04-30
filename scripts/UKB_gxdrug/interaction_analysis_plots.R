library(data.table)
library(ggplot2)
df <- fread("/Users/andrewmarderstein/Documents/Research/druggene/output/pgs_drug.txt",data.table = F,stringsAsFactors = F)
# df$pgs <- -1*df$pgs
# df$P_0.00001.R2_0.1.KB_250 <- -1*df$P_0.00001.R2_0.1.KB_250
# df$PGS000007 <- -1*df$PGS000007
cor.test(df$P_0.00001.R2_0.1.KB_250,df$PGS000007)
g <- ggplot(df,aes(x=P_0.00001.R2_0.1.KB_250,PGS000007)) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  labs(x='Michailidou et al. PGS',y='Mavaddat et al. PGS') + 
  geom_hex() + 
  # geom_point() + 
  geom_smooth(method='lm',col='darkred',se=F) +
  labs(fill='N')
png('~/Documents/Research/druggene/output/pgs_int/figure4a.png',height=2*2080/3,width=2*2630/3,res=400)
print(g)
dev.off()

mod<-glm(disease ~ (df[,'pgs'])+bmi+age+menopause+number_live_birth+one_birth+
           PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df,family = binomial(link='logit'))
combined_beta <- exp(summary(mod)$coef[2,1])
combined_low <- exp(summary(mod)$coef[2,1]-1.96*summary(mod)$coef[2,2])
combined_high <- exp(summary(mod)$coef[2,1]+1.96*summary(mod)$coef[2,2])
mod <- glm(disease ~ (df[,'P_0.00001.R2_0.1.KB_250'])+bmi+age+menopause+number_live_birth+one_birth+
             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df,family = binomial(link='logit'))
internal_beta <- exp(summary(mod)$coef[2,1])
internal_low <- exp(summary(mod)$coef[2,1]-1.96*summary(mod)$coef[2,2])
internal_high <- exp(summary(mod)$coef[2,1]+1.96*summary(mod)$coef[2,2])
mod<-glm(disease ~ (df[,'PGS000007'])+bmi+age+menopause+number_live_birth+one_birth+
           PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df,family = binomial(link='logit'))
external_beta <- exp(summary(mod)$coef[2,1])
external_low <- exp(summary(mod)$coef[2,1]-1.96*summary(mod)$coef[2,2])
external_high <- exp(summary(mod)$coef[2,1]+1.96*summary(mod)$coef[2,2])

library(survival)
mod<-coxph(Surv(days,disease) ~ (df[,'pgs'])+bmi+age+menopause+number_live_birth+one_birth+
           PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df)
combined_beta <- (summary(mod)$conf.int[1,1])
combined_low <- (summary(mod)$conf.int[1,3])
combined_high <- (summary(mod)$conf.int[1,4])
mod<-coxph(Surv(days,disease) ~ (df[,'P_0.00001.R2_0.1.KB_250'])+bmi+age+menopause+number_live_birth+one_birth+
             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df)
internal_beta <- (summary(mod)$conf.int[1,1])
internal_low <- (summary(mod)$conf.int[1,3])
internal_high <- (summary(mod)$conf.int[1,4])
mod<-coxph(Surv(days,disease) ~ (df[,'PGS000007'])+bmi+age+menopause+number_live_birth+one_birth+
             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df)
external_beta <- (summary(mod)$conf.int[1,1])
external_low <- (summary(mod)$conf.int[1,3])
external_high <- (summary(mod)$conf.int[1,4])
res <- data.frame(PGS=c('Combined','Michailidou et al.','Mavaddat et al.'))
res$OR <- c(combined_beta,internal_beta,external_beta)
res$low <- c(combined_low,internal_low,external_low)
res$high <- c(combined_high,internal_high,external_high)
g1 <- ggplot(data=res, aes(x=as.factor(PGS), y=OR, ymin=low, ymax=high,col=PGS,fill=PGS)) +
  geom_linerange(size=5,position=position_dodge(width = 0.5)) +
  geom_point(size=3, shape=21, colour="white", fill='black',stroke = 0.5,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2,col='red') +
  coord_flip() +
  labs(x="PGS type",y="Estimated risk increase per SD increase in PGS (HR)") +
  theme_bw()  +
  theme(axis.text = element_text(family = 'Helvetica'),legend.position = 'none') + 
  scale_colour_manual(values=c('darkred','darkblue','darkgreen'))
library("rcompanion")
combined <- nagelkerke(glm(disease ~ scale(df[,'pgs'])+bmi+age+
                             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df,family = binomial(link='logit')),
                       glm(disease ~ bmi+age+
                             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df,family = binomial(link='logit'))
)$Pseudo.R.squared.for.model.vs.null[3,]
internal <- nagelkerke(glm(disease ~ scale(df[,'P_0.00001.R2_0.1.KB_250'])+bmi+age+
                             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df,family = binomial(link='logit')),
                       glm(disease ~ bmi+age+
                             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df,family = binomial(link='logit'))
)$Pseudo.R.squared.for.model.vs.null[3,]
external <- nagelkerke(glm(disease ~ scale(df[,'PGS000007'])+bmi+age+
                             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df,family = binomial(link='logit')),
                       glm(disease ~ bmi+age+
                             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = df,family = binomial(link='logit'))
)$Pseudo.R.squared.for.model.vs.null[3,]
tmp <- data.frame(type=c('Combined','Michailidou et al.','Mavaddat et al.'),R2=c(combined,internal,external))
g2 <- ggplot(tmp,aes(x=type,y=R2,fill=type)) + geom_bar(stat = "identity") + scale_fill_manual(values=c('darkred','darkblue','darkgreen')) + theme_bw() + theme(panel.grid=element_blank()) +
  labs(x='PGS type',y='Nagelkerke R-squared') + guides(fill=F) + coord_flip()
library(cowplot); 
png("~/Documents/Research/druggene/output/pgs_int/figure4b-c.png",height = (1910/3)*2,width=(4010/3)*2,res = 2*200)
plot_grid(g1,g2,ncol=1)
dev.off()

# Calculate percentiles
rng <- seq(0.01,1, by = 0.01)
for(i in rng){
  if(i == min(rng)){
    cut.high <- quantile(df$pgs, probs = i)
    dat <- apply(df[which(df$pgs <= cut.high),],2,mean)
  }
  else{
    cut.low <- quantile(df$pgs,probs = i-min(rng))
    cut.high <- quantile(df$pgs, probs = i)
    ind.sub <- which(df$pgs <= cut.high & df$pgs >= cut.low)
    dat <- apply(df[ind.sub,],2,mean)
  }
  dat <- as.data.frame(matrix(dat,nrow = 1)); colnames(dat) <- colnames(df)
  dat$Percentile <- i
  if (i == min(rng)) {
    df.percentile <- dat
  } else {
    df.percentile <- rbind(df.percentile,dat)
  }
}
ggplot(df.percentile,aes(x=Percentile,y=100*disease,col=Percentile)) + geom_point() + scale_colour_gradient(low='steelblue',high='red') +
  theme_bw() + theme(panel.grid = element_blank()) + labs(x='Percentile of polygenic score',y='Breast cancer incidence (%)') + guides(col=F)

# Calculate percentiles
seper <- 0.1
rng <- seq(seper,1, by = seper)
pgs_name <- 'pgs'
pgs_list <- c('pgs','P_0.00001.R2_0.1.KB_250','PGS000007')

val <- list(); for (j in 1:3) {
  pgs_name <- pgs_list[j]
  cut.low <- quantile(df[,pgs_name],probs = 0.5 - seper/2)
  cut.high <- quantile(df[,pgs_name], probs = 0.5 + seper/2)
  ind.sub <- which(df[,pgs_name] <= cut.high & df[,pgs_name] >= cut.low)
  val[[j]] <- mean(df[ind.sub,'disease'])
}

dat.save.all <- list()
for (j in 1:3) {
  pgs_name <- pgs_list[j]
  dat.save <- c()
  for(i in rng){
    if(i == min(rng)) {
      cut.high <- quantile(df[,pgs_name], probs = i)
      dat <- mean(df[which(df[,pgs_name] <= cut.high),'disease'])
    }
    else {
      cut.low <- quantile(df[,pgs_name],probs = i-min(rng))
      cut.high <- quantile(df[,pgs_name], probs = i)
      ind.sub <- which(df[,pgs_name] <= cut.high & df[,pgs_name] >= cut.low)
      dat <- mean(df[ind.sub,'disease'])
    }
    dat.save <- c(dat.save,dat)
  }
  tmp <- data.frame(type=pgs_name,Percentile=rng,disease_nonadj=dat.save)
  tmp$disease <- tmp$disease_nonadj/val[[j]]
  dat.save.all[[j]] <- tmp
}
dat.save <- do.call(rbind,dat.save.all)
dat.save <- subset(dat.save,Percentile <= 0.4 | Percentile >= 0.7)
ggplot(dat.save,aes(x=Percentile,y=disease,col=type)) + geom_point() + geom_hline(yintercept=1,col='red',lty='dashed') +
  theme_bw() + theme(panel.grid = element_blank()) + labs(x='Percentile of polygenic score',y='Relative change in incidence')






i=1;seper=0.2
# pgs_name = 'pgs'
# pgs_name = 'P_0.00001.R2_0.1.KB_250'
pgs_name = 'PGS000007'
cut.low <- quantile(df[,pgs_name],probs = i-seper)
cut.high <- quantile(df[,pgs_name], probs = i)
ind.sub <- which(df[,pgs_name] <= cut.high & df[,pgs_name] >= cut.low)
mean(df[ind.sub,'disease'])

df$strata <- NA
ind.high <- df$pgs > quantile(df$pgs,probs=0.8)
ind.mod <- df$pgs <= quantile(df$pgs,probs=0.8) & df$pgs > quantile(df$pgs,probs=0.2)
ind.low <- df$pgs <= quantile(df$pgs,probs=0.2)
df$strata[ind.high] <- 'Top 20%'
df$strata[ind.mod] <- 'Mid 60%'
df$strata[ind.low] <- 'Bot 20%'
aggregate(df$disease,by=list(strata=df$strata),mean)
library(survival)
library(survminer)
library(ggfortify)
fit <- survfit(Surv(days, disease) ~ strata, data = df)
# fit <- survfit(Surv(days, disease) ~ strata, data = df[1:1000,])
# g1 <- autoplot(fit,conf.int = FALSE, censor = FALSE) + theme_bw() + theme(panel.grid = element_blank()) +
#   labs(x = "Days", y = "Proportion without disease")
ggsurvplot(fit,
           conf.int = FALSE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw()+theme(panel.grid = element_blank()), # Change ggplot2 theme
           # palette = c("#E7B800", "#2E9FDF"),
           fun = "event",xlab='Days')


tmp1 <- subset(tmp,S01BA==1)
fit <- survfit(Surv(days, disease) ~ pgs_strata, data = tmp1)
g1 <- autoplot(fit,conf.int = FALSE, censor = FALSE) + theme_bw() + theme(panel.grid = element_blank()) +
  labs(x = "Days", y = "Proportion without disease") +
  guides(fill=FALSE) + labs(colour = "User") + lims(y=c(0,1))
tmp0 <- subset(tmp,S01BA==0)
fit <- survfit(Surv(days, P) ~ drug, data = tmp0)
g2 <- autoplot(fit,conf.int = FALSE, censor = FALSE) + theme_bw() + theme(panel.grid = element_blank()) +
  labs(x = "Days", y = "Proportion without disease") +
  guides(fill=FALSE) + labs(colour = "User") +
  lims(y=c(0,1))
plot_grid(g2,g1,ncol = 2)


