library(data.table)
library(ggplot2)
f <- "/Users/andrewmarderstein/Documents/Research/druggene/output/snp_int/SNP_by_drug.BC.txt"
df <- fread(f,data.table = F,stringsAsFactors = F)
gwas <- fread("/Users/andrewmarderstein/Documents/Research/druggene/output/snp_int/clean_PGS000007.txt.gz",data.table = F,stringsAsFactors = F)
df$SNP.adj <- substring(df$SNP,1,as.numeric(gregexpr('_',df$SNP))-1)
df <- merge(df,gwas[,c('RSID','BETA','BP')],by.x='SNP.adj',by.y='RSID')
df$FDR <- p.adjust(df$p_int,method = 'fdr')
df.sub <- subset(df,FDR<0.5)# & p_snp <0.05)
df.sub <- subset(df,p_int<0.05)# & p_snp <0.05)
df.sub <- df.sub[order(df.sub$CHR),]
binom.test(mean(sign(df.sub$BETA)==sign(log(df.sub$beta_int)))*nrow(df.sub),nrow(df.sub),0.5)
x2 <- binom.test(mean(sign(df.sub$BETA)==-1 & sign(log(df.sub$beta_int))==1)*nrow(df.sub),nrow(df.sub),0.25)
x3 <- binom.test(mean(sign(df.sub$BETA)==1 & sign(log(df.sub$beta_int))==1)*nrow(df.sub),nrow(df.sub),0.25)
x4 <- binom.test(mean(sign(df.sub$BETA)==-1 & sign(log(df.sub$beta_int))==-1)*nrow(df.sub),nrow(df.sub),0.25)
tmp <- list(); ind = 0
for (i in c(-1,1)) {
  for (j in c(-1,1)) {
    ind <- ind + 1
    # x1 <- binom.test(mean(sign(log(df.sub$beta_snp))==i & sign(log(df.sub$beta_int))==j)*nrow(df.sub),nrow(df.sub),0.25)
    x1 <- binom.test(mean(sign(df.sub$BETA)==i & sign(log(df.sub$beta_int))==j)*nrow(df.sub),nrow(df.sub),0.25)
    tmp[[ind]] <- data.frame(marginal_sign=i,interaction_sign=j,name=paste0(ifelse(i==1,'+','-'),'/',ifelse(j==1,'+','-')),coord=as.numeric(i==j),low=x1$conf.int[1],est=as.numeric(x1$estimate),high=x1$conf.int[2],p=x1$p.value)
  }
}
tmp <- do.call(rbind,tmp)
ggplot(tmp,aes(x=name,y=est,ymin=low,ymax=high,col=name)) + geom_pointrange() + 
  geom_hline(yintercept=0.25,lty='dashed',col='red') + scale_color_brewer(palette = 'Reds') + 
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(x='Sign: Marginal effect/interaction effect',y='Proportion of nominal interactions (P<0.05)') + guides(col=F)

df.sub.sub <- subset(df.sub,BETA > 0)
ggplot(subset(df.sub.sub,BETA!=0),aes(x=log(beta_int),fill=as.factor(sign(log(beta_int))))) + geom_histogram(col='black') + # geom_density(alpha=0.2)
  guides(fill=F) + scale_fill_discrete() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Observed interaction effect',y='Count')

x1 <- binom.test(mean(sign(log(df.sub.sub$beta_int))==j)*nrow(df.sub.sub),nrow(df.sub.sub),0.5)
x1 <- binom.test(mean(sign(log(df.sub.sub$beta_int))==-1)*nrow(df.sub.sub),nrow(df.sub.sub),0.5)

tmp[[ind]] <- data.frame(marginal_sign=i,interaction_sign=j,name=paste0(ifelse(i==1,'+','-'),'/',ifelse(j==1,'+','-')),coord=as.numeric(i==j),low=x1$conf.int[1],est=as.numeric(x1$estimate),high=x1$conf.int[2],p=x1$p.value)

mean(log(df.sub.sub $beta_int) > 0)
mean(log(df.sub.sub $beta_int) < 0)

library(qqman)
manhattan(df,p='p_int',suggestiveline = -log10(max(subset(df,FDR<0.1)$p_int)))
binom.test(3,28,0.25)
df.sub
