library(data.table)
res.save <- fread('/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.both_pgs.txt',data.table = F,stringsAsFactors = F)
res.save$FDR <- p.adjust(res.save$p_int,method = 'fdr')
# fwrite(res.save,"/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.both_pgs.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
sum(res.save$FDR<0.1,na.rm = T); mean(res.save$FDR<0.1,na.rm=T)
sum(res.save$p_int<0.05,na.rm = T); mean(res.save$p_int<0.05,na.rm = T)

gwas1 <- fread("/athena/elementolab/scratch/anm2868/druggene/output/ss/PGS000007/clean_PGS000007.txt.gz",data.table = F,stringsAsFactors = F)
gwas2 <- fread("/athena/elementolab/scratch/anm2868/druggene/output/ss/michailidou/clean_michailidou.txt.gz",data.table = F,stringsAsFactors = F)

res.save$SNP.adj <- substring(res.save$SNP,1,as.numeric(gregexpr('_',res.save$SNP))-1)

df <- merge(res.save,gwas1[,c('RSID','BETA','BP')],by.x='SNP.adj',by.y='RSID',all.x=T)
df <- merge(df,gwas2[,c('RSID','BETA','BP')],by.x='SNP.adj',by.y='RSID',all.x=T)
df$BETA <- df$BETA.x
df$BETA[is.na(df$BETA)] <- -1*df$BETA.y[is.na(df$BETA)]

df.sub <- subset(df,p_int<0.05)# & p_snp <0.05)
df.sub.sub <- subset(df.sub,BETA > 0)
binom.test(mean(sign(log(df.sub.sub$beta_int))==1)*nrow(df.sub.sub),nrow(df.sub.sub),0.5)
df.sub.sub <- subset(df.sub,BETA < 0)
binom.test(mean(sign(log(df.sub.sub$beta_int))==-1)*nrow(df.sub.sub),nrow(df.sub.sub),0.5)

f.out <- '/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.SNP_names.txt'
fwrite(data.frame(subset(res.save,FDR<0.1)$SNP.adj),f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
fwrite(df,'/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.txt',quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

###########

df <- fread("/Users/andrewmarderstein/Documents/Research/druggene/output/pgs_int/SNP_by_drug.txt",data.table = F,stringsAsFactors = F)
df.sub <- subset(df,p_int<0.05)# & p_snp <0.05)
df.sub.sub <- subset(df.sub,BETA > 0)
ggplot(subset(df.sub.sub,BETA!=0),aes(x=log(beta_int),fill=as.factor(sign(log(beta_int))))) + geom_histogram(col='black') + # geom_density(alpha=0.2)
  guides(fill=F) + scale_fill_discrete() + theme_bw() + theme(panel.grid = element_blank()) + labs(x='Observed interaction effect',y='Count')

library(qqman)
df$BP <- apply(df[,c('BP.x','BP.y')],1,min,na.rm=T)
df <- subset(df,!is.na(p_int))
manhattan(df,p='p_int',suggestiveline = -log10(max(subset(df,FDR<0.1)$p_int)))


########

library(data.table)
df <- fread("/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.txt",data.table = F,stringsAsFactors = F)
df.sub <- subset(df,FDR<0.1)
df.sub$ind[df.sub$ind==1] <- 'external'
df.sub$ind[df.sub$ind==2] <- 'internal'
df.sub$ind[df.sub$ind==2] <- 'both'
colnames(df.sub)[1] <- 'rs'
df.sub <- df.sub[,c(1:9,14)]
colnames(df.sub)[c(1,ncol(df.sub))] <- c('rs','BETA.marginal_GWAS')

f <- '/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.SNP_names.query.gene_name.txt'
genes <- fread(f,data.table = F,stringsAsFactors = F,header=F)
f <- '/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.SNP_names.query.genes.txt'
ensembl <- fread(f,data.table = F,stringsAsFactors = F,header=F)

df.sub$gene <- genes[,1]
df.sub$ensembl <- ensembl[,1]
fwrite(df.sub,'/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.sig_hit.supp_tab.txt',quote = F,na = 'NA',sep = '\t',row.names = F,col.names = T)


