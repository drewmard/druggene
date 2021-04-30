library(data.table)
f <- paste0("/Users/andrewmarderstein/Documents/Research/druggene/output/snp_int/nrf2_gr_locus.txt")
# f <- paste0("/Users/andrewmarderstein/Documents/Research/druggene/output/snp_int/chr19_locus.txt")
df <- fread(f,data.table = F,stringsAsFactors = F)
attr <- fread("~/Downloads/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",data.table = F,stringsAsFactors = F)
df.mg <- merge(df,attr[,c('SAMPID','SMTSD')],by.x='gene_names',by.y='SAMPID')
df.mg$IID <- unlist(lapply(strsplit(df.mg$gene_names,'-'),function(x) paste(x[1:2],collapse = '-')))
attr2 <- fread("~/Downloads/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",data.table = F,stringsAsFactors = F)
df.mg <- merge(df.mg,attr2[,c('SUBJID','SEX')],by.x='IID',by.y='SUBJID')

res.lst <- list()
for (tis in unique(df.mg$SMTSD)) {
  df.sub <- subset(df.mg,SMTSD==tis)
  res <- cor.test(df.sub[,3],df.sub[,4])
  res.lst[[tis]] <- data.frame(tis,n=nrow(df.sub),cor=as.numeric(res$estimate),p=as.numeric(res$p.value))
}
res.save <- do.call(rbind,res.lst); row.names(res.save) <- NULL
res.save <- subset(res.save,n>70)
res.save$fdr <- p.adjust(res.save$p,method='fdr')
res.save[order(res.save$cor,decreasing = F),][1:10,]
nrow(subset(res.save,fdr<0.1))/nrow(res.save)
nrow(subset(res.save,fdr<0.1 & cor>0))/nrow(res.save)
nrow(subset(res.save,fdr<0.1 & cor<0))/nrow(res.save)
(subset(res.save,fdr>0.1))
df.sub <- subset(df.mg,SMTSD=="Breast - Mammary Tissue" & SEX==1)
cor.test(df.sub[,3],df.sub[,4])
df.sub <- subset(df.mg,SMTSD=="Breast - Mammary Tissue" & SEX==2)
cor.test(df.sub[,3],df.sub[,4])


cor.diff.test = function(x1, x2, y1, y2, method="pearson") {
  cor1 = cor.test(x1, x2, method=method)
  cor2 = cor.test(y1, y2, method=method)
  
  r1 = cor1$estimate
  r2 = cor2$estimate
  n1 = sum(complete.cases(x1, x2))
  n2 = sum(complete.cases(y1, y2))
  fisher = ((0.5*log((1+r1)/(1-r1)))-(0.5*log((1+r2)/(1-r2))))/((1/(n1-3))+(1/(n2-3)))^0.5
  
  p.value = (2*(1-pnorm(abs(fisher))))
  
  result= list(
    "cor1" = list(
      "estimate" = as.numeric(cor1$estimate),
      "p.value" = cor1$p.value,
      "n" = n1
    ),
    "cor2" = list(
      "estimate" = as.numeric(cor2$estimate),
      "p.value" = cor2$p.value,
      "n" = n2
    ),
    "p.value.twosided" = as.numeric(p.value),
    "p.value.onesided" = as.numeric(p.value) / 2
  )
  cat(paste(sep="",
            "cor1: r=", format(result$cor1$estimate, digits=3), ", p=", format(result$cor1$p.value, digits=3), ", n=", result$cor1$n, "\n",
            "cor2: r=", format(result$cor2$estimate, digits=3), ", p=", format(result$cor2$p.value, digits=3), ", n=", result$cor2$n, "\n",
            "diffence: p(one-sided)=", format(result$p.value.onesided, digits=3), ", p(two-sided)=", format(result$p.value.twosided, digits=3), "\n"
  ))
  return(result);
}


df.sub1 <- subset(df.mg,SMTSD=="Breast - Mammary Tissue" & SEX==1)
# cor.test(df.sub1[,3],df.sub1[,4])
df.sub2 <- subset(df.mg,SMTSD=="Breast - Mammary Tissue" & SEX==2) # 2 is female
# cor.test(df.sub2[,3],df.sub2[,4])
cor.diff.test(df.sub1[,3],df.sub1[,4],df.sub2[,3],df.sub2[,4])

library(ggplot2)
ggplot(subset(df.mg,SMTSD=="Breast - Mammary Tissue"),aes(x=NR3C1,y=NFE2L2,col=as.factor(SEX))) + 
  geom_point(alpha=0.1) +
  geom_smooth(method='lm',se=F) + theme_bw() + theme(panel.grid = element_blank()) +
  guides(col=F) + labs(x='GR gene expression (TPM)',y='NRF2 gene expression (TPM)') + scale_color_manual(values = c('black','red'))




















# create plots
g1 <- ggplot(df.mg,aes(x=reorder(SMTSD,`EI Score`,FUN=median),y=`EI Score`)) + 
  geom_point(col='orange',alpha=0.6) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  scale_y_continuous(labels=function(x) sprintf("%.2f", x)) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1,size = rel(0.7)),
        panel.grid = element_blank()) +
  labs(x='Tissue Type',y="EI Score"); g1
