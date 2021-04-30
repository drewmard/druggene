library(data.table)
library(ggplot2)
res <- fread("/Users/andrewmarderstein/Documents/Research/druggene/output/pgs_int/snp_drug_enrichr_ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt",data.table = F,stringsAsFactors = F)
res$Term2 <- substring(res[,1],1,as.numeric(gregexpr(' ',res[,1]))-1)
res$db <- substring(res[,1],as.numeric(gregexpr(' ',res[,1]))+1)

g1 <- ggplot(res[1:10,],aes(x=reorder(Term2,-log10(Adjusted.P.value)),y=-log10(Adjusted.P.value),fill=db)) + geom_bar(stat='identity') +
  theme_bw() + theme(panel.grid = element_blank(),axis.text.x = element_text(angle=30,hjust=1)) +
  labs(x='Transcription Factor',y='Significance (-log10 FDR)') + scale_fill_manual(values=c('steelblue4','lightblue')) + geom_hline(yintercept = 1,col='red',lty='dashed'); g1

bp <- fread("/Users/andrewmarderstein/Documents/Research/druggene/output/pgs_int/snp_drug_enrichr_GO_Biological_Process_2018.txt",data.table = F,stringsAsFactors = F)
mf <- fread("/Users/andrewmarderstein/Documents/Research/druggene/output/pgs_int/snp_drug_enrichr_GO_Molecular_Function_2018.txt",data.table = F,stringsAsFactors = F)
bp$lib <- "BP"
mf$lib <- "MF"
df.mg <- rbind(bp,mf)
df.mg <- subset(df.mg,Adjusted.P.value < 0.1)
g2 <- ggplot(df.mg,aes(x=reorder(Term,-log10(Adjusted.P.value)),y=-log10(Adjusted.P.value),fill=lib)) + geom_bar(stat='identity') +
  theme_bw() + theme(panel.grid = element_blank(),axis.text.x = element_text(angle=-90,hjust=0)) +
  labs(x='Gene ontology',y='-log10 FDR') + scale_fill_manual(values=c('chocolate4','chocolate1')); g2

library(cowplot)
plot_grid(g1,g2,ncol=1)
g2
