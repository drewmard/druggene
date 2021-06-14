# after running:

library(parallel)
library(data.table)
library(rlang)
library(httr)
library(stringr)
source('/athena/elementolab/scratch/anm2868/open_targets/scripts/open_targets_query_func.R')
query <- fread("/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.SNP_names.query.txt",data.table = F,stringsAsFactors = F,header = F,fill=T)
f <- '/home/anm2868/data/protein_coding_genes.txt'
protein_coding_genes <- fread(f,data.table = F,stringsAsFactors = F)


# need to create query file:

assign_genes_to_variants2 <- function(i,prioritize_context=TRUE) {
  print(i)
  snp <- query[i,2]
  res <- genesForVariant2(snp)
  if (length(res)==0) { # no genes for variant through v2g algorithm
    res <- nearestCodingGene(snp)
    res.df <- data.frame(rs=query[i,1],SNP=query[i,2],Gene=res$symbol,Ensembl=res$id,Score=100,func=0,stringsAsFactors=FALSE)
  } else {
    res.df <- do.call(rbind,lapply(res,function(x) {
      funcVal=ifelse(length(x$functionalPredictions)>0,x$functionalPredictions[[1]]$aggregatedScore,0)
      data.frame(rs=query[i,1],SNP=query[i,2],Gene=x$gene$symbol,Ensembl=x$gene$id,Score=x$overallScore,func=funcVal,stringsAsFactors=FALSE)}))
    if (prioritize_context) {
      res.df <- res.df[order(res.df$func,res.df$Score,decreasing = T),]
    } else {
      res.df <- res.df[order(res.df$Score,decreasing = T),]
    }
    res.df <- subset(res.df,Gene %in% protein_coding_genes[,'Gene name'])
    if (nrow(res.df)==0) {
      res <- nearestCodingGene(snp)
      res.df <- data.frame(rs=query[i,1],SNP=query[i,2],Gene=res$symbol,Ensembl=res$id,Score=100,func=0,stringsAsFactors=FALSE)
    }
  }
  return(res.df)
}

res.df.lst <- lapply(1:nrow(query),assign_genes_to_variants2)
res.df.save <- do.call(rbind,res.df.lst)
res.df.save.sub <- by(res.df.save,res.df.save$rs,head,1)
res.df.save.sub <- Reduce(rbind,res.df.save.sub)
f.out <- '/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.SNP_names.query.genes.txt'
fwrite(as.data.frame(res.df.save.sub$Ensembl),f.out,sep = '\t',quote = F,na='NA',row.names = F,col.names = F)
length(unique(res.df.save.sub$Ensembl))
#40 gwas hits mapping to 35 unique genes
f.out <- '/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.SNP_names.query.gene_name.txt'
fwrite(as.data.frame(res.df.save.sub$Gene),f.out,sep = '\t',quote = F,na='NA',row.names = F,col.names = F)

##############

# res.df.lst <- lapply(1:nrow(query),assign_genes_to_variants2,prioritize_context=FALSE)
# res.df.save <- do.call(rbind,res.df.lst)
# res.df.save.sub <- by(res.df.save,res.df.save$rs,head,1)
# res.df.save.sub <- Reduce(rbind,res.df.save.sub)
# f.out <- '/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.SNP_names.query_nocontext.genes.txt'
# fwrite(as.data.frame(res.df.save.sub$Ensembl),f.out,sep = '\t',quote = F,na='NA',row.names = F,col.names = F)
# length(unique(res.df.save.sub$Ensembl))
# #40 gwas hits mapping to 35 unique genes
# f.out <- '/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.SNP_names.query_nocontext.gene_name.txt'
# fwrite(as.data.frame(res.df.save.sub$Gene),f.out,sep = '\t',quote = F,na='NA',row.names = F,col.names = F)

########################

library(data.table)
f.out <- '/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.SNP_names.query.gene_name.txt'
res.df.save.sub <- fread(f.out,data.table = F,stringsAsFactors = F,header = F)
library(enrichR)
dbs <- listEnrichrDbs()
enriched <- enrichr(res.df.save.sub$Gene, dbs$libraryName)
res <- list()
for (data_base in dbs$libraryName) {
  res[[data_base]] <- subset(enriched[[data_base]],Adjusted.P.value<0.1)
  if (nrow(res[[data_base]])) {
    res[[data_base]]$db <- data_base
  }
}
res2 <- do.call(rbind,res)
res2[order(res2$Adjusted.P.value)[1:5],]

data_base="ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"
tmp <- enriched[[data_base]]
f.out <- paste0("/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/snp_drug_enrichr_",data_base,".txt")
fwrite(tmp,f.out,row.names = F,col.names = T,sep = '\t',na = "NA",quote=F)

data_base="DSigDB"
enriched <- as.data.frame(enrichr(res.df.save.sub[,1], data_base))
library(stringr); enriched$num_gene <- str_count(enriched[,ncol(enriched)],';')+1
tmp <- enriched[enriched$num_gene>=3,]
grep('dexa',enriched[,1])

f.out <- paste0("/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/snp_drug_enrichr_",data_base,".txt")
fwrite(tmp,f.out,row.names = F,col.names = T,sep = '\t',na = "NA",quote=F)

data_base="GO_Biological_Process_2018"
tmp <- enriched[[data_base]]
library(stringr); num_gene <- str_count(tmp[,ncol(tmp)],';')+1
tmp <- tmp[num_gene>=3,]
subset(tmp,Adjusted.P.value < 0.1)
f.out <- paste0("/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/snp_drug_enrichr_",data_base,".txt")
fwrite(tmp,f.out,row.names = F,col.names = T,sep = '\t',na = "NA",quote=F)
data_base="GO_Molecular_Function_2018"
tmp <- enriched[[data_base]]
library(stringr); num_gene <- str_count(tmp[,ncol(tmp)],';')+1
tmp <- tmp[num_gene>=3,]
subset(tmp,Adjusted.P.value < 0.1)
f.out <- paste0("/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/snp_drug_enrichr_",data_base,".txt")
fwrite(tmp,f.out,row.names = F,col.names = T,sep = '\t',na = "NA",quote=F)

