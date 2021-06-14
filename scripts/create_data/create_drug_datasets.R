conda activate DrugPGS
spack load /mjrrusu
R

library(data.table)
library(parallel)
pheno <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)
drug_data_names <- colnames(pheno)[grep("20003-",colnames(pheno))]
x <- (pheno[,drug_data_names])
unique_drug_data_names <- unique(as.vector(as.matrix(pheno[,drug_data_names])))
unique_drug_data_names <- unique_drug_data_names[!is.na(unique_drug_data_names)]
drugs_per_individual <- as.character(lapply(1:nrow(x),function(i) {paste(x[i,which(!is.na(x[i,]))],collapse = ',')}))

res <- mclapply(1:length(unique_drug_data_names),function(i) grepl(unique_drug_data_names[i],drugs_per_individual),mc.cores = 8)
names(res) <- unique_drug_data_names
res.df <- as.data.frame(res)
# res.df2 <- as.data.frame(apply(res.df,2,as.numeric))

# res <- mclapply(1:ncol(res.df),function(i) sum(res.df[,i]),mc.cores = 4)
# res <- as.numeric(res)
# ind <- order(res,decreasing = T)
# unique_drug_data_names[ind[1]]

f <- '/athena/elementolab/scratch/anm2868/DrugPGS/Wu_et_al_2019_data.csv'
Wu_et_al <- fread(f,data.table = F,stringsAsFactors = F)[,1:4]
Wu_et_al <- Wu_et_al[-nrow(Wu_et_al),]
ATC_codes <- (lapply(1:nrow(Wu_et_al),function(i) strsplit(Wu_et_al[i,3],' \\|')))
ATC_codes.unlist <- unlist(ATC_codes)

##########

level2.uniq <- unique(substring(ATC_codes.unlist,1,3))
ATC_codes.level2 <- lapply(ATC_codes,function(s) {substring(unlist(s),1,3)})
level2.res <- lapply(1:length(level2.uniq),function(i) {
  x <- Wu_et_al[grep(level2.uniq[i],ATC_codes.level2),2]
  if (length(x)>0) {
    paste0('X',x)
  } else {
    x
  }
})
names(level2.res) <- level2.uniq

level2.res2 <- mclapply(1:length(level2.res),function(i) {
  drug_codes <- unlist(level2.res[i])
  if (length(drug_codes)>1) {
    unname(apply(res.df[,drug_codes],1,sum)>0)
  } else {
    res.df[,drug_codes]
  }
},mc.cores = 4)
names(level2.res2) <- level2.uniq
level2.res2.df <- as.data.frame(level2.res2)

##########
level3.uniq <- unique(substring(ATC_codes.unlist,1,4))
ATC_codes.level3 <- lapply(ATC_codes,function(s) {substring(unlist(s),1,4)})
level3.res <- lapply(1:length(level3.uniq),function(i) {
  x <- Wu_et_al[grep(level3.uniq[i],ATC_codes.level3),2]
  if (length(x)>0) {
    paste0('X',x)
  } else {
    x
  }
})
names(level3.res) <- level3.uniq

level3.res2 <- mclapply(1:length(level3.res),function(i) {
  drug_codes <- unlist(level3.res[i])
  if (length(drug_codes)>1) {
    unname(apply(res.df[,drug_codes],1,sum)>0)
  } else {
    res.df[,drug_codes]
  }
},mc.cores = 4)
names(level3.res2) <- level3.uniq
level3.res2.df <- as.data.frame(level3.res2)
#######3

level4.uniq <- unique(substring(ATC_codes.unlist,1,5))
level4.uniq <- level4.uniq[nchar(level4.uniq)==5]
ATC_codes.level4 <- lapply(ATC_codes,function(s) {substring(unlist(s),1,5)})
level4.res <- lapply(1:length(level4.uniq),function(i) {
  x <- Wu_et_al[grep(level4.uniq[i],ATC_codes.level4),2]
  if (length(x)>0) {
    paste0('X',x)
  } else {
    x
  }
})
names(level4.res) <- level4.uniq

level4.res2 <- mclapply(1:length(level4.res),function(i) {
  drug_codes <- unlist(level4.res[i])
  if (length(drug_codes)>1) {
    unname(apply(res.df[,drug_codes],1,sum)>0)
  } else {
    res.df[,drug_codes]
  }
},mc.cores = 4)
level4.res2 <- mclapply(1:length(level4.res),function(i) {
  drug_codes <- unlist(level4.res[i])
  if (length(drug_codes)>1) {
    unname(apply(res.df[,drug_codes],1,sum)>0)
  } else {
    res.df[,drug_codes]
  }
},mc.cores = 2)
names(level4.res2) <- level4.uniq
level4.res2.df <- as.data.frame(level4.res2)

#########

res.df$eid <- pheno$eid
level2.res2.df$eid <- pheno$eid
level3.res2.df$eid <- pheno$eid
level4.res2.df$eid <- pheno$eid

res.df[,1:(ncol(res.df)-1)] <- apply(res.df[,1:(ncol(res.df)-1)],2,as.numeric)
level4.res2.df[,1:(ncol(level4.res2.df)-1)] <- apply(level4.res2.df[,1:(ncol(level4.res2.df)-1)],2,as.numeric)

level4.res2.df[,1:(ncol(level4.res2.df)-1)] <- apply(level4.res2.df[,1:(ncol(level4.res2.df)-1)],2,as.numeric)

f <- '/athena/elementolab/scratch/anm2868/DrugPGS/UKB_drug_data/orig_drug.txt'
fwrite(res.df,f,quote = F,sep = '\t',na = NA,row.names = F,col.names = T)
f <- '/athena/elementolab/scratch/anm2868/DrugPGS/UKB_drug_data/level2.txt'
fwrite(level2.res2.df,f,quote = F,sep = '\t',na = NA,row.names = F,col.names = T)
f <- '/athena/elementolab/scratch/anm2868/DrugPGS/UKB_drug_data/level3.txt'
fwrite(level3.res2.df,f,quote = F,sep = '\t',na = NA,row.names = F,col.names = T)
f <- '/athena/elementolab/scratch/anm2868/DrugPGS/UKB_drug_data/level4.txt'
fwrite(level4.res2.df,f,quote = F,sep = '\t',na = NA,row.names = F,col.names = T)

level2_cases <- apply(level2.res2.df,2,sum)
level3_cases <- apply(level3.res2.df,2,sum)
level4_cases <- apply(level3.res2.df,2,sum)

sort(level3_cases)
apply(level3.res2.df,2,sum)

