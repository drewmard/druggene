library(data.table)
library(stringr)
f <- '/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.txt'
df <- fread(f,data.table = F,stringsAsFactors = F)
# x <- strsplit(df$SNP,'_')[1:5]
df$effect_allele <- unlist(lapply(strsplit(df$SNP,'_'),function(x) x[[2]]))
df$BP <- apply(df[,c('BP.x','BP.y')],1,min,na.rm=T)
df <- df[,c('ind','CHR','BP','SNP.adj',
            'beta_snp','p_snp','beta_int','p_int','FDR','BETA.x','BETA.y')]
colnames(df) <- 
  c('Study','CHR','BP','RSID',
    'HR_SNP','P_SNP','HR_INT','P_INT','FDR_INT',
    'BETA_Mav','BETA_Mich')
fwrite(df,'/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.SUPP.txt',quote = F,na = "NA",sep = '\t',col.names = T,row.names = F)