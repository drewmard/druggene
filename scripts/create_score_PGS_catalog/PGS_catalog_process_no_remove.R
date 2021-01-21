library(data.table)
library("dplyr")

# conda activate DrugPGS

# parameters
author <- 'PGS000004'; use_col <- c(1:5)
author <- 'PGS000007'; use_col <- c(1:5)

origdir <- '/athena/elementolab/scratch/anm2868/druggene'

# initialize
ss.full <- fread(paste0(origdir,'/output/ss/',author,'/',(author),'.txt.gz'),data.table = F,stringsAsFactors = F)
impute <- fread(paste0(origdir,'/output/ukb/impute_rsids'), header = F, stringsAsFactors = F,data.table = F)
colnames(impute) <- c("LOC", "SNP", "POS", "A1", "A2", "MAF", "AX", "INFO","CHR")

dir <- paste0(origdir,'/output/ss/',author)
fullFile=paste0(origdir,'/output/ukb/ukb_mfi_chrAll_v3.txt')
f.out=paste0(origdir,'/output/ss/',author,'/for_grep')
fwrite(data.frame(paste0(ss.full$chr_name,':',ss.full$chr_position,'_')),f.out,col.names = F,row.names = F,quote = F,na = 'NA',sep = '\t')
outFile=paste0(origdir,'/output/ss/',author,'/ukb_info_rsid')
system(paste0('grep -f ',f.out,' ',fullFile,' > ',outFile))

impute <- fread(outFile, header = F, stringsAsFactors = F,data.table = F)
colnames(impute) <- c("LOC", "SNP", "POS", "A1", "A2", "MAF", "AX", "INFO","CHR")
ss <- ss.full[,use_col]
colnames(ss) <- c('CHR','BP','A2','A1','BETA')

#Make sure the imputation reference and the ss objects have the same number of SNP IDs
ss <- subset(ss,BP %in% impute$POS & CHR %in% impute$CHR)
sub_impute <- merge(impute,ss[,c('CHR','BP')],by.x=c('CHR','POS'),by.y=c('CHR','BP'))
ss <- merge(ss,sub_impute[,c('CHR','POS')],by.x=c('CHR','BP'),by.y=c('CHR','POS'))
ss$RSID <- sub_impute$SNP
# system(paste0("echo match reference to ss: ", nrow(ss), " >> ", dir, "/clean.log"))

flip_reverse <- function(ss,sub_impute) {
  flip_strand <- function(allele) {
    dplyr::case_when(
      allele == "A" ~ "T",
      allele == "C" ~ "G",
      allele == "T" ~ "A",
      allele == "G" ~ "C",
      TRUE ~ NA_character_
    )
  }
  i <- which((ss$A1 == "A" & ss$A2 == "T") | (ss$A1 == "T" & ss$A2 == "A") | (ss$A1 == "G" & ss$A2 == "C") | (ss$A1 == "C" & ss$A2 == "G"))
  if (length(i) > 0) {
    ss.tmp <- ss[-i,]
    ss.tmp2 <- ss[i,]
  } else {
    ss.tmp <- ss
  }
  ss2 <- ss.tmp
  ss3 <- ss.tmp
  ss3$A1 <- flip_strand(ss2$A1)
  ss3$A2 <- flip_strand(ss2$A2)
  ss3 <- rbind(ss2, ss3) #####
  ss4 <- ss3
  ss4$A1 <- ss3$A2
  ss4$A2 <- ss3$A1
  ss4$BETA <- -ss3$BETA
  ss4 <- rbind(ss3, ss4) ######
  if (length(i) > 0) { # for ambiguous strands, dont flip --> just reverse
    ss.tmp2.reverse <- ss.tmp2
    ss.tmp2.reverse$A1 <- ss.tmp2$A2
    ss.tmp2.reverse$A2 <- ss.tmp2$A1
    ss.tmp2.reverse$BETA <- -ss.tmp2$BETA
    ss.tmp2.reverse <- rbind(ss.tmp2, ss.tmp2.reverse) ######
    ss4 <- rbind(ss4,ss.tmp2.reverse)
  }
  # ss.matched <- merge(as.data.table(ss4), as.data.table(sub_impute[,c('CHR','POS','A1','A2')]), 
  #                  by.x = c('CHR','BP','A1','A2'),by.y=c('CHR','POS','A1','A2'), all = FALSE)
  ss.matched <- merge(ss4, sub_impute[,c('CHR','POS','A1','A2')], 
                      by.x = c('CHR','BP','A1','A2'),by.y=c('CHR','POS','A1','A2'), all = FALSE)
  return(ss.matched)
}
ss <- flip_reverse(ss,sub_impute)

ss <- ss[,c('CHR','BP','A2','A1','BETA','RSID')]
fwrite(ss, paste0(dir, "/clean_", author, ".no_remove.txt"), row.names = F, col.names = T, sep = '\t', quote = F,na = 'NA')
system(paste0("gzip -f ", dir, "/clean_", author, ".no_remove.txt"))
fwrite(as.data.frame(ss$RSID),paste0(dir, "/", 'rsids.no_remove'),row.names = F, col.names = F, sep = '\t', quote = F,na = 'NA')

