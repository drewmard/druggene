library(data.table)
library("dplyr")

# parameters
author <- 'PGS000004'; use_col <- c(1:5)
# author <- 'PGS000007'; use_col <- c(1:5)
# author <- 'PGS000008'
# author <- 'PGS000009'
# author <- 'PGS000082'; use_col <- c(2:6)

origdir <- '/athena/elementolab/scratch/anm2868/druggene'

# initialize
ss.full <- fread(paste0(origdir,'/output/ss/',author,'/',(author),'.txt.gz'),data.table = F,stringsAsFactors = F)
# impute <- fread(paste0(origdir,'/output/ukb/impute_rsids'), header = F, stringsAsFactors = F,data.table = F)
# colnames(impute) <- c("LOC", "SNP", "POS", "A1", "A2", "MAF", "AX", "INFO","CHR")

dir <- paste0(origdir,'/output/ss/',author)
system(paste0("mkdir -p ",dir))
system(paste0("echo Standardizing the summary statistics... > ", dir, "/clean.log"))
ss <- ss.full[,use_col]
system(paste0("echo original: ", nrow(ss), " >> ", dir, "/clean.log"))
colnames(ss) <- c('CHR','BP','A2','A1','BETA')

# Remove NA columns
ss <- ss[!is.na(ss$CHR),]
ss <- ss[!is.na(ss$BP),]
ss <- ss[!is.na(ss$BETA),]
system(paste0("echo removed NA cols: ", nrow(ss), " >> ", dir, "/clean.log"))
#Remove duplicate SNP IDs
ss <- ss[!(duplicated(ss[,c('CHR','BP')]) |
             duplicated(ss[,c('CHR','BP')],fromLast = TRUE)),]
system(paste0("echo removed duplicated IDs: ", nrow(ss), " >> ", dir, "/clean.log"))
# only SNV. no indels. must be A, C, G, or T.
ss <- ss[ss$A1 %in% c("A", "C", "G", "T") & ss$A2 %in% c("A", "C", "G", "T"),]
system(paste0("echo removed indels and non ACGT SNV: ", nrow(ss), " >> ", dir, "/clean.log"))
# ambiguous A/T or C/G removed
ss <- ss[!((ss$A1 == "A" & ss$A2 == "T") |
             (ss$A1 == "T" & ss$A2 == "A") |
             (ss$A1 == "G" & ss$A2 == "C") |
             (ss$A1 == "C" & ss$A2 == "G")),]
system(paste0("echo removed ambiguous A/T or C/G SNPs: ", nrow(ss), " >> ", dir, "/clean.log"))

# head(subset(ss,!(BP %in% impute$POS & CHR %in% impute$CHR)))
#Make sure the imputation reference and the ss objects have the same number of SNP IDs
ss <- subset(ss,BP %in% impute$POS & CHR %in% impute$CHR)
sub_impute <- subset(impute,POS %in% ss$BP & CHR %in% ss$CHR & 
                       A1 %in% c('A','C','G','T') &
                       A2 %in% c('A','C','G','T') &
                       !((A1 == "A" & A2 == "T") |
                           (A1 == "T" & A2 == "A") |
                           (A1 == "G" & A2 == "C") |
                           (A1 == "C" & A2 == "G")))
sub_impute <- sub_impute[!(
  duplicated(sub_impute[,c('CHR','POS')]) |
    duplicated(sub_impute[,c('CHR','POS')],fromLast = TRUE)
),]
sub_impute <- merge(sub_impute,ss[,c('CHR','BP')],by.x=c('CHR','POS'),by.y=c('CHR','BP'))
ss <- merge(ss,sub_impute[,c('CHR','POS')],by.x=c('CHR','BP'),by.y=c('CHR','POS'))
ss$RSID <- sub_impute$SNP
system(paste0("echo match reference to ss: ", nrow(ss), " >> ", dir, "/clean.log"))

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
  ss2 <- ss
  ss3 <- ss
  ss3$A1 <- flip_strand(ss2$A1)
  ss3$A2 <- flip_strand(ss2$A2)
  ss3 <- rbind(ss2, ss3) #####
  ss4 <- ss3
  ss4$A1 <- ss3$A2
  ss4$A2 <- ss3$A1
  ss4$BETA <- -ss3$BETA
  ss4 <- rbind(ss3, ss4) ######
  # ss.matched <- merge(as.data.table(ss4), as.data.table(sub_impute[,c('CHR','POS','A1','A2')]), 
  #                  by.x = c('CHR','BP','A1','A2'),by.y=c('CHR','POS','A1','A2'), all = FALSE)
  ss.matched <- merge(ss4, sub_impute[,c('CHR','POS','A1','A2')], 
                      by.x = c('CHR','BP','A1','A2'),by.y=c('CHR','POS','A1','A2'), all = FALSE)
  return(ss.matched)
}
ss <- flip_reverse(ss,sub_impute)

system(paste0("echo alleles: flipping, reversing, and removing mismatches ss: ", nrow(ss), " >> ", dir, "/clean.log"))
# ss <- ss[,c('CHR','BP','A2','A1','BETA','A2_AF','RSID')]
ss <- ss[,c('CHR','BP','A2','A1','BETA','RSID')]
fwrite(ss, paste0(dir, "/clean_", author, ".txt"), row.names = F, col.names = T, sep = '\t', quote = F,na = 'NA')
system(paste0("gzip -f ", dir, "/clean_", author, ".txt"))
fwrite(as.data.frame(ss$RSID),paste0(dir, "/", 'rsids'),row.names = F, col.names = F, sep = '\t', quote = F,na = 'NA')

