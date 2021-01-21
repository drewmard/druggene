# Original script by Scott Kulm. Adapted by Andrew Marderstein.
# https://kulmsc.github.io/pgs_book/

# arguments
author <- 'michailidou'
origdir <- '/athena/elementolab/scratch/anm2868/druggene'
#####################
# load
library(data.table)
library("dplyr")
#####################

# initialize
ss.full <- fread(paste0(origdir,'/output/ss/',author,'/','raw_',tolower(author),'.ss.gz'),data.table = F,stringsAsFactors = F)
impute <- fread(paste0(origdir,'/output/ukb/impute_rsids'), header = F, stringsAsFactors = F,data.table = F)
colnames(impute) <- c("LOC", "SNP", "POS", "A1", "A2", "MAF", "AX", "INFO","CHR")
params <- read.table(paste0(origdir,'/output/ss/',author, "/parameters"), stringsAsFactors=F)

dir <- paste0(origdir,'/output/ss/',author)
system(paste0("mkdir -p ",dir))
system(paste0("echo Standardizing the summary statistics... > ", dir, "/clean.log"))
ss <- ss.full
system(paste0("echo original: ", nrow(ss), " >> ", dir, "/clean.log"))
#####################

# preprocess: 
# Check if there are missing columns, and fill them in just so things don't go haywire below
missing_cols <- params[!(params[,1] %in% colnames(ss)),1]
if(length(missing_cols) > 0){
  if("NO_CHR" %in% missing_cols){
    ss$NO_CHR <- 0
  }
  if("NO_POS" %in% missing_cols){
    ss$NO_POS <- 0
  }
  if("NO_SE" %in% missing_cols){
    ss$NO_SE <- 0
  }
}

#If there is a name in the parameters that states CHANGE_BOTH assuming that it's not BETA but OR
#This would just be a simple exponentiation
if(nrow(params) > 8){
  if(params[9,1] == "CHANGE_BOTH"){
    #will assume that the effect name is for odds ratio
    #following the walds ratio tests can just switch things over assuming a normal distribution
    #This may not be exactly what was completed in the published GWAS, but it should be a good approximation
    #And few methods seem to use the SE column anyway
    ss$BETA <- log(ss$BETA)
    ss$SE <- abs(ss$BETA/qnorm(ss$P))
  }
}
########################################

#Construct a common ss object, with common column names based on the manually curated params
ss <-  ss[,c(which(colnames(ss) == params[1,1]),
             which(colnames(ss) == params[2,1]),
             which(colnames(ss) == params[3,1]),
             which(colnames(ss) == params[4,1]),
             which(colnames(ss) == params[5,1]),
             which(colnames(ss) == params[6,1]),
             which(colnames(ss) == params[7,1]),
             which(colnames(ss) == params[8,1]))]
colnames(ss) <- c("CHR", "BP", "RSID", "A1", "A2", "SE", "BETA", "P")
ss$A1 <- toupper(ss$A1)
ss$A2 <- toupper(ss$A2)
ss <- ss[!is.na(ss$RSID),]
ss <- ss[!is.na(ss$BETA),]
ss <- ss[!is.na(ss$P),]

system(paste0("echo removed NA cols: ", nrow(ss), " >> ", dir, "/clean.log"))

#Remove duplicate SNP IDs
ss <- ss[!(duplicated(ss[,c('CHR','BP')]) |
  duplicated(ss[,c('CHR','BP')],fromLast = TRUE)),]
# unique(ss$RSID[duplicated(ss$RSID)])

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

#Make sure the imputation reference and the ss objects have the same number of SNP IDs
ss$RSID <- tolower(ss$RSID)
# ss <- ss[ss$RSID %in% impute$SNP,]
# sub_impute <- impute[impute$SNP %in% ss$RSID,]
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

###############

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

f.out <- paste0(dir, "/clean_", tolower(author), ".txt")
fwrite(ss, paste0(dir, "/clean_", tolower(author), ".txt"), row.names = F, col.names = T, sep = '\t', quote = F,na = 'NA')
system(paste0("gzip ", dir, "/clean_", tolower(author), ".txt"))
