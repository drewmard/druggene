library(data.table)
ss <- fread('/athena/elementolab/scratch/anm2868/druggene/output/ss/schumacher/29892016-GCST006085-EFO_0001663.h.tsv.gz',data.table=F,stringsAsFactors = F)
liftover <- fread('/athena/elementolab/scratch/anm2868/druggene/output/ss/schumacher/liftOver.output',data.table = F,stringsAsFactors = F)
liftover <- unique(liftover)
# ss.sub <- ss[,c('hm_rsid','hm_chrom','hm_pos','hm_beta','p_value')]
df.mg <- merge(ss[,c('hm_rsid','hm_chrom','hm_other_allele','hm_effect_allele','hm_beta','p_value')],liftover[,c(4,3)],by.x='hm_rsid',by.y='V4')
colnames(df.mg)[colnames(df.mg)=='V3'] <- 'pos_hg19'

origdir <- '/athena/elementolab/scratch/anm2868/druggene'
author <- 'schumacher'
f <- paste0(origdir,'/output/ss/',author,'/','raw_',tolower(author),'.ss.gz')
fwrite(df.mg,f,row.names = F, col.names = T, sep = '\t', quote = F,na = 'NA')

#######################


ss.full <- df.mg

impute <- fread(paste0(origdir,'/output/ukb/impute_rsids'), header = F, stringsAsFactors = F,data.table = F)
colnames(impute) <- c("LOC", "SNP", "POS", "A1", "A2", "MAF", "AX", "INFO","CHR")

dir <- paste0(origdir,'/output/ss/',author)
system(paste0("mkdir -p ",dir))
system(paste0("echo Standardizing the summary statistics... > ", dir, "/clean.log"))
ss <- ss.full
system(paste0("echo original: ", nrow(ss), " >> ", dir, "/clean.log"))

ss <- ss[,c("hm_chrom",
            "pos_hg19",
            'hm_rsid',
            'hm_other_allele',
            'hm_effect_allele',
            'hm_beta',
            'p_value')]
colnames(ss) <- c("CHR", "BP", "RSID", "A1", "A2", "BETA", "P")

ss$A1 <- toupper(ss$A1)
ss$A2 <- toupper(ss$A2)
ss <- ss[!is.na(ss$RSID),]
ss <- ss[!is.na(ss$BETA),]
ss <- ss[!is.na(ss$P),]

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

