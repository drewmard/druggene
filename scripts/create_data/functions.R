
# required to specify:
# disease_score <- 'prostate_cancer'
# icd9_code <- c('185')
# icd10_code <- c('C61',paste0('C61',0:9))
# cancer_code <- TRUE; self_report_code <- '1044'
# doctor_diag_avail <- FALSE; doctor_diag_code <- NA
# female_only <- FALSE
# opcs_code_avail <- FALSE; opcs_code <- c('')


read_in_analysis_data <- function() {
  polyscore <- fread('/athena/elementolab/scratch/anm2868/DrugPGS/geno_data/andrew_w_eid.txt',data.table = F,stringsAsFactors = F)
  drug_data <- fread('/athena/elementolab/scratch/anm2868/DrugPGS/UKB_drug_data/level4.txt',data.table = F,stringsAsFactors = F)
}
# 
# read_in_phenotype_data <- function() {
  hesin_diag <- fread('/home/kulmsc/athena/ukbiobank/hesin/hesin_diag.txt',data.table = F,stringsAsFactors = F)
  hesin_oper <- fread('/home/kulmsc/athena/ukbiobank/hesin/hesin_oper.txt',data.table = F,stringsAsFactors = F)
  hesin <- fread('/home/kulmsc/athena/ukbiobank/hesin/hesin.txt',data.table = F,stringsAsFactors = F)
  pheno1 <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb26867.csv.gz',data.table=F,stringsAsFactors = F)
# }

create_phenotype_file <- function() {
  hesin_diag.sub <- subset(hesin_diag,diag_icd9 %in% icd9_code | diag_icd10 %in% icd10_code)
  tmp <- aggregate(hesin_diag.sub$ins_index,list(eid=hesin_diag.sub$eid),min)
  hesin_tmp <- merge(hesin[,c('eid','ins_index','epistart','admidate')],tmp,by.x=c('eid','ins_index'),by.y=c('eid','x'))
  hesin_tmp$epistart[is.na(hesin_tmp$epistart)] <- hesin_tmp$admidate[is.na(hesin_tmp$epistart)]
  hesin_tmp$case <- 1
  pheno_tmp <- pheno1[,c("eid","40000-0.0","53-0.0",'26410-0.0','26426-0.0','26427-0.0','191-0.0')]
  pheno_tmp2 <- merge(pheno_tmp,hesin_tmp,by="eid",all.x=TRUE)
  pheno_tmp2$case[is.na(pheno_tmp2$case)] <- 0
  
  if (opcs_code_avail) {
    hesin_oper.sub <- hesin_oper[grep(opcs_code,hesin_oper$oper4),]
    hesin_oper.sub$oper <- 1
    pheno_tmp2 <- merge(pheno_tmp2,hesin_oper.sub[,c('eid','opdate','oper')],by="eid",all.x=TRUE)
    pheno_tmp2$oper[is.na(pheno_tmp2$oper)] <- 0
  } else {
    pheno_tmp2$opdate <- NA
    pheno_tmp2$oper <- 0
  }
  
  censor_locat <- which(!is.na(pheno_tmp2[,c('26410-0.0','26426-0.0','26427-0.0')]),arr.ind = TRUE)
  pheno_tmp2$Hosp <- NA
  pheno_tmp2$Hosp[censor_locat[censor_locat[,2]==1,1]] <- 'HES'
  pheno_tmp2$Hosp[censor_locat[censor_locat[,2]==2,1]] <- 'PEDW'
  pheno_tmp2$Hosp[censor_locat[censor_locat[,2]==3,1]] <- 'SMR'
  tmp <- subset(hesin,eid %in% subset(pheno_tmp2,is.na(Hosp))$eid)
  tmp2 <- (aggregate(tmp$dsource,list(eid=tmp$eid),function(x) x[1]))
  pheno_tmp2$hesin <- as.character(tmp2$x[match(pheno_tmp2$eid,tmp2$eid)])
  i <- which(!is.na(pheno_tmp2$hesin))
  pheno_tmp2$Hosp[i] <- pheno_tmp2$hesin[i]
  pheno_tmp2$End <- as.Date(NA)
  pheno_tmp2$End[pheno_tmp2$Hosp=='HES'] <- as.Date("2020-03-31") # England_end
  pheno_tmp2$End[pheno_tmp2$Hosp=='PEDW'] <- as.Date("2016-02-29")# Wales_end
  pheno_tmp2$End[pheno_tmp2$Hosp=='SMR'] <- as.Date("2016-10-31") # Scotland_end
  pheno_tmp2$End[is.na(pheno_tmp2$End)] <- as.Date("2020-03-31") 
  pheno_tmp2$diagnosis <- dmy(pheno_tmp2$epistart)
  pheno_tmp2$opdat <- dmy(pheno_tmp2$opdate)
  
  df <- pheno_tmp2[,c('eid','53-0.0','case','oper','diagnosis','40000-0.0','191-0.0','End','opdat')]
  colnames(df) <- c('eid','start','case','oper','diagnosis','death','lost','end','opdate')
  df$death[df$death==''] <- NA
  df$lost[df$lost==''] <- NA
  df$final <- apply(df[,c('diagnosis','death','lost','end','opdate')],1,min,na.rm=T)
  df$days <- apply(df[,c('final','start')],1,function(x) as.numeric(difftime(as.Date(x[1]),as.Date(x[2]))))
  df$death_status <- as.numeric(!is.na(df$death))
  
  if (cancer_code) {
    self_report_field <- '20001-0'
  } else {
    self_report_field <- '20002-0'
  }
  ind <- unique(which(Reduce(`|`, 
                             lapply(self_report_code, `==`, pheno1[,grep(self_report_field,colnames(pheno1))])
  ), arr.ind = TRUE)[,1])
  self_reported <- rep(0,nrow(pheno1)); self_reported[ind] <- 1
  if (doctor_diag_avail) {
    doctor_diag <- as.numeric(pheno1[,c(doctor_diag_code)]==1)
    x <- as.numeric(apply(data.frame(self_reported,doctor_diag),1,sum,na.rm=T)>0)
    x2 <- data.frame(eid=pheno1$eid,self_reported=x)
  } else if (disease_score=='stroke') {
    doctor_diag <- as.numeric(apply(pheno1[,c('6150-0.0','6150-0.1','6150-0.2','6150-0.3')],1,function(x){3 %in% x}))
    x <- as.numeric(apply(data.frame(self_reported,doctor_diag),1,sum,na.rm=T)>0)
    x2 <- data.frame(eid=pheno1$eid,self_reported=x)
  } else if (disease_score=='asthma') {
    doctor_diag1 <- as.numeric(pheno1[,c('22127-0.0')]==1)
    doctor_diag2 <- as.numeric(apply(pheno1[,c('6152-0.0','6152-0.1','6152-0.2','6152-0.3','6152-0.4')],1,function(x){8 %in% x}))
    x <- as.numeric(apply(data.frame(self_reported,doctor_diag1,doctor_diag2),1,sum,na.rm=T)>0)
    x2 <- data.frame(eid=pheno1$eid,self_reported=x)
  } else if (disease_score=='eczema') {
    doctor_diag <- as.numeric(apply(pheno1[,c('6152-0.0','6152-0.1','6152-0.2','6152-0.3','6152-0.4')],1,function(x){9 %in% x}))
    x <- as.numeric(apply(data.frame(self_reported,doctor_diag),1,sum,na.rm=T)>0)
    x2 <- data.frame(eid=pheno1$eid,self_reported=x)
  } else {
    x2 <- data.frame(eid=pheno1$eid,self_reported)
  }
  df2 <- merge(df,x2,by='eid')
  df2$disease <- as.numeric(df2$case==1 | df2$oper==1 | df2$self_reported==1)
  
  field_codes <- c('21001-0.0','31-0.0','21022-0.0',paste0('22009-0.',1:10))
  field_names <- c('bmi','sex','age',paste0('PC',1:10))
  x <- pheno1[,c('eid',field_codes)]
  colnames(x)[-1] <- field_names
  df3 <- merge(df2,x,by='eid')
  
  library(dplyr)
  df4 <- as.data.frame(df3 %>% group_by(eid) %>% slice(which.min(days)))

  fwrite(df4,
         paste0('/athena/elementolab/scratch/anm2868/DrugPGS/UKB_pheno_data/',disease_score,'.txt'),
         quote = F,sep = '\t',na = 'NA',row.names = F,col.names = T)
  return(df4)
  
}
