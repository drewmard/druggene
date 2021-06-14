# this pulls ICD9/10 codes for diagnosis
hesin <- fread('/home/kulmsc/athena/ukbiobank/hesin/hesin.txt',data.table = F,stringsAsFactors = F)
hesin_diag <- fread('/home/kulmsc/athena/ukbiobank/hesin/hesin_diag.txt',data.table = F,stringsAsFactors = F)
hesin_diag.sub <- hesin_diag[grep('J|L|M',hesin_diag$diag_icd10),]
tmp <- aggregate(hesin_diag.sub$ins_index,list(eid=hesin_diag.sub$eid),min)
hesin_tmp <- merge(hesin[,c('eid','ins_index','epistart','admidate')],tmp,by.x=c('eid','ins_index'),by.y=c('eid','x'))
hesin_tmp$epistart[is.na(hesin_tmp$epistart)] <- hesin_tmp$admidate[is.na(hesin_tmp$epistart)]
hesin_tmp$case <- 1
pheno_tmp <- pheno1[,c("eid","40000-0.0","53-0.0",'26410-0.0','26426-0.0','26427-0.0','191-0.0')]
pheno_tmp2 <- merge(pheno_tmp,hesin_tmp,by="eid",all.x=TRUE)
pheno_tmp2$case[is.na(pheno_tmp2$case)] <- 0

pheno_tmp2$opdate <- NA
pheno_tmp2$oper <- 0

# censor end date dependent on hospital/country
censor_locat <- which(!is.na(pheno_tmp2[,c('26410-0.0','26426-0.0','26427-0.0')]),arr.ind = TRUE)
pheno_tmp2$Hosp <- NA
pheno_tmp2$Hosp[censor_locat[censor_locat[,2]==1,1]] <- 'HES'
pheno_tmp2$Hosp[censor_locat[censor_locat[,2]==2,1]] <- 'PEDW'
pheno_tmp2$Hosp[censor_locat[censor_locat[,2]==3,1]] <- 'SMR'
# if hospital/country unknown, then use hospital codes to figure it out
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

df <- pheno_tmp2[,c('eid','53-0.0','case','oper','diagnosis','40000-0.0','191-0.0','End','opdate')]
colnames(df) <- c('eid','start','case','oper','diagnosis','death','lost','end','opdate')
df$death[df$death==''] <- NA
df$lost[df$lost==''] <- NA

df2 <- df
df2$self_reported <- 0
df2$disease <- as.numeric(df2$case==1 | df2$oper==1 | df2$self_reported==1)

field_codes <- c('21001-0.0','31-0.0','21022-0.0',paste0('22009-0.',1:10))
field_names <- c('bmi','sex','age',paste0('PC',1:10))
x <- pheno1[,c('eid',field_codes)]
colnames(x)[-1] <- field_names
df3 <- merge(df2,x,by='eid')

df3$first_diagnosis <- apply(df3[,c('diagnosis','death','lost','end','opdate')],1,min,na.rm=T)
i <- which(df3$self_reported==1)
df3$first_diagnosis_sr <- df3$first_diagnosis
df3$first_diagnosis_sr[i] <- apply(data.frame(start=df3$start[i],first_diagnosis=df3$first_diagnosis[i]),1,min)
df3$first_diagnosis_days <- apply(df3[,c('first_diagnosis_sr','start')],1,function(x) as.numeric(difftime(as.Date(x[1]),as.Date(x[2]))))

fwrite(df3,'/athena/elementolab/scratch/anm2868/DrugPGS/UKB_pheno_data/corticoster_indication.txt',quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

