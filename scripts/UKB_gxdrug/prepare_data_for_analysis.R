prepare_data_for_analysis <- function(j) { 
  disease_name <- disease_lst[j]
  print(paste0('Disease ',j,'/',length(disease_lst),': ',disease_name))
  f <- paste0('/athena/elementolab/scratch/anm2868/DrugPGS/UKB_pheno_data/',disease_name,'.txt')
  df3 <- fread(f,data.table = F,stringsAsFactors = F)
  #########
  # A
  # library(dplyr)
  # df3 <- as.data.frame(df3 %>% group_by(eid) %>% slice(which.min(days)))
  # B
  # df3 <- df3[!duplicated(df3$eid),]
  #########
  df4 <- merge(df3,drug_data,by='eid')

  if (disease_name %in% c('breast_cancer','ovarian_cancer')) {
    dataf <- subset(df4,days > 0 & self_reported==0 & sex==0) # this finds all individuals not yet w/ disease
  } else if (disease_name %in% c('prostate_cancer','colorectal_cancer')) {
    dataf <- subset(df4,days > 0 & self_reported==0 & sex==1) # this finds all individuals not yet w/ disease
  } else {
    dataf <- subset(df4,days > 0 & self_reported==0) # this finds all individuals not yet w/ disease
  }
  x <- apply(dataf[,drug_names_full],2,sum)
  x <- data.frame(E=names(x),count=unname(x),stringsAsFactors = F)
  drug_names <- (subset(x,count>=1000))[,1]
  
  #######
  
  dataf.mg <- merge(dataf,polyscore[,c('eid',disease_name)],by='eid')
  if (use_percentile & remainder_popu) {
    dataf.mg[,disease_name] <- as.numeric(dataf.mg[,disease_name] > quantile(dataf.mg[,disease_name],probs = (1-percentile)))
  } else if (use_percentile & !remainder_popu) {
    polygenic_score <- dataf.mg[,disease_name]
    ind.top <- which(dataf.mg[,disease_name] > quantile(dataf.mg[,disease_name],probs=(1-percentile)))
    ind.bot <- which(dataf.mg[,disease_name] < quantile(dataf.mg[,disease_name],probs=percentile))
    dataf.mg[,disease_name] <- NA
    dataf.mg[ind.top,disease_name] <- 1
    dataf.mg[ind.bot,disease_name] <- 0
  }

  if (disease_name%in%c('breast_cancer','ovarian_cancer')) {
    tmp <- pheno2[,c('eid','2724-0.0')]
    tmp[,2] <- as.numeric(gsub(' ','',tmp[,2])); tmp[,2] <- as.numeric(tmp[,2]==1)
    # tmp[,3] <- as.numeric(gsub(' ','',tmp[,3])); tmp[,3] <- as.numeric(tmp[,3]==1); colnames(tmp) <- c('eid','menopause','oral_contraceptive')
    colnames(tmp) <- c('eid','menopause')
    dataf.mg <- merge(dataf.mg,tmp,by='eid')
    dataf.mg$menopause[is.na(dataf.mg$menopause)] <- 0
    
    tmp <- pheno2[,c('eid','2734-0.0')]
    colnames(tmp) <- c('eid','number_live_birth')
    dataf.mg <- merge(dataf.mg,tmp,by='eid')
    dataf.mg$number_live_birth <- as.numeric(gsub(' ','',dataf.mg$number_live_birth)); 
    dataf.mg$number_live_birth[dataf.mg$number_live_birth==-3] <- 0
    dataf.mg$number_live_birth[is.na(dataf.mg$number_live_birth)] <- 0
    dataf.mg$one_birth <- as.numeric(dataf.mg$number_live_birth>0)
    
  }
  
  dataf.mg$bmi[is.na(dataf.mg$bmi)] <- median(dataf.mg$bmi,na.rm = T)
  
  to_be_returned <- list(disease_name,drug_names,dataf.mg)
  return(to_be_returned)
}