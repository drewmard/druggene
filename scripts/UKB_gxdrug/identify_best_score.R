library(data.table)
library(reshape2)
library(survival)
# author <- 'schumacher'

nagelkerke <-function(fit, null=NULL, restrictNobs=FALSE){
    TOGGLE =   (class(fit)[1]=="lm"
                | class(fit)[1]=="gls"
                | class(fit)[1]=="lme"
                | class(fit)[1]=="glm"
                | class(fit)[1]=="negbin"
                | class(fit)[1]=="zeroinfl"
                | class(fit)[1]=="clm"
                | class(fit)[1]=="vglm"
                | class(fit)[1]=="betareg"
                | class(fit)[1]=="rq")
    BOGGLE =   (class(fit)[1]=="nls"
                | class(fit)[1]=="lmerMod"
                | class(fit)[1]=="glmerMod"
                | class(fit)[1]=="merModLmerTest"
                | class(fit)[1]=="lmerModLmerTest"
                | class(fit)[1]=="clmm")
    SMOGGLE =   (class(fit)[1]=="lmerMod"
                 | class(fit)[1]=="glmerMod"
                 | class(fit)[1]=="merModLmerTest"
                 | class(fit)[1]=="lmerModLmerTest"
                 | class(fit)[1]=="vglm")
    ZOGGLE  = (class(fit)[1]=="zeroinfl")
    ZOGGLE2 = (class(fit)[1]=="rq")
    NOGGLE = is.null(null)
    ERROR  = "Note: For models fit with REML, these statistics are based on refitting with ML"
    ERROR2 = "None"
    
    if(!restrictNobs & NOGGLE  & TOGGLE){null = update(fit, ~ 1)}
    if(restrictNobs  & NOGGLE  & TOGGLE){null = update(fit, ~ 1, data=fit$model)}
    
    if(restrictNobs  & !NOGGLE){null = update(null, data=fit$model)}
    
    if(NOGGLE & BOGGLE)
    {ERROR = "You need to supply a null model for nls, lmer, glmer, or clmm"}
    if((!TOGGLE) & (!BOGGLE))
    {ERROR = "This function will work with lm, gls, lme, lmer, glmer, glm, negbin, zeroinfl, nls, clm, clmm, and vglm"}
    
    SMOGGLE2 = (class(null)[1]=="lmerMod"
                | class(null)[1]=="glmerMod"
                | class(null)[1]=="merModLmerTest"
                | class(null)[1]=="lmerModLmerTest"
                | class(null)[1]=="vglm")   
    
    Y = matrix(rep(NA,2),
               ncol=1)
    colnames(Y) = ""
    rownames(Y) = c("Model:", "Null:")
    
    Z = matrix(rep(NA, 3),
               ncol=1)
    colnames(Z) = c("Pseudo.R.squared")
    rownames(Z) = c("McFadden", "Cox and Snell (ML)", 
                    "Nagelkerke (Cragg and Uhler)")
    
    X = matrix(rep(NA,4),
               ncol=4)
    colnames(X) = c("Df.diff","LogLik.diff","Chisq","p.value")
    rownames(X) = ""
    
    U = matrix(rep(NA,2),
               ncol=1)
    colnames(U) = ""
    rownames(U) = c("Model:", "Null:")
    
    if(TOGGLE | BOGGLE){
      if (!SMOGGLE){Y[1]= toString(fit$call)}
      if (SMOGGLE){Y[1]= toString(fit@call)}
    }
    
    if(TOGGLE | (BOGGLE & !NOGGLE)){
      
      if (!SMOGGLE2){Y[2]= toString(null$call)}
      if (SMOGGLE2){Y[2]= toString(null@call)}
      
      if(!ZOGGLE & !ZOGGLE2){N = nobs(fit)
      U[1,1]= nobs(fit); U[2,1]= nobs(null)}
      if(!ZOGGLE &  ZOGGLE2){N = length(fit$y)
      U[1,1]= length(fit$y); U[2,1]= length(null$y)}
      if(ZOGGLE){N = fit$n
      U[1,1]= fit$n; U[2,1]= null$n}
      
      if(U[1,1] != U[2,1]){
        ERROR2 = "WARNING: Fitted and null models have different numbers of observations"}
      
      m = suppressWarnings(logLik(fit, REML=FALSE))[1]
      n = suppressWarnings(logLik(null, REML=FALSE))[1]
      mf = 1 - m/n
      Z[1,] = signif(mf, digits=6)
      cs = 1 - exp(-2/N * (m - n))
      Z[2,] = signif(cs, digits=6)
      nk = cs/(1 - exp(2/N * n))
      Z[3,] = signif(nk, digits=6)
      
      o = n - m
      dfm = attr(logLik(fit),"df")
      dfn = attr(logLik(null),"df")
      if(class(fit)[1]=="vglm"){dfm=df.residual(fit)}
      if(class(fit)[1]=="vglm"){dfn=df.residual(null)}
      dff = dfn - dfm
      CHI = 2 * (m - n)
      P = pchisq(CHI, abs(dff), lower.tail = FALSE)
      
      X [1,1] = dff
      X [1,2] = signif(o, digits=5)             
      X [1,3] = signif(CHI, digits=5)
      X [1,4] = signif(P, digits=5)     
    }
    
    W=ERROR
    
    WW=ERROR2
    
    V = list(Y, Z, X, U, W, WW) 
    names(V) = c("Models", "Pseudo.R.squared.for.model.vs.null", 
                 "Likelihood.ratio.test", "Number.of.observations",
                 "Messages", "Warnings")
    return(V)            
  }

# Internal score processing:
author <- 'michailidou'
origdir <- '/athena/elementolab/scratch/anm2868/druggene'
iter <- 0
polyscore <- list()
for (PTHRES in c("0.05","0.005","0.0005","0.00001","5e-8")) {
  for (R2THRES in c("0.1")) {
    for (KBTHRES in c("250")) {
      iter <- iter + 1
      type <- paste0('P_',PTHRES,'.R2_',R2THRES,'.KB_',KBTHRES)
      print(paste0(iter,': ',type))
      dir <- paste0(origdir,'/output/ss/',author,'/clumped/',type)
      polyscore2 <- fread(paste0(dir,'/score_chrALL.profile'),data.table=F,stringsAsFactors = F)
      colnames(polyscore2)[2:23] <- paste0('score_',1:22)
      polyscore2[,2:ncol(polyscore2)] <- polyscore2[,2:ncol(polyscore2)]
      polyscore2$score_all <- apply(polyscore2[,-1],1,sum)
      polyscore2$score_all <- polyscore2$score_all
      polyscore2$score_all <- scale(polyscore2$score_all)[,1]
      colnames(polyscore2)[1] <- 'eid'; polyscore2[,type] <- polyscore2$score_all # this is necessary to work w/ previous script
      polyscore[[iter]] <- polyscore2[,c('eid',type)]
    }
  }
}

# External score processing:
iter <- 0
polyscore2.save <- list()
for (type in c('PGS000004','PGS000007')) {
  # for (FLAG in c('.no_remove','')) {
  for (FLAG in c('')) {
    iter <- iter + 1
    print(paste0(iter,': ',type,FLAG))
    f <- paste0(origdir,'/output/ss/',type,'/score_chrALL',FLAG,'.profile')
    polyscore2 <- fread(f,data.table=F,stringsAsFactors = F)
    colnames(polyscore2)[2:ncol(polyscore2)] <- paste0('score_',1:c(ncol(polyscore2)-1))
    polyscore2$score_all <- apply(polyscore2[,-1],1,sum)
    polyscore2$score_all <- polyscore2$score_all*-1
    polyscore2$score_all <- scale(polyscore2$score_all)[,1]
    colnames(polyscore2)[1] <- 'eid'; polyscore2[,paste0(type,FLAG)] <- polyscore2$score_all # this is necessary to work w/ previous script
    polyscore2.save[[iter]] <- polyscore2[,c('eid',paste0(type,FLAG))]
  }
}

polyscore <- do.call(cbind,polyscore)
polyscore <- polyscore[,!duplicated(colnames(polyscore))]
polyscore2 <- do.call(cbind,polyscore2.save)
polyscore2 <- polyscore2[,!duplicated(colnames(polyscore2))]
polyscore.full <- merge(polyscore,polyscore2,by='eid')

# 
polyscore <- polyscore.full[,c(1:2)]
colnames(polyscore)[2] <- 'breast_cancer'

# load
# pheno2 <- fread('/home/kulmsc/athena/ukbiobank/phenotypes/ukb33822.csv.gz',data.table=F,stringsAsFactors = F)
# level=4; drug_data <- fread(f <- paste0('/athena/elementolab/scratch/anm2868/DrugPGS/UKB_drug_data/level',level,'.txt'),data.table = F,stringsAsFactors = F)
# drug_names_full <- colnames(drug_data)[which(colnames(drug_data) != 'eid')]

# create final dataset:
disease_lst <- c('breast_cancer');j = 1
use_percentile <- FALSE; remainder_popu=TRUE; percentile <- .2
source("/athena/elementolab/scratch/anm2868/druggene/scripts/UKB_gxdrug/prepare_data_for_analysis.R")
returned <- prepare_data_for_analysis(j)
disease_name <- returned[[1]]
drug_names <- returned[[2]]
dataf.mg <- returned[[3]]
dataf.mg <- merge(dataf.mg,polyscore.full,by='eid')

## What is the scores to be tested? print to type.lst
type.lst <- c()
for (PTHRES in c("0.05","0.005","0.0005","0.00001","5e-8")) {
  for (R2THRES in c("0.1")) {
    for (KBTHRES in c("250")) {
      type <- paste0('P_',PTHRES,'.R2_',R2THRES,'.KB_',KBTHRES)
      type.lst <- c(type.lst,type)
    }
  }
}
for (type in c('PGS000004','PGS000007')) {
  for (FLAG in c('')) {
    type2 <- paste0(type,FLAG)
    type.lst <- c(type.lst,type2)
  }
}

# scale and standardize
for (score in type.lst) {
  dataf.mg[,score] <- scale(dataf.mg[,score])
}

# calculate Nagelkerke R2 for each score
disease_name <- disease_lst[j]
param.vec <- c()
null.mod <- glm(disease ~ bmi+age+menopause+
                  number_live_birth+one_birth+
                  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = dataf.mg,family = binomial(link='logit'))
for (i in 1:length(type.lst)) {
  print(paste0(i,'/',length(type.lst)))
  mod <- glm(disease ~ dataf.mg[,type.lst[i]]+bmi+age+menopause+
                 number_live_birth+one_birth+
                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = dataf.mg,family = binomial(link='logit'))
  param <- nagelkerke(fit=mod,null=null.mod)$Pseudo.R.squared.for.model.vs.null[3,]
  param.vec <- c(param.vec,param)
}
tmp <- data.frame(pgs_name=type.lst,r2=param.vec); tmp
fwrite(tmp,'/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/pgs_nagelkerke_r2.txt',quote = F,na = 'NA',sep = '\t',col.names = T,row.names = F)

# calculate internal vs external
dataf.mg$Internal.Top <- as.numeric(dataf.mg[,"P_0.00001.R2_0.1.KB_250"]>quantile(dataf.mg[,"P_0.00001.R2_0.1.KB_250"],probs = 0.9))
dataf.mg$External.Top <- as.numeric(dataf.mg[,"PGS000007"]>quantile(dataf.mg[,"PGS000007"],probs = 0.9))
# dataf.mg$External.Top <- as.numeric(dataf.mg[,"P_0.0005.R2_0.1.KB_250"]>quantile(dataf.mg[,"P_0.0005.R2_0.1.KB_250"],probs = 0.9))
table(dataf.mg$Internal.Top,dataf.mg$External.Top)
sum(dataf.mg$Internal.Top==0 & dataf.mg$External.Top==1)
sum(dataf.mg$Internal.Top==1 & dataf.mg$External.Top==0)

summary(mod)$coef

# create combined score
dataf.mg$pgs <- scale(apply(dataf.mg[,c('P_0.00001.R2_0.1.KB_250','PGS000007')],1,mean))


# save dataframe:
dataf.save <- dataf.mg[,c('eid','disease','days','pgs','S01BA','bmi','age','menopause',
                          'number_live_birth','one_birth',paste0('PC',1:10),
                          type.lst)]
fwrite(dataf.save,'/athena/elementolab/scratch/anm2868/druggene/output/pgs_drug.txt',quote = F,na = 'NA',sep = '\t',col.names = T,row.names = F)

# for interaction analysis:
dataf.save <- dataf.mg[,c('eid','disease','days','pgs','bmi','age','menopause',
                          'number_live_birth','one_birth',paste0('PC',1:10),
                          type.lst,drug_names)]
fwrite(dataf.save,'/athena/elementolab/scratch/anm2868/druggene/output/pgs_drug_names.txt',quote = F,na = 'NA',sep = '\t',col.names = T,row.names = F)

