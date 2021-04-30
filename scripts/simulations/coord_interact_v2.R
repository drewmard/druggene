library(parallel)

# initialize:
# nsnp <- 200 # number of snps
# nindiv <- 100000 # number of individuals
# add_herit=0.4 # additive heritability (percent)
# int_herit=0.8 # interaction heritability (percent)
nsnp <- 5000 # number of snps
nindiv <- 5000 # number of individuals
add_herit=0.6 # additive heritability (percent)
int_herit=0.1 # interaction heritability (percent)
ncores=8 # number of computer cores to run... scale down if you have less!!
uncoord=FALSE; beta_cor=0 # uncoordinated interactions? (TRUE) or coordinated? (FALSE)
single_regulator=TRUE; MAF_regulator <- 0.4 # if a single factor regulating/interacting genome-wide
e_factor=TRUE; e_freq <- 0.4; e_levels <- 1
scaled_genotypes=TRUE # scale genotypes
prop_int = 0.2 # proportion of possible interactions that have non zero effect
# case_control=TRUE

# if random uncoordinated interactions:
# if (!single_regulator & uncoord) {prop_int <- 0.01} # proportion of possible interactions that have non zero effect

###################################################
set.seed(3)
###################################################

# https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables/15035#15035
complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

# generate MAF for SNPs
if (single_regulator & !e_factor) {
  MAF <- as.list(c(MAF_regulator,runif(nsnp-1,0.01,0.5))) 
} else {
  MAF <- as.list(runif(nsnp,0.01,0.5)) # generate MAF for SNPs
}

create_population <- function(pgs_beta=NULL,BETA=NULL,ind=NULL,BETA.ADD=NULL) {
  
  # simulate genotypes
  if (scaled_genotypes) {
    G <- do.call(cbind,mclapply(MAF,function(x) scale(rbinom(nindiv,2,x)),mc.cores = ncores)) # scaled
  } else {
    G <- do.call(cbind,mclapply(MAF,function(x) (rbinom(nindiv,2,x)),mc.cores = ncores))
  }
  
  # simulate environmental factor
  if (e_factor) {
    E <- scale(rbinom(nindiv,e_levels,e_freq))[,1]
  }
  
  ########################
  # If creating the GWAS dataset, then need to generate the genotype-phenotype mapping:
  
  # Additive effects
  if (is.null(BETA.ADD)) {
    BETA.ADD <- matrix(rnorm(nsnp),nsnp,1)
  }
  
  # which((G[,1:num_snp_interact]*E) == sapply(1:num_snp_interact,function(j) {G[,j]*E}))
  
  # Interaction matrix
  if (single_regulator & e_factor) { # single SNP interacts pairwise with other genome wide SNPs
    tmp <- E*G[,1:(nsnp*prop_int)]
  } else if (single_regulator & !e_factor) {
    tmp <- G[,1]*G[,2:(nsnp*prop_int)]
  } else {
    # design matrix for pairwise interactions
    # tmp <- model.matrix( ~.^2, data=as.data.frame(G))[,-c(1:(nsnp+1))]
    # this step identifies which pairwise interactions have a non-zero effect
    print("ERROR: NOT IMPLEMENTED"); break
    tmp <- sapply(1:(nsnp*(nsnp-1)*propInt),function(j) {
      int_pairs <- sample(nsnp,2,replace = F)
      G[,int_pairs[1]]*G[,int_pairs[2]]
    })
    ## ERROR: NOT COMPLETE............
  }
  
  # Interaction effects
  if (is.null(BETA)) {
    if (uncoord) {
      # uncoordinated: random interaction effects
      BETA <- matrix(rnorm(ncol(tmp)),ncol(tmp),1)
    } else {
      if (beta_cor==1) {
        if (!single_regulator) {
          print("ERROR: NOT IMPLEMENTED"); break
          # single SNP interacts pairwise with other genome wide SNPs
        } else if (single_regulator & !e_factor) {
          # coordinated interaction effects are proportional to additive effects
          BETA <- BETA.ADD[-1][1:(nsnp*prop_int)]
        } else if (single_regulator & e_factor) {
          BETA <- BETA.ADD[1:(nsnp*prop_int)]
        }
      } else if (single_regulator & e_factor) {
        BETA <- matrix(complement(BETA.ADD[1:(nsnp*prop_int)],beta_cor,rnorm(ncol(tmp))),ncol(tmp),1)
      }
    }
  }
  ########################
  # this step generates the genetic components of the final phenotype value
  Y.add <- scale(as.numeric(G %*% BETA.ADD))
  Y.int <- scale(as.numeric(tmp %*% BETA))
  
  ########################
  # environmental noise
  Y.env <- scale(rnorm(nindiv))
  
  # this step generates the final phenotype value
  # additive component is scaled to explain A% phenotype variance
  # interaction component is scaled to explain B% phenotype variance
  # and environmental component is scaled to explain 1-A-B% variance
  P <- as.numeric(scale(
    Y.add*sqrt(add_herit) +
      Y.int*sqrt(int_herit) +
      Y.env*sqrt(1-add_herit-int_herit)
  ))
  
  # label genotype matrix w/ snp names snp1,...,snpX
  colnames(G) <- paste0('snp',1:ncol(G))
  
  # dataframe to be returned
  df <- data.frame(G,P)
  if (e_factor) {df$E <- E}
  
  # calculate polygenic scores for a testing data set
  if (!is.null(pgs_beta)) {
    df$pgs <- scale(as.vector(G %*% matrix(pgs_beta,nsnp,1)))
    df$odd_pgs <- scale(as.vector(G[,seq(1,nsnp,by=2)] %*% matrix(pgs_beta[seq(1,nsnp,by=2)],nsnp/2,1)))
    df$even_pgs <- scale(as.vector(G[,seq(2,nsnp,by=2)] %*% matrix(pgs_beta[seq(2,nsnp,by=2)],nsnp/2,1)))
    if (single_regulator) {
      df$pgs.1 <- scale(as.vector(G[,-1] %*% matrix(pgs_beta[-1],nsnp-1,1)))
    }
  }
  return(list(df,BETA,BETA.ADD,Y.int,Y.add))
}

# calculated betas for creating polygenic scores
genetic_association <- function(i,logistic=TRUE) {
  mod <- lm(P~df[,paste0('snp',i)],data = df)
  res <- summary(mod)$coef[2,c(1:4)]
  return(res)
}

# create GWAS population
# for (e_levels in c(2)) {
# print(e_levels)

dataf <- list(); dataf.tmp <- c(); indexi <- 0
nsim=100;
for (beta_cor in seq(0,0.5,by=0.1)) {
  indexi <- indexi+1
  for (k in 1:nsim) {
    print(k)
    
    print('Creating original population...')
    create_population_results <- create_population()
    df <- create_population_results[[1]]
    BETA <- create_population_results[[2]]
    BETA.ADD <- create_population_results[[3]]
    Y.int <- create_population_results[[4]]
    Y.add <- create_population_results[[5]]
    
    # perform GWAS
    print('Performing GWAS...')
    res.add <- mclapply(1:nsnp,genetic_association,logistic=logistic,mc.cores = ncores)
    res.add <- as.data.frame(do.call(rbind,res.add))
    colnames(res.add) <- c('BETA.MEAN','SE.MEAN','T.MEAN','P.MEAN')
    
    # save data set
    print('Creating testing population...')
    # df <- create_population(pgs_beta=res.add$BETA.MEAN,BETA=BETA,ind=ind)[[1]]
    df <- create_population(pgs_beta=res.add$BETA.MEAN,BETA=BETA,ind=ind,BETA.ADD=BETA.ADD)[[1]]
    mse <- mean((res.add$BETA.MEAN - BETA.ADD)^2)
    
    # perform even-odd test
    # mod <- lm(P~odd_pgs*even_pgs,data=df)
    # summary(mod)$coef
    # 
    # mod <- lm(P~pgs.1*snp1,data=df)
    # summary(mod)$coef
    # 
    # if (e_factor) {
    #   mod <- lm(P~pgs*E,data=df)
    #   summary(mod)$coef
    # }
    
    # mod <- lm(P~pgs.1*snp1,data=df)
    # p <- c(p,summary(mod)$coef["pgs.1:snp1","Pr(>|t|)"])
    
    mod <- lm(P~pgs*E,data=df)
    cor.res <- cor.test(BETA.ADD,BETA[1:length(BETA.ADD)])
    dataf[[k]] <- data.frame(nindiv,nsnp,prop_int,e_freq=e_freq,e_levels=e_levels,uncoord,beta_cor,add_herit,int_herit,int_beta=summary(mod)$coef["pgs:E","Estimate"],int_p=summary(mod)$coef["pgs:E","Pr(>|t|)"],r=as.numeric(cor.res$estimate),r_pval=as.numeric(cor.res$p.value),mse)
  }
  dataf.tmp[[indexi]] <- do.call(rbind,dataf)
}
dataf.save <- do.call(rbind,dataf.tmp)
# dataf.save$beta_cor <- rep(seq(0,1,by=0.1),each=nsim)
aggregate(dataf.save$int_p,list(dataf.save$beta_cor),function(x) {mean(x<.05)})

##########3
print(mean(dataf.save$int_p<0.05))

# plot(dataf.save$r,dataf.save$int_beta)
cor.test(dataf.save$r,dataf.save$int_beta)
cor.test(dataf.save$mse,dataf.save$int_beta)
cor.test(dataf.save$mse,dataf.save$r)

aggregate(dataf.save,list(dataf.save$int_herit),function(x) {mean(x<.05)})
aggregate(dataf.save,list(dataf.save$beta_cor),mean)
aggregate(dataf.save,list(dataf.save$nindiv),median)
library(data.table)
fwrite(dataf.save,"/athena/elementolab/scratch/anm2868/druggene/output/sim/results3.txt",quote = F,na = 'NA',sep = '\t',col.names = T,row.names = F)
# cor.test(subset(dataf.save,nindiv==5000)$int_beta,subset(dataf.save,nindiv==5000)$r)
cor.test(subset(dataf.save,nindiv==50000)$int_beta,subset(dataf.save,nindiv==50000)$r)

# dataf.save
# cor.test(-log10(dataf.save$int_p),-log10(dataf.save$r_pval))
# mean(subset(dataf.save,r_pval>0.05)$int_p<0.05)


# }
