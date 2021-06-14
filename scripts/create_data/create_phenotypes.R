library(data.table)
library(lubridate)

source('/athena/elementolab/scratch/anm2868/DrugPGS/DrugPGS_Interactions/scripts/read_in_phenotype_data.R')
source('/athena/elementolab/scratch/anm2868/DrugPGS/DrugPGS_Interactions/scripts/functions.R')

disease_score <- 'breast_cancer'
icd9_code <- c('174',paste0('174',0:9))
icd10_code <- c('C50',paste0('C50',0:9))
cancer_code <- TRUE; self_report_code <- '1002'
doctor_diag_avail <- FALSE; doctor_diag_code <- NA
female_only <- TRUE
opcs_code_avail <- TRUE; opcs_code <- c('B27|B28|B29')

df3 <- create_phenotype_file()
