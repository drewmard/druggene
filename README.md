# druggene

## Information about the repository

TBD

## Creating polygenic scores
Directory: scripts/create_score

#### 1. Create the list of valid SNPs in UKB.
impute_snp.sh

#### 2. QC the summary statistics file.
ss_preprocess.R

#### 3. Clump & create polygenic scores.
scoring_file.sh

### Creating polygenic scores from PGS catalog
Directory: scripts/create_score_PGS_catalog

	Note: use scripts/create_score/impute_snp.sh as a pre-processing first step

#### 2. QC the summary statistics file.
*Option 1: PGS_catalog_process.R:* 

#### 3. 
scoring_file.PGS_catalog.sh