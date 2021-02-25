# Original script by Scott Kulm. Adapted by Andrew Marderstein.
# https://kulmsc.github.io/pgs_book/

source activate PGS

# PTHRES=0.0005
# R2THRES=0.1
# KBTHRES=250
# author=michailidou
author=schumacher
genoDir=/home/kulmsc/athena/ukbiobank/imputed
brit_eid=/athena/elementolab/scratch/anm2868/druggene/output/ukb/brit_eid

myDir=/athena/elementolab/scratch/anm2868/druggene/output/ss/${author}
mkdir -p $myDir

#######################################################################################
# 1:
# Split the summary statistics file into 22 files, 1 per chr:
mkdir -p ${myDir}/chr_ss
for chr in {1..22};do

	echo "Splitting chr $chr summary statistics..."
	
	# create file w/ header row
	zcat ${myDir}/clean_${author}.txt.gz | head -1 > ${myDir}/chr_ss/${author}_${chr}.ss
	
	# append remaining rows w/ right chr
    zcat ${myDir}/clean_${author}.txt.gz | awk -v var="$chr" '$1 == var {print $0}' >> ${myDir}/chr_ss/${author}_${chr}.ss
    
    # compress
    gzip -f ${myDir}/chr_ss/${author}_${chr}.ss
done

#######################################################################################
# 2:
# Subset SNPs & indivs for genotype file. Clumping will be applied after.
for chr in {10..22};do

echo "Subset genotype files for chr$chr..."

genoName=ukbb.${chr}
	
# pull out rsids from the ss file for a chr
zcat ${myDir}/chr_ss/${author}_${chr}.ss.gz > ${myDir}/ss
cut -f5 ${myDir}/ss > ${myDir}/rsids
	
# use bgenix to subset snps from large ukb imputed file 
/home/kulmsc/bin/bgenix -g ${genoDir}/${genoName}.bgen -incl-rsids ${myDir}/rsids > ${myDir}/temp.bgen
  	
# turn into a plink file, keeping only british european individuals
plink2 --memory 12000 --threads 12 --bgen ${myDir}/temp.bgen ref-first --sample ${genoDir}/${genoName}.sample --keep-fam $brit_eid --make-bed --out ${myDir}/geno.${chr}

done

#######################################################################################
# 3:
# Performing clumping:
mkdir -p ${myDir}/clumped
for chr in {1..22};do
# for chr in {1..22};do

zcat ${myDir}/chr_ss/${author}_${chr}.ss.gz > ${myDir}/ss

for PTHRES in 0.05 0.005 0.0005 0.00001 5e-8; do
# for PTHRES in 0.05 0.005; do
for R2THRES in 0.1; do
for KBTHRES in 250; do

echo "Clumping chr $chr summary statistics..."

mkdir -p ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}
# clump snps for scoring
# plink --bfile ${myDir}/geno.${chr} --clump ${myDir}/ss --clump-snp-field RSID --clump-field P --clump-p1 $PTHRES --clump-p2 $PTHRES --clump-r2 $R2THRES --clump-kb $KBTHRES --out ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/ss.${chr}.clump
plink --bfile /athena/elementolab/scratch/kulmsc/refs/1000genomes/eur.$chr --exclude /athena/elementolab/scratch/anm2868/druggene/output/1000g/dups.$chr --clump ${myDir}/ss --clump-snp-field RSID --clump-field P --clump-p1 $PTHRES --clump-p2 $PTHRES --clump-r2 $R2THRES --clump-kb $KBTHRES --out ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/ss.${chr}.clump

# extract snp names
awk '{ print $3 }' ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/ss.${chr}.clump.clumped > ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/rsid.${chr}
done
done
done
done


#######################################################################################
# 4:
# Compute polygenic scores:
# mkdir -p ${myDir}/score
for chr in {1..22};do
# for chr in {22..22};do

echo "Scoring chr $chr..."
zcat ${myDir}/chr_ss/${author}_${chr}.ss.gz > ${myDir}/ss
genoName=ukbb.${chr}

for PTHRES in 0.05 0.005 0.0005 0.00001 5e-8; do
for R2THRES in 0.1; do
for KBTHRES in 250; do

# for PTHRES in 0.0005 0.00001 5e-8; do
# for R2THRES in 0.1; do
# for KBTHRES in 250 5000; do


mkdir -p ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}

# # use bgenix to subset snps from large ukb imputed file 
/home/kulmsc/bin/bgenix -g ${genoDir}/ukbb.${chr}.bgen -incl-rsids ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/rsid.${chr} > ${myDir}/temp.bgen
  	
# # turn into a plink file, keeping only british european individuals
plink2 --memory 12000 --threads 12 --bgen ${myDir}/temp.bgen ref-first --sample ${genoDir}/ukbb.${chr}.sample --keep-fam $brit_eid --make-bed --out ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/geno.${chr}

# plink2 --memory 12000 --threads 12 --bfile ${myDir}/geno.${chr} --keep-fam $brit_eid --extract ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/rsid.${chr} --make-bed --out ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/geno.${chr}

# perform scoring
# are you using SE or not SE in the table?
# SE:
# plink --bfile ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/geno.${chr}.v2 --keep-allele-order --score ${myDir}/ss 5 3 7 sum --out ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/score_chr${chr}
# rm ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/geno.${chr}.v2.{bed,bim,fam} # waste of space
# No SE:
plink --bfile ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/geno.${chr} --keep-allele-order --score ${myDir}/ss 5 3 6 sum --out ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/score_chr${chr}
rm ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/geno.${chr}.{bed,bim,fam}
done
done
done
done

#######################################################################################
# 5:
# Combine polygenic scores on chromosome level into one file:

for PTHRES in 0.05 0.005 0.0005 0.00001 5e-8; do
for R2THRES in 0.1; do
for KBTHRES in 250; do

echo "PTHRES=${PTHRES}, R2THRES=${R2THRES}, KBTHRES=${KBTHRES}"
chr=1
echo $chr
awk -v OFS='\t' '{ print $2, $6 }' ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/score_chr${chr}.profile > ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/score_chrALL.profile
for chr in {2..22}; do
echo $chr
awk '{ print $6 }' ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/score_chr${chr}.profile > ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/tmp
paste ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/score_chrALL.profile ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/tmp > ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/tmpout && mv ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/tmpout ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/score_chrALL.profile
done
rm ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/tmp

done
done
done

rm ${myDir}/temp.bgen