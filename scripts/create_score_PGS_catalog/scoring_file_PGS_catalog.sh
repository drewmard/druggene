source activate PGS

# author=PGS000004; use_col="6 3 5" # hit A2
author=PGS000007; use_col="6 3 5"

genoDir=/home/kulmsc/athena/ukbiobank/imputed
brit_eid=/athena/elementolab/scratch/anm2868/druggene/output/ukb/brit_eid

myDir=/athena/elementolab/scratch/anm2868/druggene/output/ss/${author}

for FLAG in ".no_remove" ""; do

zcat $myDir/clean_${author}${FLAG}.txt.gz > $myDir/ss
# Compute polygenic scores:
for chr in {1..22};do
# for chr in {22..22};do

echo "Scoring chr $chr..."

genoName=ukbb.${chr}

# use bgenix to subset snps from large ukb imputed file 
/home/kulmsc/bin/bgenix -g ${genoDir}/${genoName}.bgen -incl-rsids ${myDir}/rsids${FLAG} > ${myDir}/temp.bgen
  	
# turn into a plink file, keeping only british european individuals
plink2 --memory 12000 --threads 12 --bgen ${myDir}/temp.bgen ref-first --sample ${genoDir}/${genoName}.sample --keep-fam $brit_eid --make-bed --out ${myDir}/geno.${chr}

# plink2 --memory 12000 --threads 12 --bfile ${myDir}/geno.${chr} --keep-fam $brit_eid --extract ${myDir}/rsids --make-bed --out ${myDir}/geno.${chr}.v2

# perform scoring
# plink --bfile ${myDir}/geno.${chr} --keep-allele-order --score ${myDir}/ss 7 4 5 sum --out --out ${myDir}/score_chr${chr}${FLAG} # PGS000007
plink --bfile ${myDir}/geno.${chr} --keep-allele-order --score ${myDir}/ss $use_col sum --out ${myDir}/score_chr${chr}${FLAG}
rm ${myDir}/geno.${chr}.{bed,bim,fam} # waste of space

done
done
rm ${myDir}/temp.bgen # waste of space


#######################################################################################
# 5:
# Combine polygenic scores on chromosome level into one file:

for FLAG in ".no_remove" ""; do

chr=1
awk -v OFS='\t' '{ print $2, $6 }' ${myDir}/score_chr${chr}${FLAG}.profile > ${myDir}/score_chrALL${FLAG}.profile
for chr in {2..22}; do
echo $chr
awk '{ print $6 }' ${myDir}/score_chr${chr}${FLAG}.profile > ${myDir}/tmp
paste ${myDir}/score_chrALL${FLAG}.profile ${myDir}/tmp > ${myDir}/tmpout && mv ${myDir}/tmpout ${myDir}/score_chrALL${FLAG}.profile
done
rm ${myDir}/tmp

done
