PGS_name=PGS000007
# PGS_name=PGS000007
# PGS_name=PGS000008
# PGS_name=PGS000009
# PGS_name=PGS000082

# conda activate PGS

zcat $myDir/clean_${PGS_name}.txt.gz > $myDir/ss
for chr in {1..22};do

myDir=/athena/elementolab/scratch/anm2868/DrugPGS/DrugPGS_Interactions/output/PGS/PGS_Catalog/${PGS_name}
ss=myDir/clean_${PGS_name}.txt.gz
genoDir=/home/kulmsc/athena/ukbiobank/imputed
genoName=ukbb.${chr}
brit_eid=/athena/elementolab/scratch/anm2868/DrugPGS/DrugPGS_Interactions/output/PGS/michailidou/brit_eid

# use bgenix to subset snps from large ukb imputed file 
/home/kulmsc/bin/bgenix -g ${genoDir}/${genoName}.bgen -incl-rsids ${myDir}/rsids > ${myDir}/temp.bgen
  	
# turn into a plink file, keeping only british european individuals
plink2 --memory 12000 --threads 12 --bgen ${myDir}/temp.bgen ref-first --sample ${genoDir}/${genoName}.sample --keep-fam $brit_eid --make-bed --out ${myDir}/geno.${chr}

# perform scoring
ss=$myDir/ss
plink --bfile ${myDir}/geno.${chr} --keep-allele-order --score ${myDir}/ss 6 3 5 sum --out ${myDir}/score_chr${chr}
# rm ${myDir}/geno.${chr}.{bed,bim,fam} # waste of space
rm ${myDir}/temp.bgen

done

rm ${myDir}/geno*log
rm ${myDir}/score*{log,nopred,nosex}

chr=1
awk -v OFS='\t' '{ print $2, $6 }' ${myDir}/score_chr${chr}.profile > ${myDir}/score_chrALL.profile
for chr in {2..22}; do
echo $chr
awk '{ print $6 }' ${myDir}/score_chr${chr}.profile > ${myDir}/tmp
paste ${myDir}/score_chrALL.profile ${myDir}/tmp > ${myDir}/tmpout && mv ${myDir}/tmpout ${myDir}/score_chrALL.profile
rm ${myDir}/tmp
done

##################################

mkdir -p subset

snp=rs1058402
chr=19

myDir=/athena/elementolab/scratch/anm2868/DrugPGS/DrugPGS_Interactions/output/PGS/PGS_Catalog/subset
genoDir=/home/kulmsc/athena/ukbiobank/imputed
genoName=ukbb.${chr}
brit_eid=/athena/elementolab/scratch/anm2868/DrugPGS/DrugPGS_Interactions/output/PGS/michailidou/brit_eid

/home/kulmsc/bin/bgenix -g ${genoDir}/${genoName}.bgen -incl-rsids $snp > ${myDir}/temp.bgen
  	
# turn into a plink file, keeping only british european individuals
plink2 --memory 12000 --threads 12 --bgen ${myDir}/temp.bgen ref-first --sample ${genoDir}/${genoName}.sample --keep-fam $brit_eid --make-bed --out ${myDir}/geno.${chr}.${snp}
rm ${myDir}/temp.bgen

	