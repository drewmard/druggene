author=michailidou
genoDir=/home/kulmsc/athena/ukbiobank/imputed
brit_eid=/athena/elementolab/scratch/anm2868/druggene/output/ukb/brit_eid
myDir=/athena/elementolab/scratch/anm2868/druggene/output/ss/${author}

PTHRES=0.00001
R2THRES=0.1
KBTHRES=250


for chr in {1..22};do

echo "Scoring chr $chr..."
zcat ${myDir}/chr_ss/${author}_${chr}.ss.gz > ${myDir}/ss
genoName=ukbb.${chr}


# # use bgenix to subset snps from large ukb imputed file 
/home/kulmsc/bin/bgenix -g ${genoDir}/ukbb.${chr}.bgen -incl-rsids ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/rsid.${chr} > ${myDir}/temp.bgen
  	
# # turn into a plink file, keeping only british european individuals
plink2 --memory 12000 --threads 12 --bgen ${myDir}/temp.bgen ref-first --sample ${genoDir}/ukbb.${chr}.sample --keep-fam $brit_eid --make-bed --out ${myDir}/clumped/P_${PTHRES}.R2_${R2THRES}.KB_${KBTHRES}/geno.${chr}

done

################

author=PGS000007

genoDir=/home/kulmsc/athena/ukbiobank/imputed
brit_eid=/athena/elementolab/scratch/anm2868/druggene/output/ukb/brit_eid
myDir=/athena/elementolab/scratch/anm2868/druggene/output/ss/${author}
FLAG=""

zcat $myDir/clean_${author}${FLAG}.txt.gz > $myDir/ss

for chr in {1..22};do

echo "Scoring chr $chr..."
genoName=ukbb.${chr}

# use bgenix to subset snps from large ukb imputed file 
/home/kulmsc/bin/bgenix -g ${genoDir}/${genoName}.bgen -incl-rsids ${myDir}/rsids${FLAG} > ${myDir}/temp.bgen
  	
# turn into a plink file, keeping only british european individuals
plink2 --memory 12000 --threads 12 --bgen ${myDir}/temp.bgen ref-first --sample ${genoDir}/${genoName}.sample --keep-fam $brit_eid --make-bed --out ${myDir}/geno.${chr}

done

