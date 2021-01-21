# This script merges the original info files "ukb_mfi_chr${chr}_v3" into a single merged file.
# The single merged file has a column containing chromosome information.

# Merge into single file and append chr column
dir=/home/kulmsc/athena/ukbiobank/qc/imputed
mydir=/athena/elementolab/scratch/anm2868/DrugPGS/DrugPGS_Interactions/output/ss/ukb
outF=ukb_mfi_chrAll_v3
rm $mydir/$outF.txt
for chr in {1..22};
do
  print $chr
  f=ukb_mfi_chr${chr}_v3 # does not include CHR
  awk -v OFS='\t' -v chr=$chr '{print $0, chr}' $dir/$f.txt >> $mydir/$outF.txt
done

#Reduce the file to only contain SNPs with INFO score greater than 0.9
awk '$8 > 0.9 {print $0}' $mydir/$outF.txt > $mydir/impute_rsids

#Extract all of the duplicated RSIDs
cut -f2 $mydir/$outF.txt | sort | uniq -d > $mydir/dup_ids

#Remove the duplicated RSIDs from the full list of (INFO approved) SNPs
cat $mydir/impute_rsids | fgrep -w -v -f $mydir/dup_ids > $mydir/temp

# Final output: impute_rsids
mv $mydir/temp $mydir/impute_rsids

# Create a british-only eid file
cut -f1 /home/kulmsc/athena/doc_score/qc/fam_files/qc_brit_fam > ${myDir}/brit_eid

# make edits to $outF
f=$mydir/$outF.txt
extra_col=$mydir/extra_col
awk -v OFS=":" '{print $9,$3"_"}' $f > $extra_col
paste $f $extra_col > $f.tmp
mv $f.tmp $f