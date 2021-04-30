# for uk biobank snps

SNPLIST=/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.SNP_names.txt
output=/athena/elementolab/scratch/anm2868/druggene/output/pgs_int/SNP_by_drug.SNP_names.query.txt

rm $output
while read -r SNP; do
if [[ "$SNP" == *":"* ]]; then
SNP=$(echo $SNP | sed 's/:/_/g')
SNP=${SNP}_b37
else
SNP=${SNP},
fi
query=$(grep -m 1 $SNP /athena/elementolab/scratch/anm2868/open_targets/csv_files/variant_index.csv)
echo $query >> $output
done < $SNPLIST

