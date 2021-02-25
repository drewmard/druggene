for chr in {1..22}; do
echo $chr
cut -f 2 /athena/elementolab/scratch/kulmsc/refs/1000genomes/eur.1.bim | sort | uniq -d > /athena/elementolab/scratch/anm2868/druggene/output/1000g/dups.$chr
done