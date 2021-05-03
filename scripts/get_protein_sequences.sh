# this script goes through each gene list and then subsets the main fasta file to get the protein sequences for each antioxidant group

for i in ../data/antiox_gene_name_lists/*txt;
do
echo $i

    output_fasta=${i/.txt/}
    echo $output_fasta
    grep -w -A 1 -f $i ../data/tara_ocean_smags/SMAGs_v1_concat_single.fasta --no-group-separator > $output_fasta'.fasta'

done
