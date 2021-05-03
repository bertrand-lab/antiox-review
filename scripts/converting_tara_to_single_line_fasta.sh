awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' < ../data/tara_ocean_smags/SMAGs_v1_concat.faa > ../data/tara_ocean_smags/SMAGs_v1_concat_single.fasta

