for gtf_file in /mnt/nemahomes/shaky005/RNAseq_data/stringtie_pal/*.gtf
do
echo "$(basename $gtf_file .gtf)
/t$gtf_file" >> list_name.txt
done



for gtf_file in /mnt/nemahomes/shaky005/RNAseq_data/stringtie_rosto/*.gtf
do
echo "$(basename $gtf_file .gtf)
/t$gtf_file" >> list_name_rosto.txt
done