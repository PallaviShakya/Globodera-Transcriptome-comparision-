
## using HISAT2 for mapping 
for fwd_read_file in /mnt/nemahomes/shaky005/RNAseq_data/raw_data/*_1.fastq
do
	echo $fwd_read_file
    rev_read_file=/mnt/nemahomes/shaky005/RNAseq_data/raw_data/$(basename $fwd_read_file _1.fastq)_2.fastq
	echo $rev_read_file
    output_file=/mnt/nemahomes/shaky005/RNAseq_data/output_pal/$(basename $fwd_read_file _1.fastq).BAM
	echo $output_file
    hisat2 -p 8 --dta -x /mnt/nemahomes/shaky005/RNAseq_data/genome_index/G_pallida_D383 -1 $fwd_read_file -2 $rev_read_file | samtools sort -@ 2 -m 20G -o $output_file
done

for fwd_read_file in /mnt/nemahomes/shaky005/RNAseq_data/raw_data_rosto/*_1.fastq
do
	echo $fwd_read_file
    rev_read_file=/mnt/nemahomes/shaky005/RNAseq_data/raw_data_rosto/$(basename $fwd_read_file _1.fastq)_2.fastq
	echo $rev_read_file
    output_file=/mnt/nemahomes/shaky005/RNAseq_data/output_rosto/$(basename $fwd_read_file _1.fastq).BAM
	echo $output_file
    hisat2 -p 8 --dta -x /mnt/nemahomes/shaky005/RNAseq_data/genome_index/Grost19_08 -1 $fwd_read_file -2 $rev_read_file | samtools sort -@ 2 -m 20G -o $output_file
done

##stringtie for assembly 
##downloading stringtie

git clone https://github.com/gpertea/stringtie
cd stringtie 
make release 


for bam_files in /mnt/nemahomes/shaky005/RNAseq_data/output_pal/*.BAM
do 
	echo $bam_files
	output=/mnt/nemahomes/shaky005/RNAseq_data/stringtie_pal/$(basename $bam_files .BAM).gtf
	echo $output
    stringtie $bam_files  -P 8 -G /mnt/nemahomes/steen176/pallida/G_pallida_D383_v.0.8.1.structural_annotation_minus_scaffold_93A.renamed.gff -o $output -A /mnt/nemahomes/shaky005/RNAseq_data/stringtie_pal/pal_abundance.tab -B -e 
done


for bam_files in /mnt/nemahomes/shaky005/RNAseq_data/output_rosto/*.BAM
do
	echo $bam_files
	output=/mnt/nemahomes/shaky005/RNAseq_data/stringtie_rosto/$(basename $bam_files .BAM).gtf
	echo $output
	stringtie $bam_files  -P 8 -G /mnt/nemahomes/steen176/Annotation/G_ros19/braker_results/augustus.hints.gff3 -o $output -A /mnt/nemahomes/shaky005/RNAseq_data/stringtie_rosto/rosto_abundance.tab -B -e 
done

 python stringtie_merger.py --stringtie-dir /mnt/nemahomes/shaky005/RNAseq_data/stringtie_pal --outdir /mnt/nemahomes/shaky005/RNAseq_data/stringtie_merge_pal/ --suffix G.pallida --model /mnt/nemahomes/steen176/pallida/G_pallida_D383_v.0.8.1.structural_annotation_minus_scaffold_93A.renamed.gff
