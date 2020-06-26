
# making databases for G. rostochiensis and G. pallida using their cdna files in FASTA format


C:\Users\USER1>cd C:\

C:\>mkdir blast\db

C:\> cd blast\db\Ros

C:\blast\db\Ros>makeblastdb -in G_rostochiensis_L19_codingseq.fasta -dbtype nucl  

C:\blast\db\Ros>cd ..

C:\blast\db>cd Pal

C:\blast\db\Pal>makeblastdb -in G_pallida_D383_v.0.8.1.minus_scaff93A.cdna.fasta -dbtype nucl 




##reciprocal BLAST 

C:\blast\db\Pal>cd ..

C:\blast\db>blastn -db C:\blast\db\Pal\G_pallida_D383_v.0.8.1.minus_scaff93A.cdn
a.fasta -query C:\blast\db\Ros\G_rostochiensis_L19_codingseq.fasta -out rosto_on
_pal.blastn -outfmt "6 qseqid qacc sacc evalue bitscore length pident qstart qen
d sstart send" -word_size 10

C:\blast\db>blastn -db C:\blast\db\Ros\G_rostochiensis_L19_codingseq.fasta -out rosto_on
_pal.blastn  -query C:\blast\db\Pal\G_pallida_D383_v.0.8.1.minus_scaff93A.cdn
a.fasta -outfmt "6 qseqid qacc sacc evalue bitscore length pident qstart qen
d sstart send" -word_size 10

