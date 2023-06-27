#reference transcriptome and GTF for STAR
wget http://ftp.ensembl.org/pub/release-108/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz
gunzip -k Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz

#downloading vcf file and indexing
wget http://ftp.ensembl.org/pub/release-108/variation/vcf/saccharomyces_cerevisiae/saccharomyces_cerevisiae.vcf.gz
tabix saccharomyces_cerevisiae.vcf.


#indexing the genome using samtools
samtools faidx Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa.gz

#You can also do it using picard
#indexing the genome using picard
#java -Xmx4g -jar /home/danny/software/picard/build/libs/picard.jar CreateSequenceDictionary \
#     -R Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa.gz


#generating genome index using STAR
STAR --runThreadN 2 --runMode genomeGenerate \
     --genomeDir ~/genome/STAR \
     --genomeSAindexNbases 10 \
     --sjdbGTFfile ~/genome/Saccharomyces_cerevisiae.R64-1-1.108.gtf \
     --genomeFastaFiles ~/genome/Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa


