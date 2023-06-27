#NOTE: Run this script inside "vaibhakv_rnaseqSATR/genome/" folder

url <- "https://ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/"
base <- "Saccharomyces_cerevisiae.R64-1-1.dna.chromosome."
chrs <- c(as.character(as.roman(seq(1:2))), "Mito")

#Download 
for(chr in chrs){
  fname <- paste0(base, chr, ".fa.gz")
  #Download
  cmd <- paste0("wget ", url, fname)
  system(cmd)
}

cat("", file = 	"Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa")
for(chr in chrs){
  fname <- paste0(base, chr, ".fa.gz")
  #Extract and merge into a file
  cmd <- paste0("zcat ", fname, " >> Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa")
  system(cmd)
}

#Compress the fasta file using bgzip
cmd <- paste0("bgzip -k Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa")
system(cmd)

#Delete the chromosomes
for(chr in chrs){
  fname <- paste0(base, chr, ".fa.gz")
  cmd <- paste0("rm ", fname)
  system(cmd)
}
