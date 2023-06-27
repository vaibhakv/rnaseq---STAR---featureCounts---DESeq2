#NOTE: Run this script inside "vaibhakv_rnaseqSATR/genome/" folder

#Aligning SRA reads to the Saccharomyces Cerevisiae genome


cmdlineargs -> commandArgs(trailingOnly = TRUE)

execute <- function(x, outputfile = NA, intern = FALSE, quitOnError = FALSE){
  if(!is.na(outputfile) && file.exists(outputfile)){
    cat("Output for step exists, skipping this step\n");
    invisible("")
  }
  cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
  if(res[1] >= 1){ 
    cat("Error external process did not finish\n\n");
    if(quitOnError) q("no")
  }
}

input.dir <- "~/vaibhakv_rnaseqSTAR/data/raw"
input.base <- cmdlineargs[1] #Previously SRR13978643 
output.dir <- paste0("~/vaibhakv_rnaseqSTAR/data/output/", input.base,".aln")
genome.path <- "~/vaibhakv_rnaseqSTAR/genome/STAR"
ref.fa.gz <- "~/vaibhakv_rnaseqSTAR/genome/Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa.gz"
ref.fa <- "~/vaibhakv_rnaseqSTAR/genome/Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa"
ref.snps <- "~/vaibhakv_rnaseqSTAR/genome/saccharomyces_cerevisiae.vcf.gz"

#Creating an output folder
if(!file.exists(input.dir)){ dir.create(input.dir, recursive = TRUE) }
if(!file.exists(output.dir)){ dir.create(output.dir, recursive = TRUE) }

#DOwnloading SRA files and Compressing
setwd(input.dir)

execute(paste0("fasterq-dump -p -e 4 --split-files ", input.base), paste0(input.base, "_1.fastq"))
execute(paste0("bgzip ", input.base, "_1.fastq"), paste0(input.base, "_1.fastq.gz"))
execute(paste0("bgzip ", input.base, "_2.fastq"), paste0(input.base, "_2.fastq.gz"))

#REad Trimming
trim.files  <- c(
                  paste0(input.dir, "/", input.base,"_1.fastq.gz"),
                  paste0(input.dir, "/", input.base,"_2.fastq.gz"),
                  paste0(output.dir, "/", input.base,"_1.P.fastq.gz"),
                  paste0(output.dir, "/", input.base,"_1.U.fastq.gz"),
                  paste0(output.dir, "/", input.base,"_2.P.fastq.gz"),
                  paste0(output.dir, "/", input.base,"_2.U.fastq.gz")
                )
trim.path <- "/usr/share/java/"
trim.exec <- paste0("java -jar ", trim.path, "/dist/jar/trimmomatic-0.40-rc1.jar")
trim.opts <- paste0("ILLUMINACLIP:",trim.path,"/adapters/TruSeq3-PE-2.fa:2:30:10")
trim.opts <- paste0(trim.opts, " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
trim.cmd  <- paste0(trim.exec, " PE ", paste0(trim.files, collapse=" "), " ", trim.opts)

execute(trim.cmd, trim.files[3])

execute(paste0("gunzip -k ", trim.files[3]), gsub(".fastq.gz", ".fastq", trim.files[3]))
execute(paste0("gunzip -k ", trim.files[5]), gsub(".fastq.gz", ".fastq", trim.files[5]))

files.in <- gsub(".fastq.gz", ".fastq", trim.files[c(3,5)])

#Alignment using STAR
star.outbase <- paste0(output.dir, "/", input.base)
star.bam <- paste0(star.outbase, "Aligned.sortedByCoord.out.bam")

star.exec <- "STAR --runMode alignReads"
star.opts <- paste0("--genomeDir=", genome.path, " --outSAMtype BAM SortedByCoordinate")
star.in <- paste0("--readFilesIn ", paste0(files.in, collapse=" "))
star.out <- paste0("--outFileNamePrefix ", star.outbase)
star.cmd  <- paste0(star.exec, " ", star.in, " ", star.opts, " ", star.out)

execute(star.cmd, star.bam)

#Creating a samtools index
execute(paste0("samtools index ", star.bam), paste0(star.bam, ".bai"))

#flagstats coverage statistics
execute(paste0("samtools flagstats ", star.bam))
execute(paste0("samtools coverage ", star.bam))


#We use these sorted bam files as input into the featureCounts function of RSubread in the DESeq script


#THese steps are for visualisation of the bam files in IGV just to see how RNAseq readsz look when aligned to the reference genome.



#Remove duplicate reads using picard tools
p.bam <- paste0(star.outbase, "Aligned.sortedByCoord.RD.out.bam")
metrics.out <- paste0(star.outbase, "_metrics.txt")
p.exec <- "picard"
p.in <- paste0("-I ", star.bam)
p.out <- paste0("-O ", p.bam, " -M ", metrics.out)
p.opts <- paste0("--REMOVE_DUPLICATES true")
p.cmd <- paste0(p.exec, " MarkDuplicates ", p.opts," ", p.in, " ", p.out)

execute(p.cmd, p.bam)
execute(paste0("samtools flagstats ", p.bam))
execute(paste0("samtools coverage ", p.bam))

#Add read group and sample run, library, and name
rg.bam <- "../output/"cmdlineargs[1]".aln/"cmdlineargs[1]"Aligned.sortedByCoord.RD.RG.out.bam"
rg.opts <- paste0("PL=ILLUMINA PU=run LB=", gsub("SRR", "", input.base), " SM=", input.base)
p.cmd <- paste0(p.exec, " AddOrReplaceReadGroups I=../output/"cmdlineargs[1]".aln/"cmdlineargs[1]"Aligned.sortedByCoord.RD.out.bam", " O=../output/"cmdlineargs[1].aln/"cmdlineargs[1]"Aligned.sortedByCoord.RD.RG.out.bam", " ", rg.opts)
execute(p.cmd)

execute(paste0("samtools index ", rg.bam), paste0(rg.bam, ".bai"))