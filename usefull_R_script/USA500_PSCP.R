### Script.sh FILES ###
#Get all the genomes files(gz), put all of them and the reference in a folder.
#rename all the files if they start with number.
#scp to the cluster
rm(list=ls())
library(stringr)

setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq")
#ReFF_name = ("USA300_FPR3757.fasta")
ReFF_name = ("USA500_2395.fasta")
#ReFF_name = ("CFSAN007894.fasta")
#ReFF_name = ("JH1.fasta")
ReFF_namelabels = sapply(ReFF_name, function(s) {
  s = gsub(".fasta", "", s)#replace ".fasta" with space 
})

genome.tested <-grep("fasta",ReFF_name)# find the one match "FPR"

# run the following comand on terminal of SSH
# index reference
# module load bwa/0.6.2
# module load samtools/1.3
# bwa index ReFF
# samtools faidx ReFF
# java -jar /local/apps/picard-tools/1.88/CreateSequenceDictionary.jar R=JH1.fasta O=JH1.dict


fastq_name = scan("USA500_import_files.txt",what=character(0),nlines=1,sep="\t")

# get all the name and copy to excel and column to text and take the last part
#replace data after by ";" to get all the name and dedup the data and transpose the data

#create different kinds of labels for each sequences
for(k in 1){
  newseqlabels = sapply(fastq_name, function(s) {
    s = gsub(";", "", s)
  })
  newseqlabels.sh = sapply(fastq_name, function(s) {
    s = gsub(";", ".sh", s)
  })
  newseqlabels.sge = sapply(fastq_name, function(s) {
    s = gsub(";", ".sge", s)
  })
  newseqlabels.vcf = sapply(fastq_name, function(s) {
    s = gsub(";", ".txt", s)
  })
  newseqlabels.exonic = sapply(fastq_name, function(s) {
    s = gsub(";", ".txt.exonic_variant_function", s)
  })
  newseqlabels.exonicOutput = sapply(fastq_name, function(s) {
    s = gsub(";", ".exonic.txt", s)
  })
  newseqlabels.intergenic = sapply(fastq_name, function(s) {
    s = gsub(";", ".txt.variant_function", s)
  })
  newseqlabels.intergenicOutput = sapply(fastq_name, function(s) {
    s = gsub(";", ".intergenic.txt", s)
  })
  
  newseqlabels.invalid = sapply(fastq_name, function(s) {
    s = gsub(";", ".txt.invalid_input", s)
  })
  newseqlabels.invalidoutput = sapply(fastq_name, function(s) {
    s = gsub(";", "invalid.txt", s)
  })
  newseqlabels.invalid.exonic = sapply(fastq_name, function(s) {
    s = gsub(";", "invalid.txt.exonic_variant_function", s)
  })
  newseqlabels.invalid.exonic.Output = sapply(fastq_name, function(s) {
    s = gsub(";", "invalid_exonic.txt", s)
  })
  newseqlabels.invalid.intergenic = sapply(fastq_name, function(s) {
    s = gsub(";", "invalid.txt.variant_function", s)
  })
  newseqlabels.invalid.intergenic.Output = sapply(fastq_name, function(s) {
    s = gsub(";", "invalid_intergenic.txt", s)
  })
  
  nSeq = length(newseqlabels)
}

#first, create a template of genome mapping, then replace the variables and crate bash files for each of the genome.
#create a directory of reference and create a directory of sequence 

setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/mapping_program/USA500_2395")
for(k in genome.tested){
  Script = paste(ReFF_namelabels[k],"script.sh", sep = "")
  cat(paste("#!/bin/bash"),sep="\n",file=Script,append=TRUE)
  cat(paste("bwa aln -t 2 ../reference/ReFF ../fastq/name_R1_001.fastq.gz > name_pair_1.sai1"),sep="\n",file=Script,append=TRUE)
  cat(paste("bwa aln -t 2 ../reference/ReFF ../fastq/name_R2_001.fastq.gz > name_pair_2.sai2"),sep="\n",file=Script,append=TRUE)
  cat(paste('bwa sampe -r "@RG\\tID:IDa\\tSM:SM\\tPL:Illumina" ../reference/ReFF name_pair_1.sai1 name_pair_2.sai2 ../fastq/name_R1_001.fastq.gz ../fastq/name_R2_001.fastq.gz > name.sam'),sep="\n",file=Script,append=TRUE)
  cat(paste("samtools view -b -S -t ../reference/ReFF -o name.bam name.sam"),sep="\n",file=Script,append=TRUE)
  cat(paste("samtools sort name.bam -o name.sortedbam.bam"),sep="\n",file=Script,append=TRUE)
  cat(paste("samtools index name.sortedbam.bam"),sep="\n",file=Script,append=TRUE)
  cat(paste("java -jar /local/apps/gatk/2.4-9/GenomeAnalysisTK.jar -R ../reference/ReFF -T HaplotypeCaller -I name.sortedbam.bam -o name.raw.snps.indels.vcf -stand_call_conf 30"),sep="\n",file=Script,append=TRUE)
  cat(paste("rm name.sam"),sep="\n",file=Script,append=TRUE)
  cat(paste("rm name_pair_1.sai1"),sep="\n",file=Script,append=TRUE)
  cat(paste("rm name_pair_2.sai2"),sep="\n",file=Script,append=TRUE)
  cat(paste("rm name.bam"),sep="\n",file=Script,append=TRUE)
  cat(paste("java -jar /local/apps/gatk/2.4-9/GenomeAnalysisTK.jar -T DepthOfCoverage -R ../reference/ReFF -o name_genome_coverage -I name.sortedbam.bam"),sep="\n",file=Script,append=TRUE)
  
  
  variants = paste(ReFF_namelabels[k],"script.sh", sep = "")
  for(i in 1:nSeq) {
    x <- readLines(variants)
    y <- gsub(pattern = "name", replace = newseqlabels[i], x)
    y <- gsub(pattern = "ReFF", replace = ReFF_name[k], y)
    cat(y, file=paste(ReFF_namelabels[k],"_",newseqlabels.sh[i], sep=""), sep="\n")
    cat("Output",i,"of",nSeq,"\n")
  }
  
  for(i in 1:nSeq) {
    x <- readLines(variants)
    y <- gsub(pattern = "name", replace =newseqlabels[i], x)
    y <- gsub(pattern = "ReFF", replace = ReFF_name[k], y)
    cat(y, file=paste(ReFF_namelabels[k],"_",i,".sh",sep=""), sep="\n")
    cat("Output",i,"of",nSeq,"\n")
  }
  
  ### Script.sge FILE ###
  # create these files to run each of the bash files.
  
  
  SGE = paste(ReFF_namelabels[k],"script.sge", sep = "")
  cat(paste("#!/bin/bash"),sep="\n",file=SGE,append=TRUE)
  cat(paste("#$ -S /bin/bash"),sep="\n",file=SGE,append=TRUE)
  cat(paste("#$ -cwd"),sep="\n",file=SGE,append=TRUE)
  cat(paste("#$ -j y"),sep="\n",file=SGE,append=TRUE)
  cat(paste("#$ -M dan.ding@nyumc.org"),sep="\n",file=SGE,append=TRUE)
  cat(paste("#$ -m e"),sep="\n",file=SGE,append=TRUE)
  cat(paste("#$ -t 1-",nSeq,sep=""),sep="\n",file=SGE,append=TRUE)
  cat(paste('echo "Loading bwa/0.6.2"'),sep="\n",file=SGE,append=TRUE)
  cat(paste("module load bwa/0.6.2"),sep="\n",file=SGE,append=TRUE)
  cat(paste('echo "Loading samtools/1.3"'),sep="\n",file=SGE,append=TRUE)
  cat(paste("module load samtools/1.3"),sep="\n",file=SGE,append=TRUE)
  cat(paste('echo "${SGE_TASK_ID}.sh"'),sep="\n",file=SGE,append=TRUE)
  cat(paste("chmod +x /ifs/data/shopsinlab/Dan/PSCP_PROJECT/dep_seq/",ReFF_namelabels[1],"/",ReFF_namelabels[k],"_${SGE_TASK_ID}.sh", sep = ""),sep="\n",file=SGE,append=TRUE)
  cat(paste("./",ReFF_namelabels[k],"_${SGE_TASK_ID}.sh", sep = ""),sep="\n",file=SGE,append=TRUE)
  
}  
#for(i in 1:nSeq) {
# x <- readLines(SGE)
#y <- gsub(pattern = "name", replace = newseqlabels.sh[i], x)
#cat(y, file=paste(ReFF_namelabels[k],"_",newseqlabels.sge[i], sep=""), sep="\n")
#cat("Output",i,"of",nSeq,"\n")
#}


#submit all the SGE to the cluster.
#Qsub = paste(ReFF_namelabels[k],"_qsub.sh", sep = "")
#cat(paste("#!/bin/bash"),sep="\n",file=Qsub,append=TRUE)
#cat(paste("#$ -S /bin/bash"),sep="\n",file=Qsub,append=TRUE)
#cat(paste("#$ -cwd"),sep="\n",file=Qsub,append=TRUE)
#cat(paste("#$ -j y"),sep="\n",file=Qsub,append=TRUE)
#cat(paste("#$ -M dan.ding@nyumc.org"),sep="\n",file=Qsub,append=TRUE)
#cat(paste("#$ -m be"),sep="\n",file=Qsub,append=TRUE)
#for(i in 1:nSeq) {
#cat(paste("qsub -hard -l mem_free=40G -l h_vmem=40G -l mem_token=40G ", ReFF_namelabels[k],"_",newseqlabels.sge[i], sep=""),sep="\n",file=Qsub, append = TRUE)
#}

#}

#scp everything to the cluster
# We  have to submit job to cluster and let it run from the back, not just running any program from the cluster.
#we submit the job by using qsub each of the SGE files to the cluster by running Qsub files. Then it is run each of the 
#SGE file, then each of the bash file will be ran and finished Genome mapping.
cat(paste("scp *.s* dingd02@phoenix.med.nyu.edu:/ifs/data/shopsinlab/Dan/PSCP_PROJECT/dep_seq/USA500_2395", sep = ""))


#alignment

cat(paste("scp dingd02@phoenix.med.nyu.edu:/ifs/data/shopsinlab/Dan/PSCP_PROJECT/dep_seq/USA500_2395","/*.vcf .", sep = ""))
cat(paste("scp dingd02@phoenix.med.nyu.edu:/ifs/data/shopsinlab/Dan/PSCP_PROJECT/dep_seq/USA300_FPR3757","/*.sortedbam.ba* .", sep = ""))


### VCFtool FILE ###
#need to run it everytime when we change sequences
# Run the code under the file with vcf-query and get ouput in annovar files
setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/vcftools/src/perl")
VCF = "VCFtool_new.sh"
cat(paste("#!/bin/bash"),sep="\n",file=VCF,append=TRUE)
nSeq = length(newseqlabels)
for(i in 1:nSeq) {
  cat(paste("perl vcf-query /Users/dingd02/Desktop/PSCP_PROJECT/Dep_seq/VCF_files/",ReFF_namelabels[1],"/", newseqlabels[i],'.raw.snps.indels.vcf -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/DP\\t%INFO/AF\\n" > /Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/annovar/',newseqlabels[i], '.txt', sep=""),sep="\n",file=VCF, append = TRUE)
  cat("Output",i,"of",nSeq,"\n")
}

my_VCF <- list.files(pattern = "vcf")


for(k in 1){
  newseqlabels.new = sapply(my_VCF, function(s) {
    s = gsub(".raw.snps.indels.vcf", ";", s)
  })
  newseqlabels = sapply(newseqlabels.new, function(s) {
    s = gsub(";", "", s)
  })
  newseqlabels.sh = sapply(newseqlabels.new, function(s) {
    s = gsub(";", ".sh", s)
  })
  newseqlabels.sge = sapply(newseqlabels.new, function(s) {
    s = gsub(";", ".sge", s)
  })
  newseqlabels.vcf = sapply(newseqlabels.new, function(s) {
    s = gsub(";", ".txt", s)
  })
  newseqlabels.exonic = sapply(newseqlabels.new, function(s) {
    s = gsub(";", ".txt.exonic_variant_function", s)
  })
  newseqlabels.exonicOutput = sapply(newseqlabels.new, function(s) {
    s = gsub(";", ".exonic.txt", s)
  })
  newseqlabels.intergenic = sapply(newseqlabels.new, function(s) {
    s = gsub(";", ".txt.variant_function", s)
  })
  newseqlabels.intergenicOutput = sapply(newseqlabels.new, function(s) {
    s = gsub(";", ".intergenic.txt", s)
  })
  
  newseqlabels.invalid = sapply(newseqlabels.new, function(s) {
    s = gsub(";", ".txt.invalid_input", s)
  })
  newseqlabels.invalidoutput = sapply(newseqlabels.new, function(s) {
    s = gsub(";", "invalid.txt", s)
  })
  newseqlabels.invalid.exonic = sapply(newseqlabels.new, function(s) {
    s = gsub(";", "invalid.txt.exonic_variant_function", s)
  })
  newseqlabels.invalid.exonic.Output = sapply(newseqlabels.new, function(s) {
    s = gsub(";", "invalid_exonic.txt", s)
  })
  newseqlabels.invalid.intergenic = sapply(newseqlabels.new, function(s) {
    s = gsub(";", "invalid.txt.variant_function", s)
  })
  newseqlabels.invalid.intergenic.Output = sapply(newseqlabels.new, function(s) {
    s = gsub(";", "invalid_intergenic.txt", s)
  })
  
  nSeq = length(newseqlabels)
}

cat("cd /Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/vcftools/src/perl")
cat("chmod +x VCFtool.sh")


#check depth of coverage

setwd(paste("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/",ReFF_namelabels[genome.tested],sep=""))
for (i in 1:nSeq) {
  test = read.delim(newseqlabels.vcf[i],as.is=TRUE, header = FALSE)
  print(summary(test$V6))
}

#make column1 as 1 and create another column behind column 2 and  same as column 2
#clean-up the files
for(i in 1:nSeq) {
  VCF = read.delim(newseqlabels.vcf[i], as.is=TRUE, header = FALSE)
  VCF$V1 <- "1"
  VCF["V8"] <- NA
  VCF$V8 <- VCF$V2
  VCF <- VCF[,c(1,2,8,3,4,5,6,7)]
  write.table(VCF, file = newseqlabels.vcf[i], sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
}

setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/annovar")
### Annovar FILE ###

Anno = "Annovar.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Anno,append=TRUE)
for(i in 1:nSeq) {
  cat(paste("perl annotate_variation.pl ", newseqlabels.vcf[i]," --buildver ", ReFF_namelabels[1]," ", ReFF_namelabels[1],"/", sep=""),sep="\n",file=Anno, append = TRUE)
  cat("Output",i,"of",nSeq,"\n")
}
cat("chmod +x /Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/annovar/Annovar.sh")
cat(paste("mv *exonic* ../", ReFF_namelabels[1],"/exonic/. | mv *inval* ../", ReFF_namelabels[1],"/invalid/. | mv *variant_f* ../", ReFF_namelabels[1],"/intergenic/. | rm *.log ", sep = ""))
getwd()

# Invalid output
setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/USA300_FPR3757/invalid")
for(i in 1:nSeq) {
  invalid = read.delim(newseqlabels.invalid[i], as.is=TRUE, header = FALSE)
  invalid <- invalid[,c(1,2,3,5,4,6,7,8)]
  write.table(invalid, file = newseqlabels.invalidoutput[i], sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
}


#  mv ../USA300_FPR3757/invalid/*invalid.txt .
paste("mv ../", ReFF_namelabels[1],"/invalid/*invalid.txt .", sep = "")

# Annovar on Invalid
setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/annovar/")

AnnoInvl = "Annovar_invalid.sh"
cat(paste("#!/bin/bash"),sep="\n",file=AnnoInvl,append=TRUE)
for(i in 1:nSeq) {
  cat(paste("perl annotate_variation.pl ", newseqlabels.invalidoutput[i]," --buildver ", ReFF_namelabels[1]," ", ReFF_namelabels[1],"/", sep=""),sep="\n",file=AnnoInvl, append = TRUE)
  cat("Output",i,"of",nSeq,"\n")
}
cat("chmod +x /Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/annovar/Annovar_invalid.sh")
cat(paste("mv *exonic* ../", ReFF_namelabels[1],"/exonic/. | mv *variant_f* ../", ReFF_namelabels[1],"/intergenic/. | rm *.log | rm *_input", sep = ""))



### OPEN exonic_variant_function_files ###

getwd()
setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/USA300_FPR3757/exonic/")
for(i in 1:nSeq) {
  d1 <- read.table(newseqlabels.exonic[i], sep = "\t", fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE)
  d2 = data.frame(d1[3])
  d3 = str_split_fixed(d2$V3, ":", n = 5)
  d = cbind.data.frame(d1[,1:2], d3[,1:5], d1[,4:11])
  d = data.frame(d)
  # Work on Good quality sequences
  exonic <- subset.data.frame(d, d$V10>=50)
  exonic <- subset.data.frame(exonic, exonic$V9>=800)
  exonic <- exonic[,-10]
  exonic <- exonic[,-8]
  exonic <- exonic[,-5]
  exonic <- exonic[,-1]
  exonic[,1] <- gsub(" SNV","",exonic[,1])
  exonic[,4] <- gsub("c.","",exonic[,4])
  exonic[,5] <- gsub("p.","",exonic[,5])
  exonic[,5] <- gsub(",","",exonic[,5])
  write.table(exonic, file = newseqlabels.exonicOutput[i], sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
}

for(i in 1:nSeq) {
  print(i)
  tryCatch({
    d1 <- read.table(newseqlabels.invalid.exonic[i], sep = "\t", fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE)
    d2 = data.frame(d1[3])
    d3 = str_split_fixed(d2$V3, ":", n = 5)
    d = cbind.data.frame(d1[,1:2], d3[,1:5], d1[,4:11])
    d = data.frame(d)
    # Work on Good quality sequences
    invalidexonic <- subset.data.frame(d, d$V10>=50)
    invalidexonic <- subset.data.frame(invalidexonic, invalidexonic$V9>=800)
    invalidexonic <- invalidexonic[,-10]
    invalidexonic <- invalidexonic[,-8]
    invalidexonic <- invalidexonic[,-5]
    invalidexonic <- invalidexonic[,-1]
    invalidexonic <- invalidexonic[,c(1,2,3,4,5,6,8,7,9,10,11)]
    write.table(invalidexonic, file = newseqlabels.invalid.exonic.Output[i], sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

for(i in 1:length(nSeq)){
  x <- read.table(newseqlabels.exonicOutput[i], sep = "\t", fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE)
  y <- read.table(newseqlabels.invalid.exonic.Output[i], sep = "\t", fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE)
  z <- rbind(x,y)
  write.table(z, file = newseqlabels.exonicOutput[i], sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
}


### OPEN intergenic_files ###

setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/USA300_FPR3757/intergenic/")

for(i in 1:nSeq) {
  noexonic <- read.table(newseqlabels.intergenic[i], sep = "\t", fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE)
  noexonic <- subset.data.frame(noexonic, noexonic[,1]!="exonic")
  noexonic[,1] <- gsub("upstream", "promoter", noexonic[,1])
  noexonic <- subset.data.frame(noexonic, noexonic[,9]>=50)
  noexonic <- subset.data.frame(noexonic, noexonic[,8]>=800)
  noexonic <- noexonic[,-5]
  noexonic <- noexonic[,-3]
  write.table(noexonic, file = newseqlabels.intergenicOutput[i], sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

for(i in 1:nSeq) {
  tryCatch({
    invalid.noexonic <- read.table(newseqlabels.invalid.intergenic[i], sep = "\t", fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE)
    invalid.noexonic <- subset.data.frame(invalid.noexonic, invalid.noexonic[,1]!="exonic")
    invalid.noexonic[,1] <- gsub("upstream", "promoter", invalid.noexonic[,1])
    invalid.noexonic <- subset.data.frame(invalid.noexonic, invalid.noexonic[,9]>=50)
    invalid.noexonic <- subset.data.frame(invalid.noexonic, invalid.noexonic[,8]>=800)
    invalid.noexonic <- invalid.noexonic[,-5]
    invalid.noexonic <- invalid.noexonic[,-3]
    invalid.noexonic <- invalid.noexonic[,c(1,2,3,5,4,6,7,8)]   
    write.table(invalid.noexonic, file = newseqlabels.invalid.intergenic.Output[i], sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

for(i in 1:length(nSeq)){
  a <- read.table(newseqlabels.intergenicOutput[i], sep = "\t", fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE)
  b <- read.table(newseqlabels.invalid.intergenic.Output[i], sep = "\t", fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE)
  c <- rbind(a,b)
  write.table(c, file = newseqlabels.intergenicOutput[i], sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
}














