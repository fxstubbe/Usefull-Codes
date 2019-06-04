### Script.sh FILES ###
#Get all the genomes files(gz), put all of them and the reference in a folder.
#rename all the files if they start with number.
#scp to the cluster
rm(list=ls())
library(stringr)

setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq")
#ReFF_name = ("USA300_FPR3757.fasta")
#ReFF_name = ("USA500_2395.fasta")
#ReFF_name = ("CFSAN007894.fasta")
#ReFF_name = ("JH1.fasta")
ReFF_name=c("USA300_FPR3757.fasta","USA500_2395.fasta","CFSAN007894.fasta","JH1.fasta")
ReFF_namelabels = sapply(ReFF_name, function(s) {
  s = gsub(".fasta", "", s)#replace ".fasta" with space 
})

genome.tested <-grep("300",ReFF_name)# find the one match "FPR"

# run the following comand on terminal of SSH
# index reference
# module load bwa/0.6.2
# module load samtools/1.3
# bwa index CFSAN007894.fasta
# samtools faidx CFSAN007894.fasta
# java -jar /local/apps/picard-tools/1.88/CreateSequenceDictionary.jar R=CFSAN007894.fasta O=CFSAN007894.dict


fastq_name = scan(paste(ReFF_namelabels[genome.tested],"_part2_import_files.txt",sep=""),what=character(0),nlines=1,sep="\t")

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

setwd(paste("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/mapping_program/",ReFF_namelabels[genome.tested],"/part2",sep=""))


getwd()
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
    cat(y, file=paste(ReFF_namelabels[k],"_part2_",i,".sh",sep=""), sep="\n")
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
  cat(paste("chmod +x /ifs/data/shopsinlab/Dan/PSCP_PROJECT/dep_seq/",ReFF_namelabels[k],"/",ReFF_namelabels[k],"_part2_${SGE_TASK_ID}.sh", sep = ""),sep="\n",file=SGE,append=TRUE)
  cat(paste("./",ReFF_namelabels[k],"_part2_${SGE_TASK_ID}.sh", sep = ""),sep="\n",file=SGE,append=TRUE)
  
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
cat(paste("scp *.s* dingd02@phoenix.med.nyu.edu:/ifs/data/shopsinlab/Dan/PSCP_PROJECT/dep_seq/",ReFF_namelabels[genome.tested],sep = ""))


#alignment

cat(paste("scp dingd02@phoenix.med.nyu.edu:/ifs/data/shopsinlab/Dan/PSCP_PROJECT/dep_seq/",ReFF_namelabels[genome.tested],"/*.vcf .", sep = ""))
#cat(paste("scp dingd02@phoenix.med.nyu.edu:/ifs/data/shopsinlab/Dan/PSCP_PROJECT/dep_seq/",ReFF_namelabels[genome.tested],"/*.sortedbam.ba* .", sep = ""))

######################only run this part if some of the VCF failed....
setwd(paste("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/VCF_files/",ReFF_namelabels[genome.tested],sep=""))
getwd()
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
  
  nSeq = length(newseqlabels)
}
#####################################





### VCFtool FILE ###
#need to run it everytime when we change sequences
# Run the code under the file with vcf-query and get ouput in annovar files
setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/vcftools/src/perl")
VCF = paste("VCFtool_",ReFF_namelabels[genome.tested],".sh",sep="")
cat(paste("#!/bin/bash"),sep="\n",file=VCF,append=TRUE)
nSeq = length(newseqlabels)
for(i in 1:nSeq) {
  cat(paste("perl vcf-query /Users/dingd02/Desktop/PSCP_PROJECT/Dep_seq/VCF_files/",ReFF_namelabels[genome.tested],"/", newseqlabels[i],'.raw.snps.indels.vcf -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/DP\\t%INFO/AF\\n" > /Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/annovar/',newseqlabels[i], '.txt', sep=""),sep="\n",file=VCF, append = TRUE)
  cat("Output",i,"of",nSeq,"\n")
}


cat("cd /Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/vcftools/src/perl")



#check depth of coverage

setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/annovar/")





for(i in 1:nSeq) {
  VCF = read.delim(newseqlabels.vcf[i], as.is=TRUE, header = FALSE)
  #VCF = test2[[i]]
  VCF$V1 <- "1"
  VCF["V8"] <- NA
  VCF$V8 <- VCF$V2
  VCF <- VCF[,c(1,2,8,3,4,5,6,7)]
  write.table(VCF, file = newseqlabels.vcf[i], sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
}




test2 <- list()
summary.list <- list()
for (i in 1:nSeq) {
  test = read.delim(newseqlabels.vcf[i],as.is=TRUE, header = FALSE)
  x1 <- as.numeric(as.character(summary(test$V7))[2])/2
  x2 <- test[which(test$V7 >= x1),]
  x3 <- x2[which(x2$V8 == 1 | x2$V8 =="1.00"),]
  y1 <- as.numeric(as.character(summary(test$V7)))
  test2[[i]] <- x3[which(x3$V6 >= 999),]
  summary.list[[i]] <- y1
  final=test2[[i]];
  write.table(final, file = newseqlabels.vcf[i], sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}

df.summary <- as.data.frame(do.call(rbind, summary.list))
colnames(df.summary) <- c("Minimum read #", "1st Qu.",  "Median reads",    "Mean reads",  "3rd Qu.",    "Max read #")
rownames(df.summary) <- as.character(sapply(as.character(newseqlabels.vcf), function(s) gsub(".txt","",s)))
write.csv(df.summary, "quality_check_coverage.csv")

to.remove <- rownames(df.summary)[which(df.summary$`Mean reads` < 60)]  # number to modulate

newseqlabels.new <- newseqlabels.new[-grep(to.remove,newseqlabels.new)]

for(k in 1){
        newseqlabels.vcf = sapply(newseqlabels.new, function(s) {
           s = gsub(";", ".txt", s)
        })
        newseqlabels = sapply(newseqlabels.new, function(s) {
          s = gsub(";", "", s)
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

#make column1 as 1 and create another column behind column 2 and  same as column 2
#clean-up the files




### Annovar FILE ###

Anno = paste("Annovar_",ReFF_namelabels[genome.tested],".sh",sep="")
cat(paste("#!/bin/bash"),sep="\n",file=Anno,append=TRUE)
for(i in 1:nSeq) {
  cat(paste("perl annotate_variation.pl ", newseqlabels.vcf[i]," --buildver ", ReFF_namelabels[genome.tested]," ", ReFF_namelabels[genome.tested],"/", sep=""),sep="\n",file=Anno, append = TRUE)
  cat("Output",i,"of",nSeq,"\n")
}
cat("chmod +x /Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/annovar/Annovar.sh")
cat(paste("mv *exonic* ../", ReFF_namelabels[genome.tested],"/exonic/. | mv *inval* ../", ReFF_namelabels[genome.tested],"/invalid/. | mv *variant_f* ../", ReFF_namelabels[genome.tested],"/intergenic/. | rm *.log ", sep = ""))
getwd()

# Invalid output
setwd(paste("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/", ReFF_namelabels[genome.tested],"/invalid",sep=""))
for(i in 1:nSeq) {
  invalid = read.delim(newseqlabels.invalid[i], as.is=TRUE, header = FALSE)
  invalid <- invalid[,c(1,2,3,5,4,6,7,8)]
  write.table(invalid, file = newseqlabels.invalidoutput[i], sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
}


#  mv ../USA300_FPR3757/invalid/*invalid.txt .
cat(paste("mv ../", ReFF_namelabels[genome.tested],"/invalid/*invalid.txt .", sep = ""))

# Annovar on Invalid
setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/annovar/")

AnnoInvl = paste("Annovar_invalid_",ReFF_namelabels[genome.tested],".sh",sep="")

cat(paste("#!/bin/bash"),sep="\n",file=AnnoInvl,append=TRUE)
for(i in 1:nSeq) {
  cat(paste("perl annotate_variation.pl ", newseqlabels.invalidoutput[i]," --buildver ", ReFF_namelabels[genome.tested]," ", ReFF_namelabels[genome.tested],"/", sep=""),sep="\n",file=AnnoInvl, append = TRUE)
  cat("Output",i,"of",nSeq,"\n")
}
cat("chmod +x /Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/annovar/Annovar_invalid.sh")
cat(paste("mv *exonic* ../", ReFF_namelabels[genome.tested],"/exonic/. | mv *variant_f* ../", ReFF_namelabels[genome.tested],"/intergenic/. | rm *.log | rm *_input", sep = ""))



### OPEN exonic_variant_function_files ###

getwd()
setwd(paste("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/",ReFF_namelabels[genome.tested],"/exonic/",sep=""))
for(i in 1:nSeq) {
  d1 <- read.table(newseqlabels.exonic[i], sep = "\t", fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE)
  d2 = data.frame(d1[3])
  d3 = str_split_fixed(d2$V3, ":", n = 5)
  d = cbind.data.frame(d1[,1:2], d3[,1:5], d1[,4:11])
  d = data.frame(d)
  # Work on Good quality sequences
  exonic <- subset.data.frame(d, d$V11==1.0)
  # exonic <- subset.data.frame(exonic, exonic$V10>=50)
  exonic <- subset.data.frame(exonic, exonic$V9>=1000)
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
    invalidexonic <- subset.data.frame(d, d$V11==1.0)
    invalidexonic <- subset.data.frame(invalidexonic, invalidexonic$V10>=50)
    invalidexonic <- subset.data.frame(invalidexonic, invalidexonic$V9>=1000)
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

setwd(paste("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/",ReFF_namelabels[genome.tested],"/intergenic",sep=""))

for(i in 1:nSeq) {
  noexonic <- read.table(newseqlabels.intergenic[i], sep = "\t", fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE)
  noexonic <- subset.data.frame(noexonic, noexonic[,1]!="exonic")
  noexonic[,1] <- gsub("upstream", "promoter", noexonic[,1])
  noexonic <- subset.data.frame(noexonic, noexonic[,10]==1)
  noexonic <- subset.data.frame(noexonic, noexonic[,9]>=50)
  noexonic <- subset.data.frame(noexonic, noexonic[,8]>=1000)
  noexonic <- noexonic[,-5]
  noexonic <- noexonic[,-3]
  write.table(noexonic, file = newseqlabels.intergenicOutput[i], sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

for(i in 1:nSeq) {
  tryCatch({
    invalid.noexonic <- read.table(newseqlabels.invalid.intergenic[i], sep = "\t", fill = TRUE, quote = "", row.names = NULL, stringsAsFactors = FALSE)
    invalid.noexonic <- subset.data.frame(invalid.noexonic, invalid.noexonic[,1]!="exonic")
    invalid.noexonic[,1] <- gsub("upstream", "promoter", invalid.noexonic[,1])
    invalid.noexonic <- subset.data.frame(invalid.noexonic, invalid.noexonic[,10]==1)
    invalid.noexonic <- subset.data.frame(invalid.noexonic, invalid.noexonic[,9]>=50)
    invalid.noexonic <- subset.data.frame(invalid.noexonic, invalid.noexonic[,8]>=1000)
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


#variant.................................
library(data.table)
library(plyr)

setwd(paste("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/",ReFF_namelabels[genome.tested],"/exonic/",sep=""))

my_files <- list.files(pattern = "\\.exonic.txt")
my_names1 <- gsub(".ex.*","",my_files)
my_names <- gsub("merge_","",my_names1)
my_names
my_list <- lapply(my_files, function(i){
  fread( i, sep = "\t", header = F, data.table = F, fill = T, blank.lines.skip =TRUE)
})
names(my_list) <- my_names
column.names <- c("mutation nature", "product name", "locus tag", 
                  "nucleotide change","protein change", "genomic position", 
                  "ref", "alt","quality score", "depth of coverage", 
                  "allelic frequency")
my_list <- lapply(my_list, setNames, column.names)

for(o in 1:length(my_list)){
  b <- which(is.na(my_list[[o]]$`allelic frequency`))
  my_list[[o]]$`allelic frequency`[b] <- my_list[[o]]$`depth of coverage`[b]
  my_list[[o]]$`depth of coverage`[b] <- my_list[[o]]$`quality score`[b]
  my_list[[o]]$`quality score`[b] <- my_list[[o]]$alt[b]
  my_list[[o]]$alt[b] <- my_list[[o]]$ref[b]
  my_list[[o]]$ref[b] <- my_list[[o]]$`genomic position`[b]
  my_list[[o]]$`genomic position`[b] <- my_list[[o]]$`protein change`[b]
}
total.df = data.table::rbindlist(my_list)
gen.mat <- matrix(data = NA, nrow = length(total.df$`genomic position`), ncol = length(my_list))
colnames(gen.mat) <- my_names
total.df <- as.data.frame(cbind(gen.mat, total.df))
noDup.df <- as.data.frame(subset.data.frame(total.df, !duplicated(total.df$`genomic position`)))
#good.qul<-noDup.df
#good.qul<-good.qul[!grepl("phage",good.qul$`product name`,ignore.case = T),]
#good.qul<-good.qul[!grepl("transposase",good.qul$`product name`,ignore.case = T),]
#good.qul<-good.qul[!grepl("integrase",good.qul$`product name`,ignore.case = T),]

good.qul <- noDup.df[nchar(noDup.df$`alt`)==1,]
good.qul <- good.qul[nchar(good.qul$`ref`)==1,]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 34173,57914 ),] #SCCMEC
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 35877,36551 ),] 
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 42442,43965 ),] 
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 57915,88900 ),] #ACME
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 67966,68790 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 78432,79445 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 318970,319374 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 441504,473466 ),] #vSAalpha
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 443533,443805 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 681555,682115 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 881837,895808 ),] #SaPI5
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 881996,883216 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 1443183,1443668 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 1545912,1592050 ),] #PVL
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 1590733,1591938 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 1802471,1803079 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 1917170,1917655 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 1924777,1959376 ),] #vSAbeta
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 1993637,1994956 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 2084658,2127720 ),] #beta phage
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 2126612,2127649 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 2290418,2291620 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 2431859,2432173 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 2478538,2478926 ),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c( 2862387,2862884 ),]
atBeginning <- Sys.time()
then <- Sys.time()
for (i in 1:length(my_list)){
  for (j in 1:length(good.qul$`genomic position`)){
    if(good.qul$`genomic position`[j] %in% my_list[[i]]$`genomic position`) {
      good.qul[j,i] <- 1
    }else
      good.qul[j,i] <- 0
  }
}
print(Sys.time()-then) # Time difference of 3.65937 hours for ~364 strains and 59040 SNPs
mut.nature <- which(colnames(good.qul) == "mutation nature")
phylo.df = good.qul
ref <- which(colnames(phylo.df) == "ref")
alt <- which(colnames(phylo.df) == "alt")
atBeginning <- Sys.time()
then <- Sys.time()
for (i in 1:(mut.nature -1)){
  for (j in 1:nrow(phylo.df)){
    if(phylo.df[j,i] == 1) {
      phylo.df[j,i] <- phylo.df[j,alt]
    }else
      phylo.df[j,i] <- phylo.df[j,ref]
  }
}
print(Sys.time()-then) # Time difference of 3.009364 hours for ~364 strains and 59040 SNPs
MST.df <- good.qul[,c(mut.nature:length(good.qul),1:(mut.nature -1))]
phylo.df <- phylo.df[,c(mut.nature:length(phylo.df),1:(mut.nature -1))]
write.table(MST.df, file = "180326_PSCP_annnotation_exonic.txt", sep="\t", row.names = FALSE, quote = FALSE) 



#Intergenic

setwd(paste("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/",ReFF_namelabels[genome.tested],"/intergenic",sep=""))

my_files.interg <- list.files(pattern = "\\.intergenic")
my_interg1 <- gsub(".int.*","",my_files.interg)
my_interg <- gsub("merge_","",my_interg1)
my_list.interg <- lapply(my_files.interg, function(i){
  fread( i, sep = "\t", header=F, data.table = F)
})
names(my_list.interg) <- my_interg
column.interg <- c("mutation nature", "locus tag","genomic position", 
                   "ref", "alt","quality score", "depth of coverage", 
                   "allelic frequency")
my_list.interg <- lapply(my_list.interg, setNames, column.interg)
total.interg = data.table::rbindlist(my_list.interg)
gen.mat.interg <- matrix(data = NA, nrow = length(total.interg$`genomic position`), ncol = length(my_list.interg))
colnames(gen.mat.interg) <- my_interg
total.interg <- as.data.frame(cbind(gen.mat.interg, total.interg))
noDup.interg <- as.data.frame(subset.data.frame(total.interg, !duplicated(total.interg$`genomic position`)))
#good.interg <- noDup.interg
#good.interg <-good.interg [!grepl("phage",good.qul$`product name`,ignore.case = T),]
#good.interg <-good.interg [!grepl("transposase",good.qul$`product name`,ignore.case = T),]
#good.interg <-good.interg [!grepl("integrase",good.qul$`product name`,ignore.case = T),]

good.interg <- noDup.interg[nchar(noDup.interg$`alt`)==1,]
good.interg <- good.interg[nchar(good.interg$`ref`)==1,]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c( 34173,57914 ),] #SCCMEC
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c( 57915,88900 ),] #ACME
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c( 441504,473466 ),] #vSAalpha
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c( 881837,895808 ),] #SaPI5
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c( 1545912,1592050 ),] #PVL
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c( 1924777,1959376 ),] #vSAbeta
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c( 2084658,2127720 ),] #beta phage
atBeginning <- Sys.time()
then <- Sys.time()
for (i in 1:length(my_list.interg)){
  for (j in 1:length(good.interg$`genomic position`)){
    if(good.interg$`genomic position`[j] %in% my_list.interg[[i]]$`genomic position`) {
      good.interg[j,i] <- 1
    }else
      good.interg[j,i] <- 0
  }
}
print(Sys.time()-then) # 24.43684 mins

mut.location <- which(colnames(good.interg) == "locus tag")
phylo.interg = good.interg
ref <- which(colnames(phylo.interg) == "ref")
alt <- which(colnames(phylo.interg) == "alt")
atBeginning <- Sys.time()
then <- Sys.time()
for (i in 1:(mut.location -1)){
  for (j in 1:nrow(phylo.interg)){
    if(phylo.interg[j,i] == 1) {
      phylo.interg[j,i] <- phylo.interg[j,alt]
    }else
      phylo.interg[j,i] <- phylo.interg[j,ref]
  }
}
print(Sys.time()-then)

MST.interg <- good.interg[,c(mut.nature:length(good.interg),1:(mut.nature -1))]
phylo.interg <- phylo.interg[,c(mut.nature:length(phylo.interg),1:(mut.nature -1))]

write.table(MST.interg, file = "180326_annnotation_intergenic.txt", sep="\t", row.names = FALSE, quote = FALSE) 


final.phylo <- rbind.fill(phylo.df[c("genomic position","ref", "alt", my_names)], 
                          phylo.interg[c("genomic position","ref", "alt", my_interg)])
outgroup <- final.phylo["ref"]
#colnames(outgroup) <- c("USA300_FPR3757")
colnames(outgroup) <- ReFF_namelabels[genome.tested]
final.phylo <- cbind(final.phylo[1:3],outgroup, final.phylo[4:length(final.phylo)] )
final.annot <- rbind.fill(MST.df, MST.interg)
setwd(paste("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/",ReFF_namelabels[genome.tested],sep=""))
write.table(final.phylo, file = "180326_annnotation_phylo.txt", sep="\t", row.names = FALSE, quote = FALSE) 
write.table(final.annot, file = "180326_annnotation.txt", sep="\t", row.names = FALSE, quote = FALSE) 


getwd()
#phylogeny

phylogeny_name = "phly1.txt"
phyml_name = "phylo_PSCP"
date = "180524"
variants = read.delim(phylogeny_name, as.is=TRUE)

na.col <- which(as.numeric(sapply(variants, function(x) sum(is.na(x)))) >0)
# R sometimes replaces some column names with certain characters, so read in separately 
varHeader = scan(phylogeny_name,what=character(0),nlines=1,sep="\t") 
rmHeader <- varHeader[na.col]
if(length(na.col) >0){
  varHeader <- varHeader[-na.col]
}

if(length(na.col) >0){
  variants <- variants[,-na.col]
}



# Assumes a fasta file representing a single genome
readFastaRef = function(refFile) {
  row = scan(refFile,what=character(0),sep="\n")
  chars = substr(row,1,1)
  base = chars!=">"
  seq = paste(row[base],collapse="")
  return(toupper(unlist(strsplit(seq,""))))
}
# Read the reference genome
setwd("~/Desktop/PSCP_PROJECT/Reference/")
REF_Tested = readFastaRef(paste("~/Desktop/PSCP_PROJECT/Reference/",ReFF_name[genome.tested],sep = ""))
# Number of sites of A, C, G, T
table(REF_Tested) 

# Identify the first column containing sequence data for the samples
outgroup = which(varHeader==ReFF_namelabels[genome.tested])

# Get the genome names
seqLabels = varHeader[outgroup:length(varHeader)]

# Count the number of genomes
nSeq = length(seqLabels)

# Extract base calls & transpose so each row is a genome 
vSites = t(variants[,outgroup:ncol(variants)])

# Ensure all characters are in upper case
vSites = toupper(vSites)

# Replace ambiguous calls
vSites[!(vSites=="A" | vSites=="C" | vSites=="G" | vSites=="T")] = "N"
#variants$Reference..base[!(variants$Reference..base=="A" | variants$Reference..base=="C" | variants$Reference..base=="G" | variants$Reference..base=="T")] = "N"


# Get the total length of the reference genome
lengthRef = length(REF_Tested)
lengthSites = length(vSites[1,])


# Isolate SNPs coordinates
vPos = sapply(variants[1], function(s) {
  as.numeric(unlist(s))
})
typeof(vPos)

# Sanity check
table(variants[,2],REF_Tested[vPos])


# Preparation of nucleotide sequence data into the correct format
# Specify PHYLIP file name
phylipOutfile = paste(date,phyml_name, sep="")

# Output PHYLIP header: sequence number and sequence length
cat(paste(nSeq,lengthRef),sep="\n",file=phylipOutfile)

# For each genome, append to the PHYLIP file
for(i in 1:nSeq) {
  fullLengthRef = REF_Tested
  fullLengthRef[vPos] = vSites[i,]
  fullLengthRefCat = paste0(fullLengthRef,collapse="")
  cat(paste(seqLabels[i]," ",fullLengthRefCat,sep=""),sep="\n",file=phylipOutfile,append=TRUE)
  cat("Output",i,"of",nSeq,"\n")
}

paste("PhyML-3.1_macOS-MountainLion -i ",phylipOutfile," -b 0 -v 0 -c 1 -s BEST --no_memory_check", sep = "")

##### On terminal : 
# PATH=$PATH:/Users/dingd02/Desktop/phyml/PhyML-3.1/
# PhyML-3.1_macOS-MountainLion -i 180125_annnotation_all_ST8.phylip -b 0 -v 0 -c 1 -s BEST --no_memory_check
# With bootstraps: PhyML-3.1_macOS-MountainLion -i 161102_phyML_ha_ac.phylip -b 1000 -v 0 -c 1 -s BEST --no_memory_check 


















