#method1
setwd("/Users/dingd02/Desktop/sequences_for_practice/")
files=list.files(path=".")
nfiles=length(list.files(path="."))
name<-paste("seq",files[1:nfiles],sep="_")
file.rename(files[1:nfiles],name[1:nfiles])
list.files(path=".")

#method2
file.name <- list.files(pattern = ".fastq")
file.name2 <- sapply(file.name, function(s) paste("seq_",s, sep=""))
script = "mv.sh"
cat(paste("#!/bin/bash"),sep="\n",file=script,append=TRUE)
for(i in 1:length(file.name)){
  cat(paste("mv", file.name[i],file.name2[i]),sep="\n",file=script,append=TRUE)
}
getwd()
setwd("/Users/dingd02/Desktop/PSCP_Project/PSCP_BATCH9/")
#method3
file.name <- scan("cp.txt", what = character(), sep = "\t",nlines = 1)
script = "cp.sh"
cat(paste("#!/bin/bash"),sep="\n",file=script,append=TRUE)
for(i in 1:length(file.name)){
  cat(file.name[i],sep="\n",file=script,append=TRUE)
}


#get Q_sub partially
setwd("/Users/dingd02/Desktop/PSCP_Project/Dep_seq/mapping_program/USA300_FPR3757/New")
file.name <- scan("part2.txt", what = character(), sep = "\t",nlines = 1)
Qsub = "Part2_Qsub.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Qsub,append=TRUE)
cat(paste("#$ -S /bin/bash"),sep="\n",file=Qsub,append=TRUE)
cat(paste("#$ -cwd"),sep="\n",file=Qsub,append=TRUE)
cat(paste("#$ -j y"),sep="\n",file=Qsub,append=TRUE)
cat(paste("#$ -M dan.ding@nyumc.org"),sep="\n",file=Qsub,append=TRUE)
cat(paste("#$ -m be"),sep="\n",file=Qsub,append=TRUE)
for(i in 1:length(file.name)){
  cat(file.name[i],sep="\n",file=Qsub,append=TRUE)
}



setwd("/Users/dingd02/Desktop/PSCP_Project/PSCP_batch1")
file.name <- list.files(pattern = ".fastq")
setwd("/Users/dingd02/Desktop/MRSA_MASTER_LIST")
library(readxl)
Master_list_032215 <- read_excel("Master list_032215.xlsx",sheet = "Sheet1")
Master_list_032215_v1 <- read_excel("Master list_032215_v1.xlsx",sheet = "Sheet2")
merged.Master_list_032215 <- merge(Master_list_032215, Master_list_032215_v1, by="V1",all=T)
