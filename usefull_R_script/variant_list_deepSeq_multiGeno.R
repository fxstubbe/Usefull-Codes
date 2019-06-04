rm(list=ls())

setwd("~/Richard_copin/Projects/Shopsin/Deep_sequencing/")

ReFF_name = scan("ReFF_names.txt",what=character(0),nlines=1,sep="\t")
ReFF_namelabels = sapply(ReFF_name, function(s) {
  s = gsub(".fasta", "", s)
})


### To edit
# 170405
# USA300_FPR3757
# ReFF_namelabels[9]



setwd("/Users/dingd02/Desktop/annovar1/USA300_FPR3757/Annovation/exonic")



library(data.table)
library(plyr)

my_files <- list.files(pattern = "\\.exonic.txt")
my_names1 <- gsub(".ex.*","",my_files)
my_names <- gsub("merge_","",my_names1)
my_names
my_name_allFreq <- gsub(".ex.*","_allele_Freq",my_files)
my_name_allFreq <- gsub("merge_","",my_name_allFreq)
my_name_allFreq

#my_names <- sapply(my_names, function(s){ 
#  as.character(paste(unlist(strsplit(s, split='_', fixed=TRUE))[c(1,4)], collapse="_"))
# })



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




#for(s in 1:length(my_list)){
#  TOrm <- which(my_list[[s]]$`locus tag` == "exon1")
#  my_list[[s]] <- my_list[[s]][-c(TOrm),]
#}

total.df = data.table::rbindlist(my_list)

gen.mat <- matrix(data = NA, nrow = length(total.df$`genomic position`), ncol = length(my_list))
colnames(gen.mat) <- my_names
total.df <- as.data.frame(cbind(gen.mat, total.df))

for (i in 1:length(my_list)){
  for (j in 1:length(total.df$`genomic position`)){
    if(total.df$`genomic position`[j] %in% my_list[[i]]$`genomic position`) {
      total.df[j,i] <- 1
    }else
      total.df[j,i] <- 0
  }
}

total.df2 = data.table::rbindlist(my_list)

test <- matrix(data = NA, nrow = length(total.df2$`genomic position`), ncol = length(my_list))
colnames(test) <- my_name_allFreq
total.df2 <- as.data.frame(cbind(test, total.df2))

for (i in 1:length(my_list)){
  for (j in 1:length(total.df$`genomic position`)){
    if(total.df2$`genomic position`[j] %in% my_list[[i]]$`genomic position`) {
      total.df2[j,i] <- my_list[[i]]$`allelic frequency`[which(my_list[[i]]$`genomic position`== total.df2$`genomic position`[j])]
    }else
      total.df2[j,i] <- 0
  }
}


total.df.freq <- cbind(total.df, total.df2[,c(1:length(my_list))])






noDup.df <- as.data.frame(subset.data.frame(total.df.freq, !duplicated(total.df$`genomic position`)))
#good.qul <- noDup.df
good.qul <- noDup.df[!(noDup.df$`genomic position`) %between% c(881837,895808),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(1545912,1592050),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(2084658,2127720),]

mut.nature <- which(colnames(good.qul) == "mutation nature")
phylo.df = good.qul
ref <- which(colnames(phylo.df) == "ref")
alt <- which(colnames(phylo.df) == "alt")

for (i in 1:(mut.nature -1)){
  for (j in 1:nrow(phylo.df)){
    if(phylo.df[j,i] == 1) {
      phylo.df[j,i] <- phylo.df[j,alt]
    }else
      phylo.df[j,i] <- phylo.df[j,ref]
  }
}

MST.df <- good.qul[,c(mut.nature:length(good.qul),1:(mut.nature -1))]
phylo.df <- phylo.df[,c(mut.nature:length(phylo.df),1:(mut.nature -1))]


write.table(MST.df, file = "170405_annnotation_deepS_USA300_FPR3757_exonic.txt", sep="\t", row.names = FALSE, quote = FALSE) 



#Intergenic


setwd("/Users/dingd02/Desktop/annovar1/USA300_FPR3757/Annovation/intergenic")


my_files.interg <- list.files(pattern = "\\.intergenic")
my_interg1 <- gsub(".int.*","",my_files.interg)
my_interg <- gsub("merge_","",my_interg1)
my_interg
my_interg_allFreq <- gsub(".int.*","_allele_Freq",my_files.interg)
my_interg_allFreq <- gsub("merge_","",my_interg_allFreq)
my_interg_allFreq



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

for (i in 1:length(my_list.interg)){
  for (j in 1:length(total.interg$`genomic position`)){
    if(total.interg$`genomic position`[j] %in% my_list.interg[[i]]$`genomic position`) {
      total.interg[j,i] <- 1
    }else
      total.interg[j,i] <- 0
  }
}

total.interg2 = data.table::rbindlist(my_list.interg)

x <- matrix(data = NA, nrow = length(total.interg2$`genomic position`), ncol = length(my_list.interg))
colnames(x) <- my_interg_allFreq
total.interg2 <- as.data.frame(cbind(x, total.interg2))

for (i in 1:length(my_list.interg)){
  for (j in 1:length(total.interg$`genomic position`)){
    if(total.interg2$`genomic position`[j] %in% my_list.interg[[i]]$`genomic position`) {
      total.interg2[j,i] <- my_list.interg[[i]]$`allelic frequency`[which(my_list.interg[[i]]$`genomic position`== total.interg2$`genomic position`[j])]
    }else
      total.interg2[j,i] <- 0
  }
}


total.interg.freq <- cbind(total.interg, total.interg2[,c(1:length(my_list.interg))])



noDup.interg <- as.data.frame(subset.data.frame(total.interg.freq , !duplicated(total.interg$`genomic position`)))

#good.interg <- noDup.interg
good.interg <- noDup.interg[!(noDup.interg$`genomic position`) %between% c(881837,895808),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(1545912,1592050),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(2084658,2127720),]


mut.location <- which(colnames(good.interg) == "locus tag")
phylo.interg = good.interg
ref <- which(colnames(phylo.interg) == "ref")
alt <- which(colnames(phylo.interg) == "alt")

for (i in 1:(mut.location -1)){
  for (j in 1:nrow(phylo.interg)){
    if(phylo.interg[j,i] == 1) {
      phylo.interg[j,i] <- phylo.interg[j,alt]
    }else
      phylo.interg[j,i] <- phylo.interg[j,ref]
  }
}

MST.interg <- good.interg[,c(mut.nature:length(good.interg),1:(mut.nature -1))]
phylo.interg <- phylo.interg[,c(mut.nature:length(phylo.interg),1:(mut.nature -1))]

write.table(MST.interg, file = "170405_annnotation_deepS_USA300_FPR3757_intergenic.txt", sep="\t", row.names = FALSE, quote = FALSE) 


final.phylo <- rbind.fill(phylo.df[c("genomic position","ref", "alt", my_names)], 
                          phylo.interg[c("genomic position","ref", "alt", my_interg)])
outgroup <- final.phylo["ref"]
colnames(outgroup) <- c("USA300_FPR3757")
final.phylo <- cbind(final.phylo[1:3],outgroup, final.phylo[4:length(final.phylo)] )

final.annot <- rbind.fill(MST.df, MST.interg)





setwd("/Users/dingd02/Desktop/annovar1/USA300_FPR3757/Annovation/")
write.table(final.phylo, file = "170322_annnotation_deepS_USA300_FPR3757_phylo.txt", sep="\t", row.names = FALSE, quote = FALSE) 
write.table(final.annot, file = "170322_annnotation_deepS_USA300_FPR3757.txt", sep="\t", row.names = FALSE, quote = FALSE) 

