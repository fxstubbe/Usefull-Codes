rm(list=ls())

setwd("~/Richard_copin/Projects/Shopsin/Hasidic_project/")

ReFF_name = scan("../Deep_sequencing/ReFF_names.txt",what=character(0),nlines=1,sep="\t")
ReFF_namelabels = sapply(ReFF_name, function(s) {
  s = gsub(".fasta", "", s)
})
genome.tested <-  which(grepl("USA300_FPR3757", ReFF_namelabels) == T)

setwd("/Users/dingd02/Desktop/annovar1/exonic")

library(data.table)
library(plyr)

my_files <- list.files(pattern = "\\.exonic.txt")
empty <- list()
for (h in 1:length(my_files)) {
  info = file.info(my_files[h])
  empty[[h]] = rownames(info[info$size == 0, ])
}
names(empty) <- my_files
empty <- empty[lapply(empty, length)>0]
my_files <- my_files[- which(my_files %in% names(empty))]

my_names1 <- gsub(".ex.*","",my_files)
my_names <- gsub("merge_","",my_names1)
my_names

my_list <- lapply(my_files, function(i){
  fread( i, sep = "\t", header=F, data.table = F)
})
names(my_list) <- my_names

column.names <- c("mutation nature", "product name", "locus tag", 
                  "nucleotide change","protein change", "genomic position", 
                  "ref", "alt","quality score", "depth of coverage", 
                  "allelic frequency")
my_list <- lapply(my_list, setNames, column.names)

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

noDup.df <- as.data.frame(subset.data.frame(total.df, !duplicated(total.df$`genomic position`)))

good.qul <- noDup.df[!(noDup.df$`genomic position`) %between% c(881837,895808),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(1545912,1592050),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(2084658,2127720),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(34173,57914),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(57915,88900),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(546482,546753),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(1632109,1642734),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(110377,112806),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(128641,129447),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(2811366,2818182),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(347771,352173),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(448388,451806),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(611146,622360),]
good.qul <- good.qul[!(good.qul$`genomic position`) %between% c(1456810,1488076),]



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


write.table(MST.df, file = "180125_annnotation_all_ST8_exonic.txt", sep="\t", row.names = FALSE, quote = FALSE) 



#Intergenic


setwd(paste("~/Richard_copin/Projects/Shopsin/Hasidic_project/Annotation/",ReFF_namelabels[9],"/intergenic/", sep = ""))

my_files.interg <- list.files(pattern = "\\.intergenic")
empty <- list()
for (h in 1:length(my_files.interg.interg)) {
  info = file.info(my_files.interg.interg[h])
  empty[[h]] = rownames(info[info$size == 0, ])
}
names(empty) <- my_files.interg
empty <- empty[lapply(empty, length)>0]
my_files.interg <- my_files.interg[- which(my_files.interg %in% names(empty))]
my_interg1 <- gsub(".int.*","",my_files.interg)
my_interg <- gsub("merge_","",my_interg1)
my_interg


my_list.interg <- lapply(my_files.interg, function(i){
  fread( i, sep = "\t", header=F, data.table = F)
})
names(my_list.interg) <- my_interg

column.interg <- c("mutation location", "gene environment","genomic position", 
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

noDup.interg <- as.data.frame(subset.data.frame(total.interg, !duplicated(total.interg$`genomic position`)))


good.interg <- noDup.interg[!(noDup.interg$`genomic position`) %between% c(881837,895808),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(1545912,1592050),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(2084658,2127720),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(34173,57914),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(57915,88900),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(546482,546753),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(1632109,1642734),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(110377,112806),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(128641,129447),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(2811366,2818182),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(347771,352173),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(448388,451806),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(611146,622360),]
good.interg <- good.interg[!(good.interg$`genomic position`) %between% c(1456810,1488076),]


mut.location <- which(colnames(good.interg) == "mutation location")
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

write.table(MST.interg, file = "180125_annnotation_all_ST8_intergenic.txt", sep="\t", row.names = FALSE, quote = FALSE) 


final.phylo <- rbind.fill(phylo.df[c("genomic position","ref", "alt", my_names)], 
                          phylo.interg[c("genomic position","ref", "alt", my_interg)])
outgroup <- final.phylo["ref"]
colnames(outgroup) <- c("USA300_FPR3757")
final.phylo <- cbind(final.phylo[1:3],outgroup, final.phylo[4:length(final.phylo)] )

final.annot <- rbind.fill(MST.df, MST.interg)


setwd("~/Richard_copin/Projects/Shopsin/Hasidic_project/Annotation/USA300_FPR3757/")
write.table(final.phylo, file = "180125_annnotation_all_ST8_phylo.txt", sep="\t", row.names = FALSE, quote = FALSE) 
write.table(final.annot, file = "180125_annnotation_all_ST8.txt", sep="\t", row.names = FALSE, quote = FALSE) 




