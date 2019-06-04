############## PHYLOGENY #################
rm(list=ls())


# INPUTS
setwd("~/Desktop/R_Scripts/source")
phylogeny_name = "170322_annnotation_deepS_USA300_FPR3757_phylo.txt"
save.name="phylogeny"
date<-"02232018_"

#CODE:
variants = read.delim(phylogeny_name)


# R sometimes replaces some column names with certain characters, so read in separately 
varHeader = scan(phylogeny_name,what=character(0),nlines=1,sep="\t") 

# Assumes a fasta file representing a single genome
readFastaRef = function(refFile) {
  row = scan(refFile,what=character(0),sep="\n") #read the file into a vectors by each new line
  chars = substr(row,1,1) # get the first character from each line(vector)
  base = chars!=">" #result will be True if chars not equal ">",otherwise False
  seq = paste(row[base],collapse="") #Concatenate/combine all vectors.
  return(toupper(unlist(strsplit(seq,"")))) #split the sequence to a character vector by space(anything),
                                            #unlist:to simplifies it to produce a vector which contains all the atomic components which occur in x.
                                            #toupper: translate from lower case to upcase
}
# Read the reference genome

USA300_FPR3757 = readFastaRef("USA300_FPR3757.fasta")
# Number of sites of A, C, G, T
table(USA300_FPR3757) #count the number of each factors level.

# Identify the first column containing sequence data for the samples
varHeader
outgroup = which(varHeader=="USA300_FPR3757")

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
lengthRef = length(USA300_FPR3757) #length of the reference
lengthSites = length(vSites[1,]) #length of first row in vsite?


# Isolate SNPs coordinates
vPos = sapply(variants[1], function(s) {
  as.numeric(unlist(s))
})  #return a list same length as first column and make it single chanracter then numeric 
typeof(vPos) #determines the (R internal) type or storage mode of any object

# Sanity check
table(variants[,2],USA300_FPR3757[vPos])# ref, refposition
?table# what does table do here???


# Preparation of nucleotide sequence data into the correct format
# Specify PHYLIP file name
phylipOutfile = paste(date,save.name,".phylib", sep="")

# Output PHYLIP header: sequence number and sequence length
cat(paste(nSeq,lengthRef),sep="\n",file=phylipOutfile)
?cat
# For each genome, append to the PHYLIP file
for(i in 1:nSeq) {    #from 1 to 183 (183 genomes)
  fullLengthRef = USA300_FPR3757 # reference with 2872769
  fullLengthRef[vPos] = vSites[i,] #replace the the bp from same postion from the genome 
  fullLengthRefCat = paste0(fullLengthRef,collapse="") #combine all bp together
  cat(paste(seqLabels[i]," ",fullLengthRefCat,sep=""),sep="\n",file=phylipOutfile,append=TRUE) #combine title and all bp together
  cat("Output",i,"of",nSeq,"\n") # output status                                                             #append to the file
}



paste("PhyML-3.1_macOS-MountainLion -i ",phylipOutfile," -b 0 -v 0 -c 1 -s BEST --no_memory_check", sep = "")


library(ape)
tree<-read.tree(paste(phylipOutfile,"_phyml_tree.txt", sep = ""))
PatristicDistMatrix<-cophenetic.phylo(tree)
write.table(PatristicDistMatrix, file = paste(phylipOutfile,"_PATRISTIC.txt", sep = ""),sep="\t", quote = FALSE)


##### On terminal : 
# export PATH=$PATH:/Users/copinr01/Richard_copin/Ressources/Programs/PhyML-3.1/
# PhyML-3.1_macOS-MountainLion -i 180125_annnotation_all_ST8.phylip -b 0 -v 0 -c 1 -s BEST --no_memory_check
# With bootstraps: PhyML-3.1_macOS-MountainLion -i 161102_phyML_ha_ac.phylip -b 1000 -v 0 -c 1 -s BEST --no_memory_check 



# Specify FASTA file name for each genome, append to the FASTA file
#for(i in 1:nSeq) {
 # cat(paste(">", seqLabels[i], sep=""),sep="\n",file=fastaOutfile,append=TRUE)
  #seq = paste0(vSites[i,], collapse="")
  #cat(seq, sep="\n", file=fastaOutfile,append=TRUE)
#}


# --------------------------
#  BEAST

# Preparation of nucleotide sequence data into the correct format

# Specify NEXUS file name


nexusOutfile = paste(date,nex_name, sep="")

# Output PHYLIP header: sequence number and sequence length
cat(paste("#NEXUS"),sep="\n",file=nexusOutfile,append=TRUE)
cat(paste(""),sep="\n",file=nexusOutfile,append=TRUE)
cat(paste("Begin DATA;"),sep="\n",file=nexusOutfile,append=TRUE)
cat(paste('\t', "Dimensions ntax=",nSeq, " nchar=", lengthSites, ";", sep = ""),sep="\n",file=nexusOutfile,append=TRUE)
cat(paste('\t', "Format datatype=nucleotide gap=-;",sep = ""),sep="\n",file=nexusOutfile,append=TRUE)
cat(paste('\t', "Matrix",sep = ""),sep="\n",file=nexusOutfile,append=TRUE)
# For each genome, append to the Nexus file
for(i in 1:nSeq) {
  cat(paste(seqLabels[i]),sep="\n",file=nexusOutfile,append=TRUE)
  seq = paste0(vSites[i,], collapse="")
  cat(seq, sep="\n", file=nexusOutfile,append=TRUE)
}

cat(paste('\t', ";",sep = ""),sep="\n",file=nexusOutfile,append=TRUE)
cat(paste("End;"),sep="\n",file=nexusOutfile,append=TRUE)

# --------------------------


