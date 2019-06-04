setwd("~/Desktop/R_Scripts/")
rm(list=ls())

# file with mutations for each strains in binary
Mat = read.delim("161026_annnotation_ha_ac_ctl_mst.txt", as.is=TRUE)
varHeader = scan("161026_annnotation_ha_ac_ctl_mst.txt",what=character(0),nlines=1,sep="\t") 
dim(Mat)

# create a matrix of same dimensions
new.Mat <- matrix(1, nrow(Mat), length(Mat))
new.Mat <- data.frame(new.Mat)
colnames(new.Mat) <- varHeader

#implement new matrix with if((A1=B1), 0, 1). if A1 is different than B1 it means 1SNP difference
for(k in 1:length(Mat)){
  for (i in 1:nrow(Mat)) {
    for (j in 1:length(Mat)) {
      if (Mat[i,k] == Mat[i,j]) {
        new.Mat[i,j] <- 0
      } else{
        new.Mat[i,j] <- 1
      }
    }
  }
  assign(paste("result", k, sep = ""), new.Mat)
}

# save temp files
df_list = lapply(ls(pattern = "result[0-9]"), get)
names(df_list) <- ls(pattern = "result[0-9]")

lapply(1:length(df_list), function(i) write.csv(df_list[[i]], 
                                                file = paste0(names(df_list[i]), ".csv"),
                                                row.names = FALSE))

# names of temp files

a <- c("result")
b <- c(1:length(varHeader))
c <- paste(b, "-", sep="")
newnames <- apply(expand.grid(a, c), 1, paste, collapse="")
newnames = sapply(newnames, function(s) {
  s = gsub("-", ".csv", s)
})


#calculate the distance for each genome
mat.distance <- matrix(1, length(varHeader), length(varHeader))
mat.distance <- data.frame(mat.distance)
rownames(mat.distance) <- varHeader
colnames(mat.distance) <- varHeader

for (i in 1:length(newnames)) {  
  d <- read.csv(newnames[i], sep=",")
  for (j in 1:length(newnames)) {
    mat.distance[i,j] <- sum(d[,j])
  }
}

name.av <- c("Average_mutations")

av.mat.dist <- as.data.frame(round(apply(mat.distance, 2, function(x){ mean(x[x !=0])})))
colnames(av.mat.dist) <- name.av
mat.distance <- rbind(mat.distance, t(av.mat.dist))

write.csv(mat.distance, "161027_distance_matrix_ha_ac_ctl.csv")

mat.distance <- matrix(1, length(varHeader), length(varHeader))
mat.distance <- data.frame(mat.distance)
rownames(mat.distance) <- varHeader
colnames(mat.distance) <- varHeader

for (i in 1:length(newnames)) {  
  d <- read.csv(newnames[i], sep=",")
  for (j in 1:length(newnames)) {
    mat.distance[i,j] <- sum(d[,j])/nrow(d)
  }
}



# An ALGORITHM to compute minimum spanning trees for a weighted undirected graph is the Prim’s Algorithm
#install.packages("igraph")
#install.packages("/Users/copinr01/Desktop/igraph", repos = NULL, type="source")
library(igraph)
library(ggplot2)
library(data.table)


links <-read.table("180125_annnotation_all_ST8.phylip_PATRISTIC.txt", row.names = 1)
links <- as.matrix(links, mode = "all")

net <- graph.adjacency(links, weighted=TRUE)
#write.table(as.character(V(net)$name), file = "180125_PATRISTIC_allStrains_namesv4.txt", sep="\t", row.names = FALSE, quote = FALSE) 


nodes <-read.table("180125_PATRISTIC_allStrains_namesv4.txt", header = T, sep = "\t")
str(nodes)




#V(net)$frame.color <- c("dark red", "red", "blue")[nodes$names]
V(net)$frame.color <- "black"
V(net)$shape <- "circle"
V(net)$label <- as.character(nodes$simplified_name)
#V(net)$label <- as.character(V(net)$name)
V(net)$label.cex=.3
V(net)$label.font=1
E(net)$width <- 0.8
V(net)$size <- 1
E(net)$arrow.mode <- 0
E(net)$arrow.size <- .1

mst <- minimum.spanning.tree(net)
deg <- igraph::degree(mst, mode="out")
#deg[7] <- 6
V(mst)$size <- deg/4 +3
V(mst)$size[which(as.character(V(net)$name) == "Ha_K_22")] <- 7
#V(mst)$size[which(as.character(V(net)$name) == "C22")] <- 8
#V(mst)$size[which(as.character(V(net)$name) == "C40")] <- 8
#V(mst)$size[which(as.character(V(net)$name) == "C24")] <- 8
#V(mst)$size[which(as.character(V(net)$name) == "USA300_FPR3757")] <- 5
V(mst)$size[which(as.character(V(net)$name) == "Ha_K_40")] <- 5
V(mst)$size[which(as.character(V(net)$name) == "Ha_K_29")] <- 5
#V(mst)$size[which(as.character(V(net)$name) == "PSCP89_83_GTGTTCTA")] <- 1
V(mst)$size[which(as.character(V(net)$name) == "Ha_K_11")] <- 5
E(mst)$size[which(as.character(V(net)$name) == "Ha_K_46")] <- 5



#lay <-  layout.lgl(mst,maxiter = 25000, area = vcount(mst)^.1, maxdelta = vcount(mst),coolexp = 1000, 
#                   repulserad = 1000 *vcount(mst), cellsize = sqrt(sqrt(1)))

lay <-  layout.lgl(mst,maxiter = 25000, area = vcount(mst)^.1, maxdelta = vcount(mst),coolexp = 250, 
                   repulserad = 1000 *vcount(mst), cellsize = sqrt(sqrt(1)))
lay[100,][1] <-236.0226
lay[100,][2] <-357.9727
lay[100,][1] <- lay[100,][1]
lay[100,][2] <- lay[100,][2]-1 # go down

lay[99,][1] <-236.2082
lay[99,][2] <-358.5039
lay[99,][1] <- lay[99,][1]+0.3 # go right
lay[99,][2] <- lay[99,][2]+0.7 # go up

lay[101,][1] <-236.0634
lay[101,][2] <-358.1635
lay[101,][1] <- lay[101,][1]
lay[101,][2] <- lay[101,][2]-0.5

lay[102,][1] <-236.1725
lay[102,][2] <-358.3132
lay[102,][1] <- lay[102,][1]+0.03
lay[102,][2] <- lay[102,][2]-0.8

lay[103,][1] <-236.1725
lay[103,][2] <-358.3132
lay[103,][1] <- lay[103,][1]
lay[103,][2] <- lay[103,][2]

#tiff('new_plot.tiff', units="in", width=20, height=10, res=300)
par(mar=c(0,0,0,0)+.1)
plot(mst, layout=lay,
     edge.arrow.size=0.1,
     edge.width =1,
     edge.color="black",
     vertex.label.color="black",
     vertex.label=NA,
     #vertex.size=3,
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=F,
     asp = .6
) 

dev.off()



#V(net)$color <- palette()[nodes$neighborhood]
#V(net)$color <- c( "blue","red", "purple")[nodes$neighborhood]
#V(net)$color <- palette()[nodes$Mup_plas]
#write.table(color.table, file = "170115_Chlor_color.txt", sep="\t", quote = FALSE) 
#V(net)$color <- c("blue", "red", "dark red")[nodes$names]
#V(net)$color <- c( "white", "orange" )[nodes$Mup_plas]

setwd("~/Richard_copin/Projects/Shopsin/Hasidic_project/MST/new_figures")

#________Hasidic

V(mst)$color <- c("red","blue", "purple")[nodes$node]

tiff('new_Hasidic_names.tiff', units="in", width=25, height=10, res=300)
par(mar=c(0,0,0,0)+.1)
plot(mst, layout=lay,
     edge.arrow.size=0.1,
     edge.width =1,
     edge.color="black",
     vertex.label.color="black",
     vertex.label=NA,
     #vertex.size=3,
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=F,
     asp = .6
) 
dev.off()



#color.table <- cbind(V(mst), V(mst)$color)
#write.table(color.table, file = "170115_Hasidic.txt", sep="\t", quote = FALSE) 



#________Chlorhexidine

V(mst)$color <- c( "white","purple")[nodes$Mup_plas_chlor]

tiff('new_Chlorhexidine.tiff', units="in", width=20, height=10, res=300)
par(mar=c(0,0,0,0)+.1)
plot(mst, layout=lay,
     edge.arrow.size=0.1,
     edge.width =1,
     edge.color="black",
     vertex.label.color="black",
     vertex.label=NA,
     #vertex.size=3,
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=F,
     asp = .6
) 
dev.off()

#color.table <- cbind(V(mst), V(mst)$color)
#write.table(color.table, file = "170115_Chlor_color.txt", sep="\t", quote = FALSE) 


#________Mupirocin

V(mst)$color <- c("white","orange")[nodes$mup]

tiff('2018_new_Mupirocin.tiff', units="in", width=20, height=10, res=300)
par(mar=c(0,0,0,0)+.1)
plot(mst, layout=lay,
     edge.arrow.size=0.1,
     edge.width =1,
     edge.color="black",
     vertex.label.color="black",
     vertex.label=NA,
     #vertex.size=3,
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=F,
     asp = .6
)  
dev.off()

#color.table <- cbind(V(mst), V(mst)$color)
#write.table(color.table, file = "170115_MupAplas_color.txt", sep="\t", quote = FALSE) 



#________Pyr_Operon

V(mst)$color <- c("white","cyan", "red")[nodes$PyrR_operonMutations]

tiff('new_Pyr_Operon.tiff', units="in", width=20, height=10, res=300)
par(mar=c(0,0,0,0)+.1)
plot(mst, layout=lay,
     edge.arrow.size=0.1,
     edge.width =1,
     edge.color="black",
     vertex.label.color="black",
     vertex.label=NA,
     #vertex.size=3,
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=F,
     asp = .6
) 
dev.off()

#color.table <- cbind(V(mst), V(mst)$color)
#write.table(color.table, file = "170116_pyrOperon_color.txt", sep="\t", quote = FALSE)

#________PyrRGene CarB_PyrC compensatory mutation

V(mst)$color <- c("white","green", "red")[nodes$CarB_PyrC_operonMutations]

tiff('new_CarB_PyrC_mutations.tiff', units="in", width=20, height=10, res=300)
par(mar=c(0,0,0,0)+.1)
plot(mst, layout=lay,
     edge.arrow.size=0.1,
     edge.width =1,
     edge.color="black",
     vertex.label.color="black",
     vertex.label=NA,
     #vertex.size=3,
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=F,
     asp = .6
) 
dev.off()

#color.table <- cbind(V(mst), V(mst)$color)
#write.table(color.table, file = "170116_pyrRgene_color.txt", sep="\t", quote = FALSE)

#________PyrRGene

V(mst)$color <- c("red", "white")[nodes$PyrR_geneMutations]
#color.table <- cbind(V(mst), V(mst)$color)
#write.table(color.table, file = "170116_pyrRgene_color.txt", sep="\t", quote = FALSE)


#________Selected_strain_for_cytotox

V(mst)$color <- c("white","green")[nodes$Cytotox_testing]

tiff('tested_Cytotox.tiff', units="in", width=20, height=10, res=300)
par(mar=c(0,0,0,0)+.1)
plot(mst, layout=lay,
     edge.arrow.size=0.1,
     edge.width =1,
     edge.color="black",
     vertex.label.color="black",
     vertex.label=NA,
     #vertex.size=3,
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=F,
     asp = .6
) 
dev.off()

#color.table <- cbind(V(mst), V(mst)$color)
#write.table(color.table, file = "170116_cytotoxSelection_color.txt", sep="\t", quote = FALSE)


#________high_cytotox

V(mst)$color <- c("brown","white")[nodes$High_cytotox]

tiff('high_Cytotox.tiff', units="in", width=20, height=10, res=300)
par(mar=c(0,0,0,0)+.1)
plot(mst, layout=lay,
     edge.arrow.size=0.1,
     edge.width =1,
     edge.color="black",
     vertex.label.color="black",
     vertex.label=NA,
     #vertex.size=3,
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=F,
     asp = .6
) 
dev.off()

#color.table <- cbind(V(mst), V(mst)$color)
#write.table(color.table, file = "170116_cytotoxSelection_color.txt", sep="\t", quote = FALSE)


#________Phage11

V(mst)$color <- c("white","white", "blue")[nodes$Ph11]

tiff('new_Phage11_names.tiff', units="in", width=20, height=10, res=300)
par(mar=c(0,0,0,0)+.1)
plot(mst, layout=lay,
     edge.arrow.size=0.1,
     edge.width =1,
     edge.color="black",
     vertex.label.color="black",
     #vertex.label=NA,
     #vertex.size=3,
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=F,
     asp = .6
) 
dev.off()

#color.table <- cbind(V(mst), V(mst)$color)
#write.table(color.table, file = "170115_Hasidic_Island_color.txt", sep="\t", quote = FALSE) 


#________Tested for skin infection

V(mst)$color <- c("white","brown")[nodes$tested_skin_infection]

tiff('tested_Abscess_size.tiff', units="in", width=20, height=10, res=300)
par(mar=c(0,0,0,0)+.1)
plot(mst, layout=lay,
     edge.arrow.size=0.1,
     edge.width =1,
     edge.color="black",
     vertex.label.color="black",
     vertex.label=NA,
     #vertex.size=3,
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=F,
     asp = .6
)  
dev.off()

#color.table <- cbind(V(mst), V(mst)$color)
#write.table(color.table, file = "170115_Hasidic_Island_color.txt", sep="\t", quote = FALSE) 



#________Abscess size

V(mst)$color <- c("pink","white", "purple")[nodes$big_abscess]

tiff('new_Abscess_size.tiff', units="in", width=20, height=10, res=300)
par(mar=c(0,0,0,0)+.1)
plot(mst, layout=lay,
     edge.arrow.size=0.1,
     edge.width =1,
     edge.color="black",
     vertex.label.color="black",
     vertex.label=NA,
     #vertex.size=3,
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=F,
     asp = .6
) 
dev.off()

#color.table <- cbind(V(mst), V(mst)$color)
#write.table(color.table, file = "170115_Hasidic_Island_color.txt", sep="\t", quote = FALSE) 


#________remaining strains

V(mst)$color <- c("white","green")[nodes$New_strain]
color.table <- cbind(V(mst), V(mst)$color)
tiff('180126_newStrain_color.tiff', units="in", width=20, height=10, res=300)
par(mar=c(0,0,0,0)+.1)
plot(mst, layout=lay,
     edge.arrow.size=0.1,
     edge.width =1,
     edge.color="black",
     vertex.label.color="black",
     vertex.label=NA,
     vertex.size=3,
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=F,
     asp = .6
) 
dev.off()
write.table(color.table, file = "180126_newStrain_color.txt", sep="\t", quote = FALSE)

#________spa

V(mst)$color <- c( "yellow","purple")[nodes$spa]

tiff('spa.tiff', units="in", width=20, height=10, res=300)
par(mar=c(0,0,0,0)+.1)
plot(mst, layout=lay,
     edge.arrow.size=0.1,
     edge.width =1,
     edge.color="black",
     vertex.label.color="black",
     #vertex.label=NA,
     #vertex.size=3,
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=F,
     asp = .6
) 
dev.off()

#color.table <- cbind(V(mst), V(mst)$color)
#write.table(color.table, file = "170115_Chlor_color.txt", sep="\t", quote = FALSE) 





plot(mst, layout=lay,
     edge.arrow.size=0.1, 
     vertex.label.color="black",
     vertex.label=NA,
     edge.color="black",
     rescale=TRUE,
     #vertex.label.dist=0.2,
     #vertex.frame.color = NA
     add=FALSE
     ) 



# To consider choosing another layout
layouts <- grep("^layout\\.", ls("package:igraph"), value=TRUE)
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama", layouts)]
par(mfrow=c(3,3))

for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(mst)) 
  plot(mst, edge.arrow.mode=0, layout=l, main=layout,vertex.label=NA) }
dev.off()












# An ALGORITHM to compute minimum spanning trees for a weighted undirected graph is the Prim’s Algorithm
library(igraph)
library(ggplot2)
op <- par(oma=c(5,7,1,1))
par(op)

G <- graph.adjacency(as.matrix(mat.distance), weighted=TRUE)

a <- V(G)$name
#new.a = sapply(a, function(s) {
#s = gsub("\\_.*","",s)
#s = gsub("HasidicKid", "Ha", s)
#s = gsub("HasidicAdult", "HaAd", s)
#})

# Simplify
net <- simplify(G, remove.multiple = F, remove.loops = T)
summary(net)

## Some graphical parameters
E(net)$width <- 0.05
# E(net)$width <- E(net)$weight*20
V(net)$label <- a
V(net)$shape <- "circle"
V(net)$color <- "purple"
V(net)$frame.color <- "black"
V(net)$size <- 8
V(net)$label.cex <- 0.4
E(net)$arrow.mode <- 0
E(net)$arrow.size <- .001
E(net)$edge.color <- "red"
## MST and plot
# lay <- layout.kamada.kawai(mst)
# lay <- layout.circle(mst)
# lay <-  layout.random(mst)
mst <- minimum.spanning.tree(net)
lay <-  layout.fruchterman.reingold(mst)
# lay <-  layout.sphere(mst)

plot(net, layout=lay, vertex.label.color="white")
E(net)$width <- 2
mst <- minimum.spanning.tree(net)
par(new=T)
plot(mst, layout=lay, 
     # edge.arrow.size=1, 
     vertex.label.color="white",
     edge.color="black")
# edge.curved=.2) 
par(new=F)

plot(mst, layout=lay, 
     # edge.arrow.size=1, 
     vertex.label.color="white",
     edge.color="black") 
