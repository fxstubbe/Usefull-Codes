rm(list=ls())
library(stringr)
# JH1
setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/JH1/")

file_o=read.delim("180326_annnotation_test.txt", as.is=TRUE)
file_1=file_o[,12:length(file_o)]
file=read.delim("180326_annnotation_phylo_test.txt", as.is=TRUE)
file1=file[,5:length(file)]
difference.matrix<-apply(file_1,2,function(x)colSums(x!=file_1))
diag(difference.matrix)<-0
write.csv(difference.matrix, "JH1_MATRIX.csv")

JH1<-which(difference.matrix>0 & difference.matrix<=80,arr.ind =T)
JH1.2 <- JH1
JH1.2[,1] <- rownames(JH1.2)
rownames(JH1.2) <- NULL
JH1.2 <- data.frame(JH1.2, stringsAsFactors = F)
JH1.2$col <- as.numeric(JH1.2$col)
str(JH1.2)
JH1.3 <- JH1.2

for(i in 1:nrow(JH1.2)){
  JH1.3$col[i] <- colnames(difference.matrix)[JH1.2$col[i]]
  JH1.3$SNP_number[i] <- difference.matrix[which(rownames(difference.matrix) == JH1.2$row[i]), which(colnames(difference.matrix) == colnames(difference.matrix)[JH1.2$col[i]])]
  JH1.3$SNP_number1[i]<-difference.matrix[JH1[i,1],JH1[i,2]]
  JH1.3$order_min[i]<-pmin(JH1.3$col[i],JH1.3$row[i])
  JH1.3$order_max[i]<-pmax(JH1.3$col[i],JH1.3$row[i])
}

write.csv(JH1.3, "SNP_distance_JH1.csv")


#USA300............................................................
rm(list=ls())


setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/USA300_FPR3757/")

file=read.delim("180326_annnotation.txt", as.is=TRUE)
file1=file[,12:length(file)]
difference.matrix<-apply(file1,2,function(x)colSums(x!=file1))
Shared_snp.matrix<-apply(file1,2,function(x)colSums(x==file1 & x==1))
diag(difference.matrix)<-0
write.csv(difference.matrix, "USA300_matrix.csv")

USA300<-which(difference.matrix>=0 & difference.matrix<=80,arr.ind =T)
USA300.2 <- USA300
USA300.2[,1] <- rownames(USA300.2)
rownames(USA300.2) <- NULL
USA300.2 <- data.frame(USA300.2, stringsAsFactors = F)
USA300.2$col <- as.numeric(USA300.2$col)
str(USA300.2)
USA300.3 <- USA300.2

for(i in 1:nrow(USA300.2)){
USA300.3$col[i] <- colnames(difference.matrix)[USA300.2$col[i]]
USA300.3$SNP_number[i] <- difference.matrix[which(rownames(difference.matrix) == USA300.2$row[i]), which(colnames(difference.matrix) == colnames(difference.matrix)[USA300.2$col[i]])]
USA300.3$SNP_number1[i]<-difference.matrix[USA300[i,1],USA300[i,2]]
USA300.3$order_min[i]<-pmin(USA300.3$col[i],USA300.3$row[i])
USA300.3$order_max[i]<-pmax(USA300.3$col[i],USA300.3$row[i])
}
USA300.4<-USA300.3[(USA300.3$`row`)!=(USA300.3$`col`),]

write.csv(USA300.4, "SNP_distance_USA300.csv")


?apply


#USA500....................................................
rm(list=ls())

setwd("/Users/dingd02/Desktop/PSCP_PROJECT/dep_seq/annotation/USA500_2395/")

file=read.delim("180326_annnotation_phylo.txt", as.is=TRUE)
file1=file[,5:length(file)]
difference.matrix<-apply(file1,2,function(x)colSums(x!=file1))
diag(difference.matrix)<-0
write.csv(difference.matrix, "USA500_matrix.csv")

USA500<-which(difference.matrix>0 & difference.matrix<=80,arr.ind =T)
USA500.2 <- USA500
USA500.2[,1] <- rownames(USA500.2)
rownames(USA500.2) <- NULL
USA500.2 <- data.frame(USA500.2, stringsAsFactors = F)
USA500.2$col <- as.numeric(USA500.2$col)
str(USA500.2)
USA500.3 <- USA500.2

for(i in 1:nrow(USA500.2)){
  USA500.3$col[i] <- colnames(difference.matrix)[USA500.2$col[i]]
  USA500.3$SNP_number[i] <- difference.matrix[which(rownames(difference.matrix) == USA500.2$row[i]), which(colnames(difference.matrix) == colnames(difference.matrix)[USA500.2$col[i]])]
  USA500.3$SNP_number1[i]<-difference.matrix[USA500[i,1],USA500[i,2]]
  USA500.3$order_min[i]<-pmin(USA500.3$col[i],USA500.3$row[i])
  USA500.3$order_max[i]<-pmax(USA500.3$col[i],USA500.3$row[i])
  
}

write.csv(USA500.3, "SNP_distance_USA500.csv")



