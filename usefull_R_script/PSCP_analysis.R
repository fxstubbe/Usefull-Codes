setwd("Desktop/")
getwd()
#zip code
fm<-read.csv("pscp_g_12.csv")
?read.csv




library(zipcode)
library(ggmap)
library(plyr)
library(ggplot2)
?clean.zipcodes
fm1<-fm$'Isolate_Key'
fm1 = subset(fm, select = c("Isolate_Key","zip_code") )
fm1$zip<- clean.zipcodes(fm1$zip)
data(zipcode)
fm1<- merge(fm1,zipcode, by.x='zip', by.y='zip')
write.csv(fm1,"upload.csv")





density<- ddply(fm1, .(city), "nrow")
names(density)[2] <- "count"
FM<- merge(fm1, density)

duplicated(FM$city)
FM[duplicated(FM$city),]
unique(FM[duplicated(FM$city),])
FM<-FM[!duplicated(FM$city),]


map<-get_map(location='New York City', zoom=7, maptype='roadmap')
map<-get_map(location=c(lon = -74, lat = 40.8), zoom=10, maptype='roadmap')

ggmap(map)+geom_point(aes(x=longitude, y=latitude, size=(count)), data=FM, alpha=.5)
?get_map


#shared ward

getwd()
rm(list=ls())
setwd("/Users/dingd02/Desktop")
Clade<-"A"
fm<-read.csv("transfer_test_clade.csv")
#fm<-fm[(fm$`DEPARTMENT_NAME`)=="TH 9 PICU" |(fm$`DEPARTMENT_NAME`)=="TH 17 EAST" ,]
fm<-fm[(fm$`clade`)==Clade,]

fm$pt_id=as.character(fm$pt_id)

fm$transfer_in<-as.Date(fm$transfer_in, format = "%m/%d/%Y")
fm$transfer_out=as.Date(fm$transfer_out, format = "%m/%d/%Y")
start_time<-min(fm$transfer_in)
#fm<-fm[(fm$`transfer_in`)>"2016-01-01",]
#fm<-fm[(fm$`transfer_in`)<"2016-07-01",]

#fm1<-read.csv("genome_test_clade.csv")
fm1<-read.csv("pt_cultures_all1.csv")
fm1<-fm1[(fm1$`clade`)==Clade,]
#fm1 = subset(fm1, select = c("pt_id","Date_Of_Culture","clade") )

fm1$collection_date=as.Date(fm1$collection_date, format = "%m/%d/%Y")
#fm1<-fm1[(fm1$`collection_date`)>"2016-01-01",]
fm1$pt_id=as.character(fm1$pt_id)
library(ggplot2)
color1=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#c51b7d","#003c30","#67001f")
ggplot()+ geom_segment(data=fm, aes(x=transfer_in, xend=transfer_out, y=pt_id, yend=pt_id,colour=factor(DEPARTMENT_NAME)), size=10)+xlab("Duration")+
  geom_point(data =fm1, aes(x = collection_date, y=pt_id,shape=factor(MRSA_Result1)),size=2)+ggtitle("Ward-Centric(17 EAST)")+ylab("Patient ID")+scale_color_manual(values=color1)









library(ggplot2)
ggplot(fm, aes(colour=DEPARTMENT_NAME))+ geom_segment(aes(x=transfer_in, xend=transfer_out, y=pt_id, yend=pt_id), size=10)+xlab("Duration")+
  geom_point(data =fm1, aes(x = collection_date, y=pt_id,shape=factor(MRSA_Result1)),colour="black",size=3)+ggtitle("Ward-Centric(9 PICU)")+ylab("Patient ID")+scale_color_brewer(palette="Paired")


library(ggplot2)
ggplot()+ geom_segment(data=fm,aes(x=transfer_in, xend=transfer_out, y=pt_id, yend=pt_id,colour=factor(DEPARTMENT_NAME)), size=10)+xlab("Duration")+
  geom_point(data =fm1, aes(x = collection_date, y=pt_id,shape=factor(MRSA_Result1)),colour="black",size=3)+ggtitle("Ward-Centric(9 PICU)")+ylab("Patient ID")+scale_color_brewer(palette="Paired")


















library(randomcoloR)
n <- 12
palette <- distinctColorPalette(n)
color =c("#fa9fb5","#c51685","#e7e1ef","#c994c7","#dd1c77","#2ca25f","#99d8c9","#8856a7","#9ebcda","#43a2ca","#e34a33","#fdbb84","#ffeda0","#636363")


pie(color)









?plot

?geom_dotplot()
?xlab
?ggplot
