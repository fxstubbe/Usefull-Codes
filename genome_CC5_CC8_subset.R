###
# Plotting % of the palsmid overtime
###

#---------------------------------------------------------------------------------------------------------------------------------------#
# 0) Packages
#---------------------------------------------#

library("tidyverse")
library("purr")
library("lubridate")

#---------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------#
# 1.1) Paths & Inputs CC5
#---------------------------------------------#

#Get the CC5 matrix
master_CC5 <- read_delim("/Users/stubbf02/Fx_Stubbe/ressources/genomes/staph/Phylogeny/Matrix/CC5_snp_matrix_nodup.txt", delim = "\t")

#Get the CC5 features
features_CC5 <- read_delim("/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/pBSRC1_Mapping/CC5/all/tables/summary/iTOL_feature_nodup.txt", 
                           delim = "\t")
features_CC5$isolate_key <- features_CC5 %>% pull(isolate_key) %>% gsub(" ", "_", .)


# 1.2) Paths & Inputs CC8
#---------------------------------------------#

master_CC8 <- read_csv("/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/SPREADSHEETS/Master_CC8_NYU.csv")

#Get the CC8 features
setwd("/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/pBSRC1_Mapping/CC8/nodup_only/Tables/Summary/")
features_CC8 <- read_delim("/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/pBSRC1_Mapping/CC8/nodup_only/Tables/Summary/summary_index.txt", 
                       delim = "\t", col_names = c("isolate_key", "pBSRC1", "start", "stop", "length", "cov", "tot_cov"))
features_CC8 <- features_CC8 %>% select(isolate_key, tot_cov) %>% distinct(.) 

#---------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------#
# 2) Wrangling CC5
#---------------------------------------------#

##parse dates from isolate key
features_CC5 <- features_CC5 %>% mutate(
  date = str_extract(features_CC5$isolate_key, "\\_.*") %>% str_remove(., "_") %>% ymd(.) )%>%  drop_na()
features_CC5$cluster <- "CC5"


# 2) Wrangling CC8
#---------------------------------------------#

##parse dates from isolate key (features)
features_CC8 <- features_CC8 %>% mutate(
  date = str_extract(features_CC8$isolate_key, "\\_.*") %>% str_remove(., "_") %>% ymd(.) ) %>%  drop_na()

##get the religion from master & add it to features
features_CC8 <- master_CC8 %>% select(isolate_key, Religion, CITY) %>% full_join(features, by = "isolate_key")
features_CC8$cluster <- "CC8"

##Subset features by date
#sub_features <- features %>% filter(date >= as.Date("2017-03-01")) 

#---------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------#
# 3.1) MATRIX - CC5 < "2017-03-01"
#---------------------------------------------#

subset_CC5 <- features_CC5 %>% filter(date < as.Date("2017-03-01")) %>% select(isolate_key, pBSRC1,  MupA, QuacA) 

#Get the tree matrix
bin1 <- as.tibble( cbind( master[1:4],master_CC5[ which( colnames(master) %in% subset_CC5$isolate_key ) ] ) )

#Get the iTOL plot for the matrix
subset_CC5$isolate_key <- subset_CC5 %>% pull(isolate_key) %>% gsub("_", " ", .)

#write files          
write_delim(bin1, "/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/pBSRC1_Mapping/CC5/bin2/CC5_bin1_subset.txt", delim = "\t")
write_delim(subset_CC5, "/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/pBSRC1_Mapping/CC5/bin2/CC5_bin1_iTOL_feature.txt", delim = "\t")

# 3.2) MATRIX - CC5 >= 2017-03-01
#---------------------------------------------#

subset_CC5 <- features_CC5 %>% filter(date >= as.Date("2017-03-01")) %>% select(isolate_key, pBSRC1,  MupA, QuacA) 

#Get the tree matrix
bin2 <- as.tibble( cbind( master[1:4],master_CC5[ which( colnames(master) %in% subset_CC5$isolate_key ) ] ) )

#Get the iTOL plot for the matrix
subset_CC5$isolate_key <- subset_CC5 %>% pull(isolate_key) %>% gsub("_", " ", .)

#write files          
write_delim(bin2, "/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/pBSRC1_Mapping/CC5/bin2/CC5_bin2_subset.txt", delim = "\t")
write_delim(subset_CC5, "/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/pBSRC1_Mapping/CC5/bin2/CC5_bin2_iTOL_feature.txt", delim = "\t")


#---------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------#
# 3.1) Plot - Genome Distribution
#---------------------------------------------#

features_CC8 %>% full_join(features_CC5) %>% select("isolate_key", "date", "cluster") %>% drop_na() %>%
  mutate(year = year(date), month = month(date, label = T)) %>% 
  group_by(year, month, cluster) %>% tally() %>%
  ggplot(aes(x = month, y = n, group = cluster)) +
  geom_line(aes(linetype=cluster))+
  geom_point(aes(shape = cluster))+
  #scale_color_manual(values=c("#00BFC4","#F8766D")) +
  facet_wrap(~year, ncol = 3) +
  labs(title = "Nodup genomes",
       subtitle = "Genomes distribution plotted by year",
       y = "#genomes",
       x = "") +
  ylim(0,15) + 
  theme_bw(base_size = 15) +
  theme(legend.position="bottom")


# 3.2) Plot 
#---------------------------------------------#

## Plot the average tot_coverage by month over 3 years
features %>% drop_na() %>%
  mutate(year = year(date), month = month(date, label = T)) %>% 
  group_by(year, month) %>% 
  summarise(
    avg_cov = mean(tot_cov, na.rm = TRUE),
    n = n()) %>% 
  ggplot(aes(x = month, y = avg_cov, group = 1)) +
  #geom_bar(stat = "identity", fill = "darkorchid4") +
  geom_line( )+
  geom_point() + 
  facet_wrap(~year, ncol = 3) +
  labs(title = "Nodup CC8 genomes",
       subtitle = "average pBSRC1 coverage",
       y = "% pBSRC1 average coverage",
       x = "") +
  ylim(0,100) +
  theme_bw(base_size = 15)

## Plot the average tot_coverage by month over 3 years (with religion info)
features %>% drop_na() %>%
  mutate(year = year(date), month = month(date, label = T)) %>% 
  group_by(year, month, Religion) %>% 
  summarise(
    avg_cov = mean(tot_cov, na.rm = TRUE),
    n = n()) %>% 
  ggplot(aes(x = month, y = avg_cov, group = Religion)) +
  geom_line(aes(linetype=Religion, color = Religion))+
  geom_point(aes(color = Religion))+
  scale_color_manual(values=c("#00BFC4","#F8766D")) +
  facet_wrap(~year, ncol = 3) +
  labs(title = "Nodup CC8 genomes",
       subtitle = "average pBSRC1 coverage grouped by Religion",
       y = "% pBSRC1 average coverage",
       x = "") +
  ylim(0,100) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom", legend.title = element_blank())

## Plot the average tot_coverage by month over 3 years (with religion info)
features <- features %>% mutate(CITY = replace(CITY, CITY != "Brooklyn", "Other"))

master_group <- features %>% drop_na() %>%
  mutate(year = year(date), month = month(date, label = T)) %>% group_by(year, month, Religion, CITY) %>% tally() %>% mutate(prevalence = 0)
master_group

master_group2 <- features %>% drop_na() %>%
  mutate(year = year(date), month = month(date, label = T)) %>% group_by(year, month, Religion) %>% tally()
master_group2

for(row in 1:nrow(master_group) ){
  temp <- master_group[row,]
  tot <- master_group2 %>% filter(year == temp$year[1] & month == temp$month[1] & Religion == temp$Religion[1])
  master_group$prevalence[row] <- (temp$n/tot$n)*100
}

master_group %>% filter(CITY == "Brooklyn") %>%
  ggplot(aes(x = month, y = prevalence, group = Religion)) +
  geom_line(aes(linetype=Religion, color = Religion))+
  geom_point(aes(color = Religion))+
  scale_color_manual(values=c("#00BFC4","#F8766D")) +
  facet_wrap(~year, ncol = 3) +
  labs(title = "Nodup CC8 genomes",
       subtitle = "prevalence of Brooklyn zip code",
       y = "Prevalence (%)",
       x = "") +
  ylim(0,100) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none", legend.title = element_blank())
