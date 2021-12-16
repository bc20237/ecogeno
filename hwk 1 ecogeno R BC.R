#bring in the data table downloaded from QIIME viewer (taxonomy.qzv)
setwd("/Users/blairchristensen/Documents/GitHub/ecogeno")

# Read in list of positions
tl2 <- read.table("taxlevel2.txt",sep=",", header=T)
#install some stuff I didnt have bc I havent used R in like 4 years
#install.packages("ggplot2")
#install.packages("devtools")
library(ggbiplot)
library(devtools)
install_github("vqv/ggbiplot")
library(tidyverse)
library(ggplot2)
#changed the name to something shorter albeit less descriptive
tl2<-taxlevel2
require(plyr)
#compute total abundance and add to dataset
tl2<-ddply(tl2, "index", transform, total_abundance = sum(Spirochaetes+Bacteroidetes+Proteobacteria+unknownA+Fusobacteria+unknownB+Tenericutes+Firmicutes+Actinobacteria+Cyanobacteria+TM6+GN02+Chlamydiae+Thermi+Fibrobacteres+H.178+ZB3+Chlorobi+Lentisphaerae+Planctomycetes+Nitrospirae+WWE1+OP8+Verrucomicrobia+SAR406+FBP+Chloroflexi+Acidobacteria+WPS.2+Gemmatimonadetes+TM7+SBR1093+WS3+OP3))
#compute props for top 10 groupings and add to dataset
#idk if this was necessary but in the documentation for qiime i couldnt figure out if the values for each were relative or absolute
tl2<-ddply(tl2, "index", transform, Spirochaetesp=(Spirochaetes/total_abundance), Bacteroidetesp=(Bacteroidetes/total_abundance),Proteobacteriap=(Proteobacteria/total_abundance),unknownAp=(unknownA/total_abundance),Fusobacteriap=(Fusobacteria/total_abundance),unknownBp=(unknownB/total_abundance),Tenericutesp=(Tenericutes/total_abundance),Firmicutesp=(Firmicutes/total_abundance),Actinobacteriap=(Actinobacteria/total_abundance),Cyanobacteriap=(Cyanobacteria/total_abundance))
#check out that it worked and also get an idea of what the spread of the data is
summary(tl2)
#made separate samples of sick v healthy
healthysamp<-tl2[1:67,]
sicksamp<-tl2[68:85,]
#summary for all data
summary(tl2)
#descriptive stats for healthy
summary(healthysamp)
#descriptive stats for sick 
summary(sicksamp)
#pca for the proportional abundance of top 10 groups
tl2.pca<-prcomp(tl2[,c(44:53)], center=TRUE, scale.=TRUE)
#look at the table output
summary(tl2.pca)
str(tl2.pca)
#plot the PCAs, by health status site and site/health status just to check it out
#health status
ggbiplot(t12.pca, obs.scale = 1, var.scale = 1, 
         groups=tl2[,42],
   ellipse = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
#site
ggbiplot(tl2.pca, obs.scale = 1, var.scale = 1, 
         groups=tl2[,41],
         ellipse = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
#this one is for site health status and is the one i was most interested in for my analysis
ggbiplot(tl2.pca, obs.scale = 1, var.scale = 1, 
         groups=tl2[,42],
         ellipse = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
 #proteobacteria is in highest bundance in the sick groupings so lets check that out
#make a mixed effects model comparing prop of proteobact to animal health just for fun
#reference group is HH
geno_lm=lm(Proteobacteriap ~ site.animal.health, data= tl2)
summary(geno_lm)
#basic box plot of proteobact abundance grouped by site animal health 
gg<-ggplot(data=tl2, mapping=aes(x=animal.health, y=Proteobacteriap))+geom_boxplot(shape=1)
gg +
  labs(x='animal health', y='proteobacteria,p ')
#basic boxplot with outliers for Spirochaetes
gg<-ggplot(data=tl2, mapping=aes(x=animal.health, y=Spirochaetesp))+geom_boxplot(shape=1)
gg +
  labs(x='animal health', y='Spirochaetes, p')
#lets do a basic 1 way anova for proteobacteria to animal health
pro.way <-aov(Proteobacteriap~animal.health, data=tl2)
#look at the results
summary(pro.way)
#lets do a basic 1 way anova for spirochaetes to animal health 
spo.way <-aov(Spirochaetesp~animal.health, data=tl2)
#look at the results
summary(spo.way)
