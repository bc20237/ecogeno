#hwk 3 
library(ggplot2)
library(gridExtra)
library(ggbiplot)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(gridExtra) 
library(ggplot2)

setwd("/Users/blairchristensen/Desktop/ecological genomics/popgenomics")

# Read in list of positions
snps2 <- read.table("Chr19.kept.sites",sep="\t", header=T)

# Read in the local ancestry frequencies from Plink
AF2 <- read.table("Chr19_LAI_freq.afreq", skip=1,sep="\t",header=F)  

# Note the skip=1 here.  
# This skips the first line, since Plink's header line doesn't play well with R.  
# We'll define our own header line below.

names(AF2) = c("CHROM",  "ID",   "REF",  "ALT",  "ALT_FREQS",    "OBS_CT")

AF3 <- cbind(snps2,AF2)

str(AF3) # How does it look?

# A simple plot:

pa <- ggplot(AF3[,-3],aes(x=POS,y=ALT_FREQS)) +
  geom_line(size=0.25, color="blue") + 
  xlab("Position (bp) along chromosome") +
  ylab("Frequency P. trichocarpa ancestry")

pa
# Get the list of admixed individuals:
Admixed <- read.table("Admixed.Inds",header=F)

# Get the meta data:
meta <- read.table("Combined_Transect_Sampling_Data_2020.txt", sep="\t",header=T)

# merge them together:
meta_admx <- merge(meta, Admixed, by.x="ID", by.y="V1")
str(meta_admx)  

# Read in the Admixture coefficients for KBals that we made from the K=5 file:
KBals <- read.table("Admixed_KBals", sep="\t", header=F)
names(KBals) = c("ID","KBals")

# Second merge:
meta_admx_KBals <- merge(meta_admx,KBals,by="ID")


# Bring in phenotype data:
pheno <- read.table("VT_Garden_Phenotypes_2021.txt",sep="\t",header=T)
clim <- read.table("climatedata.txt",sep="\t",header=T)
# Merge pheno data with meta and KBals:
meta_admx_KBals_clim <- merge(meta_admx_KBals,clim,by="ID")
meta_admx_KBals_cp <- merge(meta_admx_KBals_clim,pheno,by="ID")
# This is the average date of the last freezing event in spring, after which temperatures stay above 0C for the rest of the growing season.
plotmean <- ggplot(meta_admx_KBals_cp,aes(x=KBals,y=mean_finalFreeze, color=Transect.x, ellipse = TRUE)) +
  geom_point(size=2) +
  xlab("Proportion P. balsamifera ancestry") +
  ylab("final freeze") 

plotmean

# This is the average number of growing degree days (a measure of spring warming) that accumulate from Jan01 to the date of last freeze in spring.
plotgdd <- ggplot(meta_admx_KBals_cp,aes(x=KBals,y=mean_cGDDfreeze, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("mean_cGDDfreeze") 

plotgdd

#This is the mean number of chilling degree days across the year.  Chilling degree days are thought to be a cue that plants use to determine how much winter they’ve passed through, which can be used to decide when to break dormancy (flush).
plotmed_DD0 <- ggplot(meta_admx_KBals_cp,aes(x=KBals,y=med_DD0, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("med_DD0") 
#individuals with trichocarp ancestry take longer in spring to break bud, higher balsamifera== earlier
#populations that have evolved under warmer climates require more heat to bud
#natural pops in warm know that they can get a few warm days even in winter
plotmed_DD0
# Bud flush
plotBudflush <- ggplot(meta_admx_KBals_cp,aes(x=KBals,y=FLUSH, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Bud flush") 
#individuals with trichocarp ancestry take longer in spring to break bud, higher balsamifera== earlier
#populations that have evolved under warmer climates require more heat to bud
#natural pops in warm know that they can get a few warm days even in winter
plotBudflush

grid.arrange(plotmed_DD0, plotmean, plotgdd, plotBudflush, nrow = 4)
#bud set for trait
budflush <- read.table("plink2.FLUSH.glm.linear",skip=1,sep="\t",header=F)
names(budflush) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
budflush <- budflush[which(budflush$TEST=="ADD"),]
budflush2 <- cbind(snps, budflush[,-c(1,2)])
budflush2$outlier = ifelse(budflush2$P<quantile(budflush2$P,0.001),2,1)
#regressions for each of these plots
# linear models testing trait ~ genome-wide admixture association
#mean final freeze
summary(lm(mean_finalFreeze~KBals + Transect.x, data=meta_admx_KBals_cp))

summary(lm(mean_cGDDfreeze~KBals + Transect.x, data=meta_admx_KBals_cp))

summary(lm(med_DD0~KBals + Transect.x, data=meta_admx_KBals_cp))

lm(formula = med_DD0 ~ KBals + Transect.x, data = meta_admx_KBals_cp)

# What about the effects of local ancestry within the genome, after controlling for genome-wide ancestry effects as a covariate in the GLM model?
######  Bring in Association results from Plink   ######

ffreeze <- read.table("plink2.mean_finalFreeze.glm.linear",skip=1,sep="\t",header=F)
names(ffreeze) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
ffreeze2 <- ffreeze[which(ffreeze$TEST=="ADD"),]
head(ffreeze2)
#as you insert balsif genome wide it tends to decrease  
# Define association outliers as the upper 1% of p-values
#define snps
snps<-read.table("Chr19.kept.sites",sep="\t", header=T)
#########
ffreeze2 <- cbind(snps, ffreeze2[,-c(1:2)])
ffreeze2$outlier = ifelse(ffreeze2$P<quantile(ffreeze2$P,0.01),2,1) #upper .1% 

# Get outliers for a given trait association:

ffreeze_outliers <- ffreeze2[which(ffreeze2$outlier==2),c(2,3,9)]
gdd_outliers <- gdd2[which(gdd2$outlier==2),c(2,3,9)]
dd0_outliers <- dd02[which(dd02$outlier==2),c(2,3,9)]
flush_outliers <- ffreeze2[which(ffreeze2$outlier==2),c(2,3,9)]
#plot ancestry
# Read in list of positions
snps <- read.table("Chr19.kept.sites",sep="\t", header=T)

# Plot freq of LAI along chr
AF <- read.table("Chr19_LAI_freq.afreq", skip=1,sep="\t",header=F)
names(AF) = c("CHROM",  "ID",   "REF",  "ALT",  "ALT_FREQS",    "OBS_CT")
str(AF)

AF2 <- cbind(snps,AF)

windows <- seq(1,max(AF2$POS),5e4)
AF_windows <- numeric()

for(i in 1:length(windows)){
  tmp=AF2[which(AF2$POS>windows[i] & AF2$POS<windows[i+1]),"ALT_FREQS"]
  ancfreq=mean(tmp)
  AF_windows[i] = ancfreq
}

AF3 <- as.data.frame(cbind(windows,AF_windows))
names(AF3) = c("window","AvgAncFreq")
#
upper = mean(AF3$AvgAncFreq,na.rm=T) + 2*sd(AF3$AvgAncFreq,na.rm=T)
lower = mean(AF3$AvgAncFreq,na.rm=T) - 2*sd(AF3$AvgAncFreq,na.rm=T)

outliers_upper = AF3[which(AF3$AvgAncFreq>upper),]
outliers_lower = AF3[which(AF3$AvgAncFreq<lower),]

# Print the outlier regions out
outliers_upper
outliers_lower

# And finally, make the 4-panel plot with the trait associations
pta <- ggplot(AF3[,-3],aes(x=window,y=AvgAncFreq)) +
  geom_line(size=0.8, color="blue") + 
  xlab("Position (bp) along chromosome") +
  ylab("Frequency P. trichocarpa ancestry") +
  geom_hline(yintercept=mean(AF2$ALT_FREQS), color = "red") + 
  geom_hline(yintercept=upper, linetype="dashed", color = "red") + 
  geom_hline(yintercept=lower, linetype="dashed", color = "red") +
  ggtitle("Chr19: Local ancestry")

pta



# Get the betas from each trait and look at pleiotropy between traits
betas <- cbind(ffreeze2[,c(1:3,9)],gdd2[,9],dd02[,9],budflush2[,9])
names(betas) = c("CHROM","POS","ID","beta_ffreeze","beta_gdd","beta_dd0", "beta_flush")
str(betas)

cor(betas[,4:6],betas[4:6])

plot(beta$beta_ffreeze,betas$beta_flush)


p1 <- ggplot(ffreeze2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=ffreeze2$outlier, color=ffreeze2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("final freeze")

p1

####### GDD #########
gdd <- read.table("plink2.mean_cGDDfreeze.glm.linear",skip=1,sep="\t",header=F)
names(gdd) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
gdd <- gdd[which(gdd$TEST=="ADD"),]
gdd2 <- cbind(snps, gdd[,-c(1,2)])
gdd2$outlier = ifelse(gdd2$P<quantile(gdd2$P,0.01),2,1)

p2 <- ggplot(gdd2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=gdd2$outlier, color=gdd2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("gdd")

p2

#########  med DD0  #########
dd0 <- read.table("plink2.med_DD0.glm.linear",skip=1,sep="\t",header=F)
names(dd0) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
dd0 <- dd0[which(dd0$TEST=="ADD"),]
dd02 <- cbind(snps, dd0[,-c(1,2)])
dd02$outlier = ifelse(dd02$P<quantile(dd02$P,0.01),2,1)

p3 <- ggplot(dd02,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=dd02$outlier, color=dd02$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("dd0")

p3
budflush <- read.table("plink2.FLUSH.glm.linear",skip=1,sep="\t",header=F)
names(budflush) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
budflush <- budflush[which(budflush$TEST=="ADD"),]
budflush2 <- cbind(snps, budflush[,-c(1,2)])
budflush2$outlier = ifelse(budflush2$P<quantile(budflush2$P,0.01),2,1)
p4 <- ggplot(budflush2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=budflush2$outlier, color=budflush2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("bud flush")

p4
#grid.arrange(p1, p2, p3, nrow = 3)

p5 <- ggplot(betas,aes(x=beta_ffreeze,y=beta_flush)) +
  geom_point(color="darkgray") + 
  xlab("first freeze") +
  ylab("bud flush") +
  ggtitle("Correlation of first freeze and bud flush")
p5

p6 <- ggplot(betas,aes(x=beta_gdd,y=beta_flush)) +
  geom_point(color="darkgray") + 
  xlab("gdd") +
  ylab("bud flush") +
  ggtitle("Correlation of gdd and bud flush")

p6

p7 <- ggplot(betas,aes(x=beta_dd0,y=beta_flush)) +
  geom_point(color="darkgray") + 
  xlab("dd0") +
  ylab("bud flush") +
  ggtitle("Correlation of dd0 and bud flush")

p7

#genomic ranges
#How can we get a test for association between RAiSD sweep regions and regions of low Fst?
# Bring in the Fst outputs we generated last week and identify windows with low Fst
# Bring in the RAiSD outputs and identify very high values of the u-stat (candidates for selection)
# Use GenomicRanges() to find their overlaps
# Test for significance by randomizing the values of meanFst among the windows many times, and estimate how much overlap there is in the randomized distribution
####### Randomization test for Fst and sweep outliers #########

fst <- read.table("Bals_Tricho_All.windowed.weir.fst", sep="\t",header=T) # Import the Fst results

cent <- read.table("Chr19_centromere.txt", sep="\t",header=F) # Import the centromere coordinates

fst <- fst[-which(fst$BIN_START>cent$V2 & fst$BIN_END<cent$V3),] # Mask Fst windows in the centromere region


# Calculate Genomic Ranges for outlier bins

CHR="Chr19"  # Customize to your chromosome.  Make sure syntax is exact!
#define ya windows:
#POS : which is the SNP being tested, also the center of the window
#P which is the pvalue for the association test
#outlier : outlier status
budflushGR <- GRanges(CHR,IRanges(budflush2$POS-2.5e4,budflush2$POS+2.5e4),POS=budflush2$POS, P=budflush2$P, outlier=budflush2$outlier)
#final freeze
ffreezeGR <- GRanges(CHR,IRanges(ffreeze2$POS-2.5e4,ffreeze2$POS+2.5e4),POS=ffreeze2$POS, P=ffreeze2$P, outlier=ffreeze2$outlier)
#dd0
dd0GR <- GRanges(CHR,IRanges(dd02$POS-2.5e4,dd02$POS+2.5e4),POS=dd02$POS, P=dd02$P, outlier=dd02$outlier)
#gdd
gddGR <- GRanges(CHR,IRanges(gdd2$POS-2.5e4,gdd2$POS+2.5e4),POS=gdd2$POS, P=gdd2$P, outlier=gdd2$outlier)

#reduced windows to capture outlier regions
#budflush
budflushGRout <- unlist(reduce(split(budflushGR, ~outlier)))
budflushGRout$outlier <- names(budflushGRout)
budflushGRCand <- subset(budflushGRout, outlier==2)

budflushGRCand # Print the candidate regions
#ffreeze
ffreezeGRout <- unlist(reduce(split(ffreezeGR, ~outlier))) 
ffreezeGRout$outlier <- names(ffreezeGRout)
ffreezeGRCand <- subset(ffreezeGRout, outlier==2)

ffreezeGRCand # Print the candidate regions

#gdd
gddGRout <- unlist(reduce(split(gddGR, ~outlier)))
gddGRout$outlier <- names(gddGRout)
gddGRCand <- subset(gddGRout, outlier==2)

gddGRCand # Print the candidate regions
#dd0
#gdd
dd0GRout <- unlist(reduce(split(dd0GR, ~outlier)))
dd0GRout$outlier <- names(dd0GRout)
dd0GRCand <- subset(dd0GRout, outlier==2)

dd0GRCand # Print the candidate regions

overlap_BF_ffreeze <- subsetByOverlaps(budflushGRCand, ffreezeGRCand)
length(overlap_BF_ffreeze)
#1
overlap_BF_ffreeze # Print the overlapping regions
#overlap flush gdd
overlap_BF_gdd <- subsetByOverlaps(budflushGRCand, gddGRCand)
length(overlap_BF_gdd)
#1
overlap_BF_gdd
#dd0
overlap_BF_dd0 <- subsetByOverlaps(budflushGRCand, dd0GRCand)
length(overlap_BF_dd0)
#1
overlap_BF_dd0
#experiments
overlaps_cli <- subsetByOverlaps(dd0GRCand, gddGRCand)
length(overlaps_cli)
head(overlaps_cli)
# Import the GFF annotation file and make a transcript database
txdb <- makeTxDbFromGFF("Ptrichocarpa_533_v4.1.gene.gff3.gz", format="gff3")

txdb

# How many chromosomes are present?
head(seqlevels(txdb))
#19
# Subset the database for just your chromosome of interest
seqlevels(txdb) <- CHR # subset for just your chromosome

# Reduce the transcript database to just the non-redundant gene names, instead of multiple entries for all the variant transcript types per gene
genes <- unlist(reduce(transcriptsBy(txdb, by="gene"))) 
genes$geneID <- names(genes)
#Now we’ll use GenomicRanges() just like before to find the genes that overlap in the intervals of our candidate regions of overlapping lowFst/highRAiSD
candGenesffreeze <- subsetByOverlaps(genes, overlap_BF_ffreeze)
candGenesgdd <- subsetByOverlaps(genes, overlap_BF_gdd)
candGenesdd0 <- subsetByOverlaps(genes, overlap_BF_dd0)


write.table(candGenesffreeze$geneID, paste0("candGenesffreeze",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")

write.table(candGenesgdd$geneID, paste0("candGenesgdd",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")
write.table(candGenesdd0$geneID, paste0("candGenesdd0",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")

grid.arrange(p1, p2, p3, p4,  pta, nrow = 5) #p2 gdd outliers,p3 dd0 outliers, pta local ancestry

candGenesffreeze <- subsetByOverlaps(genes, ffreezeGRCand)
candGenesgdd <- subsetByOverlaps(genes, gddGRCand)
candGenesdd0 <- subsetByOverlaps(genes, dd0GRCand)
candGenesBF <- subsetByOverlaps(genes, budflushGRCand)
write.table(candGenesffreeze$geneID, paste0("candGenesffreeze",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")
write.table(candGenesgdd$geneID, paste0("candGenesgdd",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")
write.table(candGenesdd0$geneID, paste0("candGenesdd0",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")
write.table(candGenesBF$geneID, paste0("candGenesBF",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")
# linear models testing trait ~ genome-wide admixture association
summary(lm(mean_finalFreeze~KBals + Transect.x, data=meta_admx_KBals_cp))
# Residual standard error: 8.861 on 388 degrees of freedom
# Multiple R-squared:  0.2849,	Adjusted R-squared:  0.272 
# F-statistic: 22.08 on 7 and 388 DF,  p-value: < 2.2e-16
summary(lm(med_DD0~ KBals+ Transect.x, data=meta_admx_KBals_cp))
# Residual standard error: 17.01 on 388 degrees of freedom
# Multiple R-squared:  0.6679,	Adjusted R-squared:  0.6619 
# F-statistic: 111.5 on 7 and 388 DF,  p-value: < 2.2e-16
summary(lm(mean_cGDDfreeze~KBals + Transect.x, data=meta_admx_KBals_cp))
# Residual standard error: 91.3 on 388 degrees of freedom
# Multiple R-squared:  0.5671,	Adjusted R-squared:  0.5593 
# F-statistic:  72.6 on 7 and 388 DF,  p-value: < 2.2e-16
#relationship between mean gdd freeze and bud flush
summary(lm(mean_cGDDfreeze~FLUSH + KBals+ Transect.x, data=meta_admx_KBals_cp))
#relationship between med dd0 and bud flush
summary(lm(med_DD0~FLUSH + KBals+ Transect.x, data=meta_admx_KBals_cp))
#relationship between final freeze and bud flush controlling for kbals and transect
summary(lm(mean_finalFreeze~FLUSH + KBals+ Transect.x, data=meta_admx_KBals_cp))

grid.arrange(p5, p6, p7, nrow = 3)

