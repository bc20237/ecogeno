#hwk 3 
library(ggplot2)
library(gridExtra)
library(ggbiplot)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(gridExtra) 

setwd("/Users/blairchristensen/Desktop/ecological genomics/popgenomics")
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

grid.arrange(plotmed_DD0, plotmean, plotgdd, nrow = 3)
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
# Define association outliers as the upper 0.1% of p-values
#define snps
snps<-read.table("Chr19.kept.sites",sep="\t", header=T)
#########  rust  #########
ffreeze2 <- cbind(snps, ffreeze2[,-c(1:2)])
ffreeze2$outlier = ifelse(ffreeze2$P<quantile(rust2$P,0.001),2,1) #upper .1% 

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
gdd2$outlier = ifelse(gdd2$P<quantile(gdd2$P,0.001),2,1)

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
dd02$outlier = ifelse(dd02$P<quantile(dd02$P,0.001),2,1)

p3 <- ggplot(dd02,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=dd02$outlier, color=dd02$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("dd0")

p3

grid.arrange(p1, p2, p3, nrow = 3)
# Get outliers for a given trait association:

ffreeze_outliers <- ffreeze2[which(ffreeze2$outlier==2),c(2,3,9)]
gdd_outliers <- gdd2[which(gdd2$outlier==2),c(2,3,9)]
dd0_outliers <- dd02[which(dd02$outlier==2),c(2,3,9)]
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
p1 <- ggplot(AF3[,-3],aes(x=window,y=AvgAncFreq)) +
  geom_line(size=0.8, color="blue") + 
  xlab("Position (bp) along chromosome") +
  ylab("Frequency P. trichocarpa ancestry") +
  geom_hline(yintercept=mean(AF2$ALT_FREQS), color = "red") + 
  geom_hline(yintercept=upper, linetype="dashed", color = "red") + 
  geom_hline(yintercept=lower, linetype="dashed", color = "red") +
  ggtitle("Chr19: Local ancestry")

p1


grid.arrange(p1, p2, p3, p4, nrow = 4)
# Get the betas from each trait and look at pleiotropy between traits
betas <- cbind(ffreeze2[,c(1:3,9)],gdd2[,9],dd02[,9],budflush2[,9])
names(betas) = c("CHROM","POS","ID","beta_ffreeze","beta_gdd","beta_dd0", "beta_flush")
str(betas)

cor(betas[,4:6],betas[4:6])

plot(beta$beta_ffreeze,betas$beta_flush)

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

# Define the genomic ranges of the Fst bins
fstGR <- GRanges(CHR,IRanges(fst[,"BIN_START"],fst[,"BIN_END"]))

# Define what should be a "low" value of Fst.  Here, we'll try the lowest 10% of windows that don't incude the centromere.  Can play with this if you want (choosing the last value in the quantile function call). 
#if we make it too small, there are hardly any instances to look at
#if too big, then wash the signal out
fstThreshold = quantile(fst$MEAN_FST,0.1)

# Call outliers (=1) or non-outliers (=2) based on fst Threshold
fstGR$outlier <- ifelse(fst$MEAN_FST<fstThreshold,1,2)

# Grab just the outlier Fst windows
fstCand <- subset(fstGR, outlier==1)

raisd <- read.table(paste0("RAiSD_Report.",CHR,".",CHR), sep="\t",header=F) # Import the RAiSD results

# Define RAiSD genomic ranges based on first and last SNP within the RAiSD windows of 50 SNPs each.   How does this relate to the Fst window size?
raisdGR <- GRanges(CHR,IRanges(raisd[,2],raisd[,3]))

# Define RAiSD outliers, based on the upper 1% of SNPs
raisdThreshold = quantile(raisd$V7,0.99)

raisdGR$outlier <- ifelse(raisd$V7>raisdThreshold,"outlier","non")
#Since the RAiSD results are 1 per SNP, and adjacent SNPs have highly similar values (because of linkage), we want to merge runs of adjacent outlier positions together into outlier windows.
raisdGR_out <- unlist(reduce(split(raisdGR, ~outlier)))
raisdGR_out$outlier <- names(raisdGR_out)
raisdCand <- subset(raisdGR_out, outlier=="outlier")
# Use GenomicRanges to pull out the RAiSD outlier windows that overlap the lowest Fst windows
overlap <- subsetByOverlaps(raisdCand, fstCand)

length(overlap) # Number of the RAiSD sweep candidate loci that overlap with the lowest n% of Fst windows
# How many windows overlap between RAiSD and Fst?
#   10 for chromosome 19
#   How much is a lot? How few is a little?
#   
#   Step 4: Test for significance of observed overlapping bins by randomizing the Fst values and recalcing the overlap
# A little code I wrote to permute Fst values randomly among the 50kb windows, and re-calculate the overlap with the RAiSD outliers. 
# Can set number of permutation replicates (NumPerm)
# 1000 permutations seems to give a decent randomized distribution

NumPerm=1000

Perm = NA
for(i in 1:NumPerm){
  FstSamp = cbind(fst[,c(2,3)], fst[sample(nrow(fst)),6])
  names(FstSamp) = c("Start", "Stop", "FST")
  FstRand = FstSamp[which(FstSamp[3]<fstThreshold),]  
  FstRand_ranges <- GRanges(seqnames=CHR, ranges=IRanges(FstRand[,1],FstRand[,2]))
  FstRand_ranges_red <- reduce(FstRand_ranges)
  Perm[i] = sum(countOverlaps(raisdCand, FstRand_ranges_red))
}

# Plot the random distribution with a blue line marking the observed overlap
hist(Perm, col="gray", main="Randomization test", xlab="Overlapping windows of Fst and RAiSD outliers for Chr19")
abline(v=length(overlap), col="blue", lwd=3)

# Calculate the p-pvalue based on the rank of the observed value in the randomized distirbution.
# Note the 1-tailed test here (alpha = 0.05)

p_value = 1-ecdf(Perm)(length(overlap))
p_value # 0.063

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
candGenes <- subsetByOverlaps(genes, overlap)

write.table(candGenes$geneID, paste0("candGenes",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")



