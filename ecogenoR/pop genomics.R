#read in datasets manually
library(maps)  # If you don't have this library, you'll need to run 
install.packages("maps")
library(plotrix) # If you don't have this library, you'll need 
install.packages("plotrix")
install.packages("gridExtra")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicFeatures")

library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(gridExtra) # If needed, install.packages("gridExtra")


setwd("/Users/blairchristensen/Desktop/ecological genomics/popgenomics") #set the path to where your downloaded files are on your laptop

meta <- read.table("Combined_Transect_Sampling_Data_2020.txt", sep="\t",header=T)  
#meta<-metadata
Qscores <- read.table("poplar_hybrids.LDpruned.5.Q", sep=" ",header=F)
#Qscores<-poplar
names(Qscores) <- c("K1","K2","K3","K4","K5")  # Customize to your level of K!

tiff("Admix.tiff", width = 10, height = 7, units = "in",
     res = 300)
map("world",xlim=c(-160,-100),ylim=c(40,70),fill=T, col="lightgrey")

title("Admixture K=5") # Change for your value of K
map.axes(cex.axis=1.5)
for(i in 1:nrow(meta)){
  floating.pie(meta$Longitude[i], meta$Latitude[i], 
               c(Qscores$K1[i],Qscores$K2[i],Qscores$K3[i], Qscores$K4[i],Qscores$K5[i]),
               col=c("yellow","green","blue","red", "purple","orange"),
               radius=0.5)
}

# Customize the above to your level of K wherever you see the ellipses (...)!

dev.off()

# Don't forget to save your R script!
library(ggplot2)

Qscores_new <- Qscores # Rename to the new file!
names(Qscores_new) <- c("K1","K2","K3", "K4", "K5")  # Customize to the level of K!

# Calculate Shannon Diversity across the K different Q-scores per individual

K=5  # Change X to the level of K we're investigating

tmp=numeric()

for(i in 1:nrow(Qscores_new)){
  for(j in 1:K){
    tmp[j] = Qscores_new[i,j]*log(Qscores_new[i,j])
  }
  Qscores_new$ShDiv[i] = -1*sum(tmp)
}
#bring in heterozygosity data
het <- read.table("het_maf05.het", sep="\t", header=T) 

str(het) # What's in this dataframe? 
# data dict
#INDV = Individual sample ID
# O.HOM. = Number of sites genome-wide that were observed to be homozygous (either 0|0 or 1|1)
# E.HOM. = Number of sites genome-wide predicted to be homozygous based on Hardy-Weinberg expectations (=1-2pq across all loci)
# N_SITES = Number of SNPs included in the calculation
# F = Inbreeding coefficient (=1-obsHet/expHet)
# Combine the meta data, heterozygosity, and admixture data

het2 <- cbind(meta,het,Qscores_new) #bind the het results with the meta data

# How does F vary within each transect?
#bimodal, non gaussian, skewed data
ggplot(het2, aes(x=Transect, y=F, color=Transect)) +
  geom_dotplot(binaxis='y', binwidth=0.01, stackdir='center')
# Plot admixture diversity spatially
ggplot(het2, aes(x=Longitude, y=Latitude, color=ShDiv)) +
  geom_point(size=4, shape=20)

# And finally, are more admixed individuals more heterozygous in their genomes?

ggplot(het2, aes(x=F, y=ShDiv, color=Transect)) +
  geom_point(size=4, shape=20)
#higher heterozygosity == lower values of F
#more admixed at k=5, tgrending towards complete random mating between the species
#chilcotin and jasper transcest especially
#lower and lower values of ancestry diversity --> start to see more homozygosity/less admixture in the genome
#look at phenotypic/performance differences and the effect of environment
cor.test(het2$F,het2$ShDiv)

pi <- read.table("Chr19.windowed.pi",sep="\t", header=T)
str(pi)

fst <- read.table("Bals_Tricho_All.windowed.weir.fst", sep="\t",header=T)
str(fst)

cent <- read.table("Chr19_centromere.txt", sep="\t",header=F)
centromere = mean(c(cent$V2,cent$V3))

raisd <- read.table("RAiSD_Report.Chr19.Chr19", sep="\t",header=F)
str(raisd)

p1 <- ggplot(pi,aes(x=BIN_START,y=PI/mean(PI))) +
  geom_line(size=0.25, color="blue") + 
  geom_point(aes(x=centromere,y=1, size=100), show.legend=F) +
  xlim(400000,6000000) +
  ggtitle("Chomosome 19: Nucleotide diversity and Fst in 50kb sliding windows") +
  xlab("") +
  ylab("Scaled nucleotide diversity")

p2 <- ggplot(fst,aes(x=BIN_START,y=MEAN_FST/mean(MEAN_FST))) +
  geom_line(size=0.25, color="red") +
  geom_point(aes(x=centromere,y=1, size=100), show.legend=F) +
  xlim(10000000,12000000) + 
  ylab("Scaled Fst")

p3 <- ggplot(raisd,aes(x=V1,y=V7/mean(V7))) +
  geom_point(size=0.25, color="black") +
  xlim(0,max(raisd$V3)) + 
  xlab("Position along chromosome (bp)") +
  ylab("RAiSD u-stat")

grid.arrange(p1, p2, nrow = 2)
#Nov 1 2021
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

# Define RAiSD genomic ranges based on first and last SNP within the RAiSD windows of 50 SNPs each.   How deos this relate to the Fst window size?
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
