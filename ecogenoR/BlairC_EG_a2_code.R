#assignment 2, comparing F3 gen to F1

#Blair


## Import or install the libraries that we're likely to need
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("DESeq2")
#install.packages("DESeq2")
install.packages("wesanderson")
install.packages("dplyr")
install.packages("ggpubr")
install.packages("vsn")
install.packages("eulerr")
install.packages("ggplot2")
install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot", force=TRUE)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(tidyverse)
library(ggplot2)
library(vsn)
library(ggbiplot)
#import cts matrix
countsT<- read.table("DE_counts_F3.txt", header=TRUE, row.names = 1)
#countsT<- DE_counts_F3
head(countsT)
dim(countsT)
# 24362    16
#rounding the decimal output from salmon so that the output data works with DESeq2
RcountsT<-round(countsT)  
head(RcountsT)
#import the sample descript table
conds<- read.table("RT_tonsa_F3_samples.txt", header=TRUE, row.names = 1)
#conds<- RT_tonsa_F3_samples
#conds<- ddply(conds, "index", transform,  )
head(conds)
#how many reads from each sample
mean(colSums(RcountsT)) #17727744
barplot(colSums(RcountsT), names.arg=colnames(RcountsT), cex.names=0.5, las=3) 
abline(h=mean(colSums(RcountsT)), col="blue", lwd=2) #shows the mean
mean(rowSums(RcountsT)) #11220.54
median(rowSums(RcountsT)) #2144
apply(RcountsT,2,mean) #rows
# AAAA_F3_REP1 AAAA_F3_REP2 AAAA_F3_REP3 AAAA_F3_REP4 
# 701.1270     723.2602     635.8753     734.2962 
# AAHH_F3_REP1 AAHH_F3_REP2 AAHH_F3_REP3 AAHH_F3_REP4 
# 754.6570     597.0389     733.5649     621.1106 
# HHAA_F3_REP1 HHAA_F3_REP2 HHAA_F3_REP3 HHAA_F3_REP4 
# 804.2911     718.5034     751.9464     717.2275 
# HHHH_F3_REP1 HHHH_F3_REP2 HHHH_F3_REP3 HHHH_F3_REP4 
# 620.4844     678.3669     734.9327     693.8526 
apply(RcountsT,1,mean) #columns
hist(apply(RcountsT,1,mean), xlim=c(0,2000), breaks=10000)
#negative binomial (counts) GLM
dds<-DESeqDataSetFromMatrix(countData=RcountsT, colData=conds,
                            design=~line+environment+line:environment)
dim(dds) #25279    16
#filtering with too few reads
#this overwrites the original object
dds<-dds[rowSums(counts(dds))>160]
dim(dds)
dds<-DESeq(dds)
resultsNames(dds)
# [1] "Intercept"                  "line_combined_vs_ambient"  
# [3] "environment_HH_vs_AA"       "linecombined.environmentHH"
#PCA to look at patterns of expression
#normalize the data
vsd<-vst(dds,blind=FALSE)
#save pca output as a dataset
data<-plotPCA(vsd,intgroup=c("line","environment"), returnData=TRUE)
percentVar<-round(100*attr(data,"percentVar"))
percentVar #47 11
ggplot(data, aes(PC1,PC2,color=environment,shape=line, groups=data[,2],ellipse = TRUE))+
  geom_point(size=4,alpha=0.85)+
  xlab(paste0('PC1:',percentVar[1],"% Variance"))+
  ylab(paste0('PC2:',percentVar[2],"% Variance"))
#order and summarize results from specific contrasts
resIxn<-results(dds,alpha=0.05)
resIxn<-resIxn[order(resIxn$padj),]
head(resIxn)
# log2 fold change (MLE): linecombined.environmentHH 
# Wald test p-value: linecombined.environmentHH 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN142181_c0_g4    384.448        3.11803  0.398309   7.82817
# TRINITY_DN143012_c0_g4    467.286       -3.31564  0.475434  -6.97391
# TRINITY_DN131723_c0_g1   1926.655        2.60636  0.387730   6.72209
# TRINITY_DN142181_c0_g18   364.882        3.11752  0.468239   6.65797
# TRINITY_DN145818_c5_g1    297.741        1.89854  0.295160   6.43225
# TRINITY_DN135177_c0_g1   5854.210        2.38301  0.396707   6.00699
# pvalue        padj
# <numeric>   <numeric>
#   TRINITY_DN142181_c0_g4  4.95009e-15 1.12976e-10
# TRINITY_DN143012_c0_g4  3.08244e-12 3.51753e-08
# TRINITY_DN131723_c0_g1  1.79140e-11 1.36283e-07
# TRINITY_DN142181_c0_g18 2.77627e-11 1.58407e-07
# TRINITY_DN145818_c5_g1  1.25730e-10 5.73907e-07
# TRINITY_DN135177_c0_g1  1.89001e-09 7.18927e-06
summary(resIxn)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 271, 1.1%
# LFC < 0 (down)     : 60, 0.24%
# outliers [1]       : 5, 0.02%
# low counts [2]     : 2451, 9.7%
# (mean count < 26)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
#######################
############################################## TEST FOR EFFECT OF ENVIRONMENT
####################### Likelhood ratio test

dds <- DESeqDataSetFromMatrix(countData = RcountsT, colData = conds, 
                              design = ~ line + environment)

dds <- DESeq(dds, test="LRT", reduced=~line)
# List the results you've generated
resultsNames(dds)
# [1] "Intercept"                "line_combined_vs_ambient"
# [3] "environment_HH_vs_AA"   
# Order and list and summarize results from specific contrasts
resEnv <- results(dds, alpha = 0.05)
resEnv <- resEnv[order(resEnv$padj),]
head(resEnv)
# log2 fold change (MLE): environment HH vs AA 
# LRT p-value: '~ line + environment' vs '~ line' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN138549_c1_g2    582.410       -1.94341  0.173407  118.0825
# TRINITY_DN138549_c2_g12   773.349       -2.01757  0.203760   91.4210
# TRINITY_DN150696_c2_g3    297.068        1.31754  0.163636   63.1253
# TRINITY_DN123676_c0_g2    179.431       -2.51746  0.309813   59.1190
# TRINITY_DN131329_c1_g1    213.660       -1.23500  0.158361   59.4117
# TRINITY_DN105043_c0_g1    101.714       -3.94548  0.471847   57.1227
# pvalue        padj
# <numeric>   <numeric>
#   TRINITY_DN138549_c1_g2  1.66325e-27 3.96651e-23
# TRINITY_DN138549_c2_g12 1.16138e-21 1.38483e-17
# TRINITY_DN150696_c2_g3  1.93963e-15 1.54188e-11
# TRINITY_DN123676_c0_g2  1.48423e-14 7.07917e-11
# TRINITY_DN131329_c1_g1  1.27903e-14 7.07917e-11
# TRINITY_DN105043_c0_g1  4.09446e-14 1.62741e-10
summary(resEnv)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 513, 2%
# LFC < 0 (down)     : 315, 1.2%
# outliers [1]       : 16, 0.063%
# low counts [2]     : 491, 1.9%
# (mean count < 19)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
resEnv <- resEnv[!is.na(resEnv$padj),]
degsEnv <- row.names(resEnv[resEnv$padj < 0.05,]) 
summary(degsEnv)
# Length     Class      Mode 
# 828 character character 
#######################
##############################################  TEST FOR EFFECT OF LINE
#######################

dds <- DESeqDataSetFromMatrix(countData = RcountsT, colData = conds, 
                              design = ~ environment + line)

dds <- DESeq(dds, test="LRT", reduced=~environment)
resultsNames(dds)
# [1] "Intercept"                "environment_HH_vs_AA"    
# [3] "line_combined_vs_ambient"
resLine <- results(dds, alpha = 0.05)
resLine <- resLine[order(resLine$padj),]
head(resLine)
# log2 fold change (MLE): line combined vs ambient 
# LRT p-value: '~ environment + line' vs '~ environment' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN132194_c0_g1  1229.949        1.53901  0.154825   94.6814
# TRINITY_DN121089_c0_g2   632.855        1.50181  0.163775   80.8863
# TRINITY_DN134798_c0_g1   152.989       -2.05171  0.224453   77.6336
# TRINITY_DN129890_c0_g4   153.602       -3.22080  0.341257   76.5916
# TRINITY_DN147342_c0_g4   151.496        1.86674  0.237816   58.6985
# TRINITY_DN134960_c1_g9  2869.532        1.94038  0.245758   57.2414
# pvalue        padj
# <numeric>   <numeric>
#   TRINITY_DN132194_c0_g1 2.23631e-22 5.43021e-18
# TRINITY_DN121089_c0_g2 2.39089e-19 2.90279e-15
# TRINITY_DN134798_c0_g1 1.24040e-18 1.00398e-14
# TRINITY_DN129890_c0_g4 2.10231e-18 1.27620e-14
# TRINITY_DN147342_c0_g4 1.83783e-14 8.92522e-11
# TRINITY_DN134960_c1_g9 3.85471e-14 1.56000e-10
summary(resLine)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 808, 3.2%
# LFC < 0 (down)     : 837, 3.3%
# outliers [1]       : 16, 0.063%
# low counts [2]     : 981, 3.9%
# (mean count < 21)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
resLine <- resLine[!is.na(resLine$padj),]
resLine
# out of 24282 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 808, 3.3%
# LFC < 0 (down)     : 837, 3.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 21)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
degsline <- row.names(resLine[resLine$padj < 0.05,])


#######################
##############################################  TEST FOR INTERACTION
#######################

dds <- DESeqDataSetFromMatrix(countData = RcountsT, colData = conds, 
                              design = ~ environment + line + environment:line)

dds <- DESeq(dds, test="LRT", reduced=~environment + line)
resultsNames(dds)

resInt <- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),]
head(resInt)
# log2 fold change (MLE): environmentHH.linecombined 
# LRT p-value: '~ environment + line + environment:line' vs '~ environment + line' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN142181_c0_g4    384.448        3.11803  0.398309   59.3373
# TRINITY_DN143012_c0_g4    467.286       -3.31564  0.475434   46.3684
# TRINITY_DN131723_c0_g1   1926.655        2.60636  0.387730   43.7887
# TRINITY_DN142181_c0_g18   364.882        3.11752  0.468239   42.7430
# TRINITY_DN145818_c5_g1    297.741        1.89854  0.295160   40.8415
# TRINITY_DN135177_c0_g1   5854.210        2.38301  0.396707   35.1139
# pvalue        padj
# <numeric>   <numeric>
#   TRINITY_DN142181_c0_g4  1.32838e-14 2.96667e-10
# TRINITY_DN143012_c0_g4  9.79853e-12 1.09415e-07
# TRINITY_DN131723_c0_g1  3.65817e-11 2.72326e-07
# TRINITY_DN142181_c0_g18 6.24250e-11 3.48535e-07
# TRINITY_DN145818_c5_g1  1.65090e-10 7.37393e-07
# TRINITY_DN135177_c0_g1  3.10977e-09 1.15751e-05
summary(resInt)
# out of 25279 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 235, 0.93%
# LFC < 0 (down)     : 48, 0.19%
# outliers [1]       : 5, 0.02%
# low counts [2]     : 2941, 12%
# (mean count < 28)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
resInt <- resInt[!is.na(resInt$padj),]
# out of 22333 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 235, 1.1%
# LFC < 0 (down)     : 48, 0.21%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 28)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
degsInt <- row.names(resInt[resInt$padj < 0.05,])
# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2802, 12%
# LFC < 0 (down)     : 1052, 4.3%
# outliers [1]       : 9, 0.037%
# low counts [2]     : 945, 3.9%
# (mean count < 20)
#check/visualize that our stats are working for us
### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validation that the normalization, model is working)
head(resInt)
# log2 fold change (MLE): environmentHH.linecombined 
# LRT p-value: '~ environment + line + environment:line' vs '~ environment + line' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN142181_c0_g4    384.448        3.11803  0.398309   59.3373
# TRINITY_DN143012_c0_g4    467.286       -3.31564  0.475434   46.3684
# TRINITY_DN131723_c0_g1   1926.655        2.60636  0.387730   43.7887
# TRINITY_DN142181_c0_g18   364.882        3.11752  0.468239   42.7430
# TRINITY_DN145818_c5_g1    297.741        1.89854  0.295160   40.8415
# TRINITY_DN135177_c0_g1   5854.210        2.38301  0.396707   35.1139
# pvalue        padj
# <numeric>   <numeric>
#   TRINITY_DN142181_c0_g4  1.32838e-14 2.96667e-10
# TRINITY_DN143012_c0_g4  9.79853e-12 1.09415e-07
# TRINITY_DN131723_c0_g1  3.65817e-11 2.72326e-07
# TRINITY_DN142181_c0_g18 6.24250e-11 3.48535e-07
# TRINITY_DN145818_c5_g1  1.65090e-10 7.37393e-07
# TRINITY_DN135177_c0_g1  3.10977e-09 1.15751e-05
d <-plotCounts(dds, gene="TRINITY_DN131561_c0_g1", intgroup = (c("line","environment")), returnData=TRUE)
d
# count     line environment
# AAAA_F3_REP1  6175.492  ambient          AA
# AAAA_F3_REP2  3123.259  ambient          AA
# AAAA_F3_REP3  5636.901  ambient          AA
# AAAA_F3_REP4  6431.794  ambient          AA
# AAHH_F3_REP1  6717.726  ambient          HH
# AAHH_F3_REP2  5270.274  ambient          HH
# AAHH_F3_REP3 11278.976  ambient          HH
# AAHH_F3_REP4   420.528  ambient          HH
# HHAA_F3_REP1  4467.336 combined          AA
# HHAA_F3_REP2  2448.682 combined          AA
# HHAA_F3_REP3  2121.085 combined          AA
# HHAA_F3_REP4  2130.135 combined          AA
# HHHH_F3_REP1 10138.211 combined          HH
# HHHH_F3_REP2  3792.778 combined          HH
# HHHH_F3_REP3 12570.410 combined          HH
# HHHH_F3_REP4  5906.772 combined          HH
#aesthetic mapping
p <-ggplot(d, aes(x=environment, y=count, color=line, shape=line, group=line)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p

#######################
############################################## PLOT OVERLAPPING DEGS IN VENN DIAGRAM
#######################

library(eulerr)

# Total
lde<-length(degsEnv)  # 828
ldl<-length(degsline)  # 1645
ldi<-length(degsInt)  # 283


# Intersections
int_el<-length(intersect(degsEnv,degsline))  # 141
int_ei<-length(intersect(degsEnv,degsInt))  # 14
int_il<-length(intersect(degsInt,degsline))  # 32

intEL <- intersect(degsEnv,degsline)
mid<-length(intersect(degsInt,intEL)) # 7

# Number unique
#env<-lde-int_ei-int_el-mid
#line<-ldl-int_el-int_ei-mid
#int<-ldi-int_ei-int_el-mid
env<-828-14-141-7 # 666
line<-1645-141-32-7 # 1465
int<-283-14-32-7 # 230
#venn diagram for just f3
fit1 <- euler(c("Env" = 666, "Line" = 1465, "Interaction" = 230, "Env&Line" = 141, "Env&Interaction" = 14, "Line&Interaction" = 32, "Env&Line&Interaction" = 7))

plot(fit1,  lty = 1:3, quantities = TRUE)

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))
#venn diagram for env f1 v. f3
fit2 <- euler(c("EnvF1" = 394, "EnvF1&EnvF3" = 54, "EnvF3" = 230))

plot(fit2,  lty = 3:1, quantities = TRUE)

plot(fit2, quantities = TRUE, fill = "transparent",
     lty = 3:1,
     labels = list(font = 4))
#venn diagram for line f1 v f3
fit3 <- euler(c("LineF1" = 177, "LineF1&LineF3" = 49, "LineF3" = 1540))

plot(fit3,  lty = 1:3, quantities = TRUE)

plot(fit3, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))
#venn diagram for interaction f1 v. f3
fit4 <- euler(c("IntF1" = 3749, "IntF1&IntF3" = 105, "IntF3" = 178))

plot(fit4,  lty = 1:3, quantities = TRUE)

plot(fit4, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))
# Heatmap of top 20 genes sorted by pvalue

library(pheatmap)

# By int

topgenes <- head(rownames(resInt),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)

# By line

topgenes <- head(rownames(resLine),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)
# By environment

topgenes <- head(rownames(resEnv),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("line","environment")])
pheatmap(mat, annotation_col=df)