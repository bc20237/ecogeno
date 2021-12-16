#Practice in clss RT tonsa data 10/6/2021
#Blair


## Import or install the libraries that we're likely to need
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

BiocManager::install("vsn")
#install.packages("DESeq2")
install.packages("wesanderson")
install.packages("dplyr")
install.packages("ggpubr")
install.packages("vsn")
install.packages("eulerr")
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)
library(eulerr)
#import cts matrix
countsT<- read.table("DE_counts_F1.txt", header=TRUE, row.names = 1)
head(countsT)

# AAAA_F1_REP1 AAAA_F1_REP2 AAAA_F1_REP3
# TRINITY_DN100115_c0_g1    310.00378    373.26532    506.80746
# TRINITY_DN100154_c0_g1   1093.44804    674.06929    628.66253
# TRINITY_DN100163_c0_g1    256.29813    176.65295    105.83147
# TRINITY_DN100287_c0_g1   1153.42268    834.78395    834.71318
# TRINITY_DN100291_c0_g1     19.55793     22.04806     44.86396
# TRINITY_DN100299_c0_g1    891.43208    455.89363    530.38247
# AAAA_F1_REP4 AAHH_F1_REP1 AAHH_F1_REP2
# TRINITY_DN100115_c0_g1    518.52031     112.9249    214.20766
# TRINITY_DN100154_c0_g1    771.10416    1056.0632    802.33107
# TRINITY_DN100163_c0_g1    170.63344     219.6631    139.31818
# TRINITY_DN100287_c0_g1    984.53094     961.8786    993.47256
# TRINITY_DN100291_c0_g1     39.26675      14.1519     17.97648
# TRINITY_DN100299_c0_g1    455.29677     713.5501    688.11495
# AAHH_F1_REP3 AAHH_F1_REP4 HHAA_F1_REP1
# TRINITY_DN100115_c0_g1    191.94707     85.35907     74.00049
# TRINITY_DN100154_c0_g1    854.45206    997.05517    846.02406
# TRINITY_DN100163_c0_g1    175.25611    198.64220    184.13711
# TRINITY_DN100287_c0_g1    846.15446    978.32033    847.49286
# TRINITY_DN100291_c0_g1     27.03752     16.77744     38.13279
# TRINITY_DN100299_c0_g1    824.05111    710.25195    956.11714
# HHAA_F1_REP2 HHAA_F1_REP3 HHAA_F1_REP4
# TRINITY_DN100115_c0_g1     37.63530    121.96639     60.56209
# TRINITY_DN100154_c0_g1   1103.74913   1039.51204    858.10746
# TRINITY_DN100163_c0_g1    210.26998    269.19308    193.66068
# TRINITY_DN100287_c0_g1   1003.01746   1041.52537   1114.88586
# TRINITY_DN100291_c0_g1     17.52053     20.06726     18.92643
# TRINITY_DN100299_c0_g1    812.00280    690.43073    728.69995
# HHHH_F1_REP1 HHHH_F1_REP2 HHHH_F1_REP3
# TRINITY_DN100115_c0_g1     99.68891    370.12500    240.19452
# TRINITY_DN100154_c0_g1    833.49675    633.80126    890.89787
# TRINITY_DN100163_c0_g1    164.15827    176.85745    128.37936
# TRINITY_DN100287_c0_g1    948.84906    580.98588    918.25166
# TRINITY_DN100291_c0_g1     14.79638     13.38407     14.27775
# TRINITY_DN100299_c0_g1    948.48384    505.33286    726.54067
# HHHH_F1_REP4
# TRINITY_DN100115_c0_g1    264.49383
# TRINITY_DN100154_c0_g1    987.25237
# TRINITY_DN100163_c0_g1    146.72343
# TRINITY_DN100287_c0_g1    941.48736
# TRINITY_DN100291_c0_g1     19.08502
# TRINITY_DN100299_c0_g1    908.40217
dim(countsT)
# 24362    16 
#rounding the decimal output from salmon so that the output data works with DESeq2
RcountsT<-round(countsT) 
head(RcountsT)
#import the sample descript table
conds<- read.table("RT_tonsa_F1_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
#conds<- ddply(conds, "index", transform,  )
head(conds)
#how many reads from each sample
mean(colSums(RcountsT))
# 18166143
barplot(colSums(RcountsT), names.arg=colnames(RcountsT), cex.names=0.5, las=3) 
abline(h=mean(colSums(RcountsT)), col="blue", lwd=2) #shows the mean
mean(rowSums(RcountsT)) #11930.81
median(rowSums(RcountsT)) #2226
apply(RcountsT,2,mean) #rows
apply(RcountsT,1,mean) #columns
hist(apply(RcountsT,1,mean), xlim=c(0,2000), breaks=10000)
#negative binomial (counts) GLM
dds<-DESeqDataSetFromMatrix(countData=RcountsT, colData=conds,
                            design=~line+environment+line:environment)
dim(dds)
#filtering with too few reads
#24362    16
#this overwrites the original object
dds<-dds[rowSums(counts(dds))>160]
dim(dds)
#24362    16
dds<-DESeq(dds)
resultsNames(dds)
#[1] "Intercept"                  "line_Combined_vs_Ambient" 
#[3] "environment_HH_vs_AA"       "lineCombined.environmentHH"
#PCA to look at patterns of expression
#normalize the data
vsd<-vst(dds,blind=FALSE)
#save pca output as a dataset
data<-plotPCA(vsd,intgroup=c("line","environment"), returnData=TRUE)
percentVar<-round(100*attr(data,"percentVar"))
ggplot(data, aes(PC1,PC2,color=environment,shape=line))+
  geom_point(size=4,alpha=0.85)+
  xlab(paste0('PC1:',percentVar[1],"% Variance"))+
  ylab(paste0('PC2:',percentVar[2],"% Variance"))
#order and summarize results from specific contrasts
resIxn<-results(dds,alpha=0.05)
resIxn<-resIxn[order(resIxn$padj),]
head(resIxn)
# log2 fold change (MLE): lineCombined.environmentHH 
# Wald test p-value: lineCombined.environmentHH 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat
# <numeric>      <numeric> <numeric> <numeric>
#   TRINITY_DN115950_c0_g1   2245.97        4.05237  0.358490   11.3040
# TRINITY_DN131561_c0_g1   3375.99        4.64570  0.439847   10.5621
# TRINITY_DN137662_c0_g1  16743.23        4.90200  0.474583   10.3291
# TRINITY_DN149842_c8_g4  25971.82        4.27274  0.420809   10.1536
# TRINITY_DN129565_c0_g3  24258.76        4.30553  0.426037   10.1060
# TRINITY_DN129401_c0_g5  11712.31        4.46355  0.446094   10.0059
# pvalue        padj
# <numeric>   <numeric>
#   TRINITY_DN115950_c0_g1 1.25396e-29 2.99446e-25
# TRINITY_DN131561_c0_g1 4.46620e-26 5.33264e-22
# TRINITY_DN137662_c0_g1 5.20658e-25 4.14444e-21
# TRINITY_DN149842_c8_g4 3.19275e-24 1.90607e-20
# TRINITY_DN129565_c0_g3 5.19661e-24 2.48190e-20
# TRINITY_DN129401_c0_g5 1.43650e-23 5.71728e-20
summary(resIxn)
# out of 24362 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2839, 12%
# LFC < 0 (down)     : 1053, 4.3%
# outliers [1]       : 9, 0.037%
# low counts [2]     : 473, 1.9%
# (mean count < 18)
#######################
############################################## TEST FOR EFFECT OF ENVIRONMENT
####################### Likelhood ratio test

dds <- DESeqDataSetFromMatrix(countData = RcountsT, colData = conds, 
                              design = ~ line + environment)

dds <- DESeq(dds, test="LRT", reduced=~line)
# List the results you've generated
resultsNames(dds)
# [1] "Intercept"                "line_Combined_vs_Ambient"
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
# out of 24353 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3290, 14%
# LFC < 0 (down)     : 1028, 4.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 12)

resEnv <- resEnv[!is.na(resEnv$padj),]

degsEnv <- row.names(resEnv[resEnv$padj < 0.05,]) 
summary(resEnv)
summary(degsEnv)
#######################
##############################################  TEST FOR EFFECT OF LINE
#######################

dds <- DESeqDataSetFromMatrix(countData = RcountsT, colData = conds, 
                              design = ~ environment + line)

dds <- DESeq(dds, test="LRT", reduced=~environment)
resultsNames(dds)

resLine <- results(dds, alpha = 0.05)
resLine <- resLine[order(resLine$padj),]
head(resLine)

summary(resLine)
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2996, 13%
# LFC < 0 (down)     : 1180, 5%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 20)

resLine <- resLine[!is.na(resLine$padj),]

degsline <- row.names(resLine[resLine$padj < 0.05,])


#######################
##############################################  TEST FOR INTERACTION
#######################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, 
                              design = ~ environment + line + environment:line)

dds <- DESeq(dds, test="LRT", reduced=~environment + line)
resultsNames(dds)

resInt <- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),]
head(resInt)

summary(resInt)


resInt <- resInt[!is.na(resInt$padj),]

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
#TRINITY_DN115950_c0_g1
#TRINITY_DN138549_c1_g2
#TRINITY_DN131561_c0_g1 
d <-plotCounts(dds, gene="TRINITY_DN131561_c0_g1", intgroup = (c("line","environment")), returnData=TRUE)
d
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
lde<-length(degsEnv)  # 448
ldl<-length(degsline)  # 226
ldi<-length(degsInt)  # 3854


# Intersections
int_el<-length(intersect(degsEnv,degsline))  # 37
int_ei<-length(intersect(degsEnv,degsInt))  # 44
int_il<-length(intersect(degsInt,degsline))  # 34

intEL <- intersect(degsEnv,degsline)
mid<-length(intersect(degsInt,intEL)) # 7

# Number unique
env<-lde-int_ei-int_el-mid
line<-ldl-int_el-int_ei-mid
int<-ldi-int_ei-int_el-mid
#448-44-37-7 # 360
#226-37-34-7 # 148
#3854-44-34-7 # 3769

fit1 <- euler(c("Env" = 360, "Line" = 148, "Interaction" = 3769, "Env&Line" = 37, "Env&Interaction" = 44, "Line&Interaction" = 34, "Env&Line&Interaction" = 7))

plot(fit1,  lty = 1:3, quantities = TRUE)

plot(fit1, quantities = TRUE, fill = "transparent",
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