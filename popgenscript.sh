#!/bin/bash

# Rename value in <> to your chromosome number! Include the zero only if your chromosome # is <10

myChr=Chr19

# myChr=Chr19  # go to the right directory

cd /data/project_data/PopGenomics

# Run VCFtools to subset the big vcf file for just your chromosome

vcftools --gzvcf poplar_hybrids.maf05.vcf.gz \
--chr $myChr \
--out shared/$myChr \
--recode
#kept 207584 our of 4803704
# Extract the centromere coordinates for your chromosome so you can exclude those regions from your sweep analysis

grep $myChr poplar_centromeres.txt > shared/${myChr}_centromere.txt # grab the centromere location for your chromosome

cd shared/

mkdir ${myChr}_sweeps  # make a new directory for your chromosome analyses

mv *${myChr}* ${myChr}_sweeps # clean up the space by moving all files into your the directory you just made

cd ${myChr}_sweeps

# Test for selective sweeps

RAiSD -n $myChr \
-I data/project_data/PopGenomics/shared/${myChr}.recode.vcf \
-f -t -R -P -D -A 0.99 \
-X ${myChr}_centromere.txt

# Estimate nucleotide diversity (pi) in sliding windows of 50kb

vcftools --vcf /data/project_data/PopGenomics/shared/${myChr}.recode.vcf \
--chr $myChr \
--window-pi 50000 \
--out $myChr

# First, need to subset the metadata file for just those individuals with balsamifera ancestry
# We can do this using an interactive R session at the commandline. 
# An alternative is to put these R commands in a script, save it with the ".r" extension, 
# and at the commandline type "Rscript myscript.r"

R # Opens an interactive R session within Unix...
Qscores <- read.table("../poplar_hybrids.LDpruned.5.Q", sep=" ",header=F)
names(Qscores) = c("K1","K2","K3","K4","K5")

meta <- read.table("../../Combined_Transect_Sampling_Data_2020.txt",sep="\t",header=T)

merged <- cbind(meta,Qscores)
str(merged)

Bals_Inds <- merged[which(merged$K4>0.5),1]  
length(Bals_Inds) # Should net you 188 individuals

Tricho_Inds <- merged[which(merged$K4<=0.5),1]
length(Tricho_Inds) # Should net you 388 individuals

# Write out your Bals and Tricho lists as tab-delimited text files
write.table(Bals_Inds, "Bals_Inds.txt", quote=F, row.names=F, col.names=F)

write.table(Tricho_Inds, "Tricho_Inds.txt", quote=F, row.names=F, col.names=F)

quit()

# When prompted with: "Save workspace image? [y/n/c]"  choose: n

vcftools --vcf /data/project_data/PopGenomics/shared/${myChr}.recode.vcf \
--weir-fst-pop Bals_Inds.txt \
--weir-fst-pop Tricho_Inds.txt \
--fst-window-size 50000 \
--out Bals_Tricho_All

#11/8/2021
#vim looking at loter output
#get out of vim, :q
#2 chunks of data, one chunk from B, one from T
#Our pipeline:

#1Convert haploid 0/1 calls to diploid 0/1/2 calls
#2Collate 0/1/2 data across all admixed individuals into a single data matrix
#3Convert matrix to Plink2.0 format (*.bed file)
#4Estimate frequencies of local ancestry along chromosomes
#setting variables such as the Chromosome ID, how many SNP positions we’re analyzing on that chromosome, and the number of individuals.
## As LOTER files get output, can convert these into 0/1/2 encoding using the bash tool 'datamash'
# steps 1 and 2
## From within your LAI/ directory:

CHR="Chr19"

echo $CHR # Is it right?

Nsites=`tail -n +2 ${CHR}.kept.sites | wc -l | sed 's/\s/\t/' | cut -f1` # calculates and stores the number of SNP sites for your chromosome

echo $Nsites # For Chr19, 207584 sites

Ninds=`wc -l Admixed.Inds | sed 's/\s/\t/' | cut -f1` # calculates and stores the number of admixed individuals you previously identified

echo $Ninds  # Should be 442 individuals
touch ${CHR}_matrix.diploid

for file in Loter_out/*.txt
do
datamash --field-separator=" " sum 1-${Nsites} <$file >>${CHR}_matrix.diploid
done
sed 's/\s/\t/g' ${CHR}_matrix.diploid | cut -f1-${Nsites} | datamash transpose >${CHR}_matrix.diploid.tr
#(3) Convert matrix to Plink2.0 format
seq -f "snp%02g" 1 $Nsites >sites

printf 'A\n%.0s' $(seq $Nsites) >allele1  # Create a dummy column of 'A' the length of your Nsites file

printf "T\n%.0s" $(seq $Nsites) >allele2 # Create a dummy column of 'T' the length of your Nsites file

mkdir Plink

paste sites allele1 allele2 ${CHR}_matrix.diploid.tr >Plink/${CHR}_matrix.diploid.tr.forPlink
#This step will make the fam file for our samples:
cat /data/project_data/PopGenomics/Combined_Transect_Sampling_Data_2020.txt | \
cut -f1-2 | \
grep -w -f Admixed.Inds - | \
cut -f2 | \
paste - Admixed.Inds >FID_IID

printf '0\t0\t0\t-9\n%.0s' $(seq $Ninds) >dummy

paste FID_IID dummy  >Plink/${CHR}_fam.forPlink
#Now, we’re ready to run the conversion into Plink format (.bed)
## This runs the Plink conversion from allele dosages to bed format

cd Plink/ 

plink2 --import-dosage ${CHR}_matrix.diploid.tr.forPlink noheader \
--fam ${CHR}_fam.forPlink \
--make-bed \
--out ${CHR}_Admixed_FAI
#we’ll simply have Plink calculate the local ancestry frequencies at each SNP position:
plink2 --bfile ${CHR}_Admixed_FAI --freq --out ${CHR}_LAI_freq
#You’ll want 2 files for this (replacing “ChrXX” with your Chr ID):
#scp
#ChrXX_SNP.kept.sites
#ChrXX_LAI_freq.afreq
#assumptions:
#hybridization, mosaic ancestry, phenotypic variation differing among them
#pipeline for admixture mapping
#Wrangle the phenotype data from the common gardens to get just the data for the admixed individuals we ran Loter on
#Get genome-wide ancestry for those same individuals (based off the K=5 model, subsetting for just the balsamifera Q value)
#Run Admixture mapping in Plink
#Bring the results into R and visualize


#getting phenotype data ready
CHR="Chr19"

echo $CHR  # Does it look right? --yes

cd /data/project_data/PopGenomics/shared/${CHR}_sweeps/LAI/
#transect, ID, rust, infection rate, bud flush, bud set

# Get phenotype data from meta file for the admixed samples:

tail -n +2 /data/project_data/PopGenomics/VT_Garden_Phenotypes_2021.txt | \
grep -w -f Admixed.Inds - >Admixed.pheno

printf "#FID\tIID\tRUST\tFLUSH\tSET\n%.0s" >Plink/pheno.forPlink

cat Admixed.pheno >>Plink/pheno.forPlink
#step 2
# Get K=2 ADMIX to use as covariate; Need to use tail -n +2 to skip the header in the metadata file before pasting to the Q file

tail -n +2 /data/project_data/PopGenomics/Combined_Transect_Sampling_Data_2020.txt | \
cut -f1 | \
paste - /data/project_data/PopGenomics/shared/poplar_hybrids.LDpruned.5.Q | \
grep -w -f Admixed.Inds - | \
sed 's/\s/\t/g' | \
cut -f 1,5 >Plink/Admixed_KBals

cat /data/project_data/PopGenomics/Combined_Transect_Sampling_Data_2020.txt | \
cut -f1-2 | \
grep -w -f Admixed.Inds - | \
sed 's/\s/\t/g' | \
cut -f2 >Plink/Transect
#would just choose the least common as default, need to specify the reference group
# Create the cov file with KBals

printf "#FID\tIID\tK2Q\n%.0s" >Plink/cov.forPlink
paste Plink/Transect Plink/Admixed_KBals >>Plink/cov.forPlink
#step 3
cd Plink/

plink2 --bfile ${CHR}_Admixed_FAI \
--pheno pheno.forPlink \
--covar cov.forPlink \
--glm omit-ref
#within plink
#-rw-r--r--. 1 bchrist4 users 36400787 Nov  9 19:01 plink2.FLUSH.glm.linear
#-rw-r--r--. 1 bchrist4 users 37591739 Nov  9 19:01 plink2.RUST.glm.linear
#-rw-r--r--. 1 bchrist4 users 35911442 Nov  9 19:01 plink2.SET.glm.linear
#needed files for scp
# Meta-data:
/data/project_data/PopGenomics/Combined_Transect_Sampling_Data_2020.txt

# Phenotypes:
/data/project_data/PopGenomics/VT_Garden_Phenotypes_2021.txt

# Admixed.Inds:
/data/project_data/PopGenomics/shared/${CHR}_sweeps/LAI/Admixed.Inds

# SNPs analyzed for your chromosome:
/data/project_data/PopGenomics/shared/${CHR}_sweeps/LAI/${CHR}.kept.sites

# Frequencies of Local Ancestry:
/data/project_data/PopGenomics/shared/${CHR}_sweeps/LAI/Plink/${CHR}_LAI_freq.afreq

# Genome-wide balsamifera ancestry:
/data/project_data/PopGenomics/shared/${CHR}_sweeps/LAI/Plink/Admixed_KBals