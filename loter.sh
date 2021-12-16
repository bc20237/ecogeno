CHR="ChrXX"  # Be sure to customize to your chromosome number!!! 19

echo $CHR  # Does it look right?

# Then make some new folders to store future results in:

mkdir LAI

cd LAI/

mkdir Admixed
# Interactive R session at the commandline:

R

# Import K=5 Admixture run
Qscores <- read.table("/data/project_data/PopGenomics/shared/poplar_hybrids.LDpruned.5.Q", sep=" ",header=F)
names(Qscores) = c("K1","K2","K3","K4","K5")

# Import meta-data
meta <- read.table("/data/project_data/PopGenomics/Combined_Transect_Sampling_Data_2020.txt",sep="\t",header=T)
#1.
# Combine them 
merged <- cbind(meta,Qscores)
str(merged)

for(i in 1:nrow(merged)){
if(merged$K4[i]>=0.99){
    merged$Anc[i]="Bals"
    }else if (sum(c(merged$K1[i],merged$K2[i],merged$K5[i]))>=0.99){
    merged$Anc[i]="Tricho"
    } else if(merged$K3[i]>0.5){
    merged$Anc[i]="PopSpp"
    }else{
    merged$Anc[i]="Admx"
    }
}

table(merged$Anc)
# Admx   Bals PopSpp Tricho 
#   442     46      8     80 
#Now, we can write these sample IDs out to separate files to use with VCFtools to group our samples.
Bals_Inds_Ref <- merged[merged$Anc=="Bals",1]  
length(Bals_Inds_Ref) # Should net you 46 individuals
write.table(Bals_Inds_Ref, "Balsam.Inds", quote=F, row.names=F, col.names=F)

Tricho_Inds_Ref <- merged[merged$Anc=="Tricho",1]
length(Tricho_Inds_Ref) # Should net you 80 individuals
write.table(Tricho_Inds_Ref, "Tricho.Inds", quote=F, row.names=F, col.names=F)

Admixed_Inds <- merged[merged$Anc=="Admx",1]
length(Admixed_Inds) # Should net you 442 individuals
write.table(Admixed_Inds, "Admixed.Inds", quote=F, row.names=F, col.names=F)

quit() # choose 'n' when prompted

#2
#VCF tools --> parse out into individuals with high balsm ancestry and high tricho ancestry
# First we grab just the Balsam and Tricho reference individuals

vcftools --gzvcf /data/project_data/PopGenomics/shared/${CHR}_sweeps/${CHR}.recode.vcf --keep Balsam.Inds --recode --stdout | gzip -c >poplar_hybrids.maf05.${CHR}.BalsRef.vcf.gz

vcftools --gzvcf /data/project_data/PopGenomics/shared/${CHR}_sweeps/${CHR}.recode.vcf --keep Tricho.Inds --recode --stdout | gzip -c >poplar_hybrids.maf05.${CHR}.TrichoRef.vcf.gz
#We’ll also want to export a list of all the SNP positions in our VCF files for bringing into R later and merging with the LAI outputs
vcftools --gzvcf poplar_hybrids.maf05.${CHR}.BalsRef.vcf.gz --kept-sites --out ${CHR}
#go into screen to do the loter stuff
screen
CHR="Chr19"
while read ID
do
  vcftools --gzvcf /data/project_data/PopGenomics/shared/${CHR}_sweeps/${CHR}.recode.vcf \
  --indv $ID \
  --recode \
  --stdout | gzip -c   >Admixed/poplar_hybrids.maf05.${CHR}.${ID}.vcf.gz
done < Admixed.Inds
#reads in our master VCF file that is already subsetted by our chromosome of interest
#subsets by just the individual of interest (--indv)
#writes a new file based on the subsetting (--recode)
#passes that file to “stdout” which is used to send it to another program on the other side of the “pipe” |
#uses gzip to compress the new output file to save space
#uses the < Admixed.Inds to feed the loop sample ID’s of the admixed individuals, one line at a time.
#exit screen using ctrl a d
#run loter
#open another screen bc it takes a long time
# First, make a new dir to store the results:
mkdir Loter_out
CHR="Chr19"
#bash script
while read ID
do
for file in Admixed/poplar_hybrids.maf05.${CHR}.${ID}.vcf
do
loter_cli -r poplar_hybrids.maf05.${CHR}.BalsRef.vcf.gz poplar_hybrids.maf05.${CHR}.TrichoRef.vcf.gz \
-a $file \
-f vcf \
-o Loter_out/${ID}_LAI.txt \
-n 1 \
-pc -v
done
done < Admixed.Inds
#screen -r reattach to screen 
#-o directs the analysis output to the Loter_out directory you made, and names it by the ID name
#-n controls how many cpu’s to use in the analysis (just 1 for now)
#-pc implements phase correction on the inferred ancestry
#-v prints “verbose” progress updates to the screen as the analysis proceeds
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
grep -w -f Admixed.Inds - | \ #searching for ids that match those 442 admixed id', match the whole word and the search command from this input file
cut -f2 | \
paste - Admixed.Inds >FID_IID

printf '0\t0\t0\t-9\n%.0s' $(seq $Ninds) >dummy

paste FID_IID dummy  >Plink/${CHR}_fam.forPlink
#can head this to check out whats up
#transect, id, a few dummy variables, phenotype etc
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

