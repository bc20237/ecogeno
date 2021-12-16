vcftools --gzvcf poplar_hybrids.vcf.gz <OPTIONS>

vcftools --gzvcf poplar_hybrids.vcf.gz --maf <float> --chr <chromosome ID> --out ~/myresults/pi_chrX_mafX_winX

vcftools --gzvcf poplar_hybrids.vcf.gz --maf <float> --het --out ~/myresults/het_mafX

cd /data/project_data/PopGenomics

zcat poplar_hybrids.vcf.gz | head -n 11

##fileformat=VCFv4.2
##filedate=20210624
##source="beagle.27Jul16.86a.jar (version 4.1)"
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated ALT Allele Frequencies">
##INFO=<ID=AR2,Number=1,Type=Float,Description="Allelic R-Squared: estimated squared correlation between most probable REF dose and true REF dose">
##INFO=<ID=DR2,Number=1,Type=Float,Description="Dosage R-Squared: estimated squared correlation between estimated REF dose [P(RA) + 2*P(RR)] and true REF dose">
##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="estimated ALT dose [P(RA) + P(AA)]">
##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Genotype Probability">
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  201 202  ......

# Run Admixture 

K=<your value>

admixture -j3 --cv $FILE2.bed $K >log${K}.out

grep "cv" logKx.out #--> cross validation and the accuracy of the model

#individual heterozygosity (x axis for our plot of diversity of ancestry scores v genome wide heterozygosity)
vcftools --gzvcf /data/project_data/PopGenomics/poplar_hybrids.maf05.vcf.gz \
--het \
--out ~/myresults/het_maf05
#copy this to my 
scp bchrist4@pbio381.uvm.edu:~/myresults/het_maf05 .

#cp the k=5 file 

# Rename value in <> to your chromosome number! Include the zero only if your chromosome # is <10

myChr=<Chr0X>  

# myChr=Chr02  # I used Chr02 for my test run...


cd /data/project_data/PopGenomics

# Run VCFtools to subset the big vcf file for just your chromosome

vcftools --gzvcf poplar_hybrids.maf05.vcf.gz \
--chr $myChr \
--out shared/$myChr \
--recode
#install python and loter package on local server
git clone https://github.com/bcm-uga/Loter.git
cd Loter/python-package
python setup.py install --user

cd /data/project_data/PopGenomics/shared/Chr19_sweeps/