!#/bin/bash
cd /data/project_data/PopGenomics/shared/Chr19_sweeps/LAI
#run loter
#open another screen bc it takes a long time
# First, make a new dir to store the results:
mkdir Loter_out
screen
CHR="Chr19"
#bash script
while read ID
do
for file in Admixed/poplar_hybrids.maf05.${CHR}.${ID}.vcf.gz
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

#ctrl c to kill process





 