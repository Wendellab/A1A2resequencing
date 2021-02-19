#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --cpus-per-task=70 
#SBATCH --mem=200G 
#SBATCH --time=7-02:30:00
#SBATCH --output=job.Fvcf.%J.out 
#SBATCH --error=job.Fvcf.%J.err
#SBATCH -J "F.jointGeno" 
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu 

ref=longicalyx.fasta
ext=F1
thr=70
name=all.samples
missing=37

mkdir -p sam sortBam dedupBam realignBam gvcf joblog metrics timelog
mv *.sam sam/
mv job.* joblog
mv *.sort.bam* sortBam
mv *.dedup.bam* dedupBam
mv *.score* metrics
mv *.metric metrics
mv *.timelog timelog

module load sentieon-genomics
time ( sentieon driver -t $thr -r $ref --algo GVCFtyper $name.$ext.vcf *.gVCF ) &> $name.$ext.vcf.timelog

mv *.gVCF* gvcf

ml vcftools
vcftools --vcf $name.$ext.vcf --remove-indels --recode --recode-INFO-all --out $name.$ext.SNPs
vcftools --vcf $name.$ext.vcf --keep-only-indels --recode --recode-INFO-all --out $name.$ext.indels
vcftools --vcf $$name.$ext.SNPs.recode.vcf --out $name.$ext.SNPs.10pct --max-missing-count 37 --recode
vcftools --vcf $name.$ext.SNPs.10pct.recode.vcf --out $name.$ext.SNPs.10pct --depth 

module load snphylo
sed -e 's/^F0//g' -e 's/^F//g' $name.$ext.SNPs.10pct.recode.vcf > $name.$ext.numOnly.vcf
snphylo.sh -v $name.$ext.numOnly.vcf -c 5 -P $name.$ext.phylo -b -B 10000


module load r-devtools/1.12.0-py2-r3.4-47jzlan
module load r-gdsfmt/1.14.1-py2-r3.4-s3s72d7
module load r-snprelate/1.12.2-py2-r3.4-dzjxl6t
module load r-stringr/1.2.0-py2-r3.4-ncynzzd
module load r-ggplot2/2.2.1-py2-r3.4-l3jmw7u
module load r-dplyr/0.7.3-py2-r3.4-p2gyosx

Rscript PCAforSNP.r $name.$ext.SNPs.10pct.recode.vcf

ml 
