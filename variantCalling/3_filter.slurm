#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --cpus-per-task=20 
#SBATCH --mem=200G 
#SBATCH --time=7-02:30:00
#SBATCH --output=job.Ffilter.%J.out 
#SBATCH --error=job.Ffilter.%J.err
#SBATCH -J "F.filter" 
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu 

ext=F1
name=all.samples

ml bcftools
# paste <(bcftools query -l all.samples.F1.SNPs.recode.vcf | sed 's/ /\\ /g') <(bcftools query -l all.samples.F1.SNPs.recode.vcf | sed 's/ //g') > rename.samples
bcftools reheader --threads 20 -s rename.samples $name.$ext.SNPs.recode.vcf > $name.$ext.SNPs.rename.vcf
bcftools reheader --threads 20 -s rename.samples $name.$ext.indels.recode.vcf > $name.$ext.indels.rename.vcf

ml vcftools
vcftools --vcf $name.$ext.SNPs.rename.vcf --out $name.$ext.SNPs.good --keep keep.samples --recode
vcftools --vcf $name.$ext.indels.rename.vcf --out $name.$ext.indels.good --keep keep.samples --recode


vcftools --vcf $name.$ext.SNPs.good.recode.vcf --out $name.$ext.SNPs.good.10pct --max-missing-count 12 --recode
vcftools --vcf $name.$ext.SNPs.good.10pct.recode.vcf --out $name.$ext.SNPs.good.10pct --depth 

vcftools --vcf $name.$ext.indels.good.recode.vcf --out $name.$ext.indels.good.10pct --max-missing-count 12 --recode




vcftools --vcf $name.$ext.SNPs.good.10pct.recode.vcf --out $name.$ext.SNPs.good.10pct.filtered --max-alleles 2 --min-meanDP 10 --max-meanDP 100 --maf 0.05 --max-non-ref-af 0.99 --recode
vcftools --vcf $name.$ext.SNPs.good.10pct.filtered.recode.vcf --out $name.$ext.SNPs.good.4phylo --max-missing-count 0 --recode
#vcftools --vcf all.samples.F1.10pct.gt10.pruned.recode.vcf --out all.samples.noF1.4PCA --recode --remove remove.F1

ml bedtools2
bedtools intersect -a $name.$ext.SNPs.good.10pct.filtered.recode.vcf -b longicalyx.exon.gff3 >  $name.$ext.SNPs.geneOnly.vcf

cat <(head -n41 $name.$ext.SNPs.good.10pct.filtered.recode.vcf) $name.$ext.SNPs.geneOnly.vcf > $name.$ext.SNPs.geneOnly.fixed.vcf
vcftools --vcf $name.$ext.SNPs.geneOnly.fixed.vcf --out $name.$ext.SNPs.geneOnly.nohybrid --recode --remove remove.hybrid


module load snphylo
sed -i 's/>F/>Chr/g' $name.$ext.10pct.gt10.recode.vcf
sed 's/Chr//g' $name.$ext.10pct.gt10.recode.vcf > $name.$ext.filter1010.numOnly.vcf
snphylo.sh -v $name.$ext.filter1010.numOnly.vcf -c 5 -P $name.$ext.filter1010.phylo -b -B 10000

ml 
