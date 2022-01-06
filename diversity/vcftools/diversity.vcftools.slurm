#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --cpus-per-task=70 
#SBATCH --mem=200G 
#SBATCH --time=2-02:30:00
#SBATCH --output=job.snpeff.%J.out 
#SBATCH --error=job.snpeff.%J.err
#SBATCH -J "snpeff" 
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu 

ml vcftools/0.1.14-xp36ajy


vcftools --vcf both.snpeff.vcf --out both.piDiversity --window-pi 100000 --window-pi-step 20000
vcftools --vcf A1.snpeff.vcf --out A1.piDiversity --window-pi 100000 --window-pi-step 20000
vcftools --vcf A2.snpeff.vcf --out A2.piDiversity --window-pi 100000 --window-pi-step 20000

vcftools --vcf both.snpeff.vcf --weir-fst-pop keep.A1 --weir-fst-pop keep.A2 --fst-window-size 100000 --fst-window-step 20000 --out A1vA2.fst

ml bedtools2
bedtools intersect -a both.snpeff.vcf -b exon.gtf -wa > both.snpeff.exon.vcf
head -n46 both.snpeff.vcf > header.both.snpeff.vcf
cat header.both.snpeff.vcf both.snpeff.exon.vcf > both.snpeff.exon.header.vcf
vcftools --vcf both.snpeff.exon.header.vcf --out both.exon.piDiversity --window-pi 100000 --window-pi-step 20000

bedtools intersect -a both.snpeff.vcf -b longicalyx. -wa > both.snpeff.intron.vcf
cat header.both.snpeff.vcf both.snpeff.intron.vcf > both.snpeff.intron.header.vcf
vcftools --vcf both.snpeff.intron.header.vcf --out both.intron.piDiversity --window-pi 100000 --window-pi-step 20000

zcat ../longicalyx.gff3.gz | grep gene | cut -f1,4,5,9 | sed -e 's/;.*//g' -e 's/-.*//g'| sort | uniq > longicalyx.gene.bed
bedtools subtract -a both.snpeff.vcf -b longicalyx.gene.bed > both.snpeff.intergenic.vcf
cat header.both.snpeff.vcf both.snpeff.intergenic.vcf > both.snpeff.intergenic.header.vcf
vcftools --vcf both.snpeff.intergenic.header.vcf --out both.intergenic.piDiversity --window-pi 100000 --window-pi-step 20000








bedtools intersect -a A1.snpeff.vcf -b exon.gtf -wa > A1.snpeff.exon.vcf
head -n46 A1.snpeff.vcf > header.A1.snpeff.vcf
cat header.A1.snpeff.vcf A1.snpeff.exon.vcf > A1.snpeff.exon.header.vcf
vcftools --vcf A1.snpeff.exon.header.vcf --out A1.exon.piDiversity --window-pi 100000 --window-pi-step 20000

bedtools intersect -a A1.snpeff.vcf -b longicalyx.intron.bed -wa > A1.snpeff.intron.vcf
cat header.A1.snpeff.vcf A1.snpeff.intron.vcf > A1.snpeff.intron.header.vcf
vcftools --vcf A1.snpeff.intron.header.vcf --out A1.intron.piDiversity --window-pi 100000 --window-pi-step 20000

bedtools subtract -a A1.snpeff.vcf -b longicalyx.gene.bed > A1.snpeff.intergenic.vcf
cat header.A1.snpeff.vcf A1.snpeff.intergenic.vcf > A1.snpeff.intergenic.header.vcf
vcftools --vcf A1.snpeff.intergenic.header.vcf --out A1.intergenic.piDiversity --window-pi 100000 --window-pi-step 20000



bedtools intersect -a A2.snpeff.vcf -b exon.gtf -wa > A2.snpeff.exon.vcf
head -n46 A2.snpeff.vcf > header.A2.snpeff.vcf
cat header.A2.snpeff.vcf A2.snpeff.exon.vcf > A2.snpeff.exon.header.vcf
vcftools --vcf A2.snpeff.exon.header.vcf --out A2.exon.piDiversity --window-pi 100000 --window-pi-step 20000

bedtools intersect -a A2.snpeff.vcf -b longicalyx.intron.bed -wa > A2.snpeff.intron.vcf
cat header.A2.snpeff.vcf A2.snpeff.intron.vcf > A2.snpeff.intron.header.vcf
vcftools --vcf A2.snpeff.intron.header.vcf --out A2.intron.piDiversity --window-pi 100000 --window-pi-step 20000

bedtools subtract -a A2.snpeff.vcf -b longicalyx.gene.bed > A2.snpeff.intergenic.vcf
cat header.A2.snpeff.vcf A2.snpeff.intergenic.vcf > A2.snpeff.intergenic.header.vcf
vcftools --vcf A2.snpeff.intergenic.header.vcf --out A2.intergenic.piDiversity --window-pi 100000 --window-pi-step 20000





