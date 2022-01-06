#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --time=4-02:30:00
#SBATCH --output=job.pixy.%J.out
#SBATCH --error=job.pixy.%J.err
#SBATCH --partition=speedy,biocrunch
#SBATCH -J "pixy"
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu

ref=longicalyx.fasta
ext=pixy
thr=20
name=all.samples

module load sentieon-genomics
sentieon driver -t $thr -r $ref --algo GVCFtyper --emit_mode all $name.pixy.vcf gvcf/*.gVCF

ml vcftools bcftools
vcftools --vcf $name.$ext.vcf --remove-indels --max-missing-count 12 --min-meanDP 10 --max-meanDP 100 --maf 0.05 --recode --recode-INFO-all --out $name.$ext.filtered
bgzip $name.$ext.filtered.recode.vcf
tabix $name.$ext.filtered.recode.vcf.gz

vcftools --vcf $name.$ext.vcf --max-maf 0 --recode --out $name.$ext.vcf.invariant
bgzip $name.$ext.invariant.recode.vcf
tabix $name.$ext.invariant.recode.vcf.gz

bcftools concat --allow-overlaps $name.$ext.invariant.recode.vcf.gz $name.$ext.filtered.recode.vcf.gz -O z -o $name.$ext.combined.vcf.gz
tabix $name.$ext.filtered.recode.vcf.gz
bcftools query -l $name.$ext.filtered.recode.vcf.gz > populations.bcftools

bcftools view -O z -o all.gene.vcf.gz -R longicalyx.gene.bed --threads 20 all.samples.pixy.combined.vcf.gz
bcftools view -O z -o all.exon.vcf.gz -R longicalyx.exon.bed --threads 20 all.samples.pixy.combined.vcf.gz
bcftools view -O z -o all.intron.vcf.gz -R longicalyx.intron.bed --threads 20 all.samples.pixy.combined.vcf.gz
bcftools view -O z -o all.intergenic.vcf.gz -R longicalyx.intergenic.bed --threads 20 all.samples.pixy.combined.vcf.gz

tabix all.gene.vcf.gz 
tabix all.exon.vcf.gz 
tabix all.intron.vcf.gz  
tabix all.intergenic.vcf.gz

module purge

ml miniconda3
source activate pixy
pixy --stats pi fst dxy --vcf all.samples.pixy.combined.vcf.gz --populations populations.bcftools --window_size 10000 --n_cores 20 --output_prefix all.total
pixy --stats pi fst dxy --vcf all.gene.vcf.gz --populations populations.bcftools --window_size 10000 --n_cores 20 --output_prefix all.gene
pixy --stats pi fst dxy --vcf all.exon.vcf.gz --populations populations.bcftools --window_size 10000 --n_cores 20 --output_prefix all.exon
pixy --stats pi fst dxy --vcf all.intron.vcf.gz --populations populations.bcftools --window_size 10000 --n_cores 20 --output_prefix all.intron
pixy --stats pi fst dxy --vcf all.intergenic.vcf.gz --populations populations.bcftools --window_size 10000 --n_cores 20 --output_prefix all.intergenic

