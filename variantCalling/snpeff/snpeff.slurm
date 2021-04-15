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

module load snpeff/2017-11-24-3ueyysf

# gather files
mkdir snpeff && cd /work/LAS/jfw-lab/corrinne/Adom/F.reference/snpeff/
cp ../longicalyx.fasta .
cp ../longicalyx.gff3 .
gzip longicalyx.fasta

# fix the gff3 names and make the gtf
sed 's/,maker[^;,]*//g' longicalyx.gff3 > longicalyx.fixed.gff3

ml genometools/1.5.9-py2-awelg6n
gt gff3_to_gtf longicalyx.gff3 > longicalyx.gtf
gzip longicalyx.gtf

# make the snpeff database for the longicalyx genome
mkdir longicalyx
cp longicalyx.gtf.gz longicalyx/genes.gtf.gz
cp longicalyx.fasta.gz longicalyx/sequences.fa.gz
echo "longicalyx.genome : longicalyx" >> snpEff.config

snpEff build -v longicalyx -gtf22 -dataDir . -c snpEff.config

# annotate effects for SNPs using snpeff
snpEff ann -t -dataDir . -c snpEff.config longicalyx ../all.samples.F1.SNPs.good.10pct.filtered.recode.vcf > snpeff.vcf
snpEff ann -t -dataDir . -no-downstream -no-upstream -no-intergenic -c snpEff.config longicalyx ../all.samples.F1.SNPs.good.10pct.filtered.recode.vcf > snpeff.geneOnly.vcf

############# get only certain features from gff and use those to prune the vcf
############zcat ../longicalyx.gff3.gz | awk ' { OFS="\t" }  ( $3 ~ "CDS" || $3 ~ "UTR" )  { print $1,$4,$5,$3 }' > longicalyx.utrcds.bed
############zcat ../longicalyx.gff3.gz | awk ' { OFS="\t" }  $3 ~ "mRNA"  { print $1,$4,$5,$3 }' > longicalyx.mRNA.bed
############zcat ../longicalyx.gff3.gz | awk ' { OFS="\t" }  ( $3 ~ "CDS" )  { print $1,$4,$5,$3 }' > longicalyx.exon.bed
############
############bedtools subtract -a longicalyx.mRNA.bed -b longicalyx.utrcds.bed > longicalyx.intron.bed
############
############bedtools intersect -wa -a snpeff.vcf -b longicalyx.exon.bed > exons.snps.gff
############bedtools intersect -wa -a snpeff.vcf -b longicalyx.intron.bed > intron.snps.gff

# extract each species and run snpeff
ml bcftools
bcftools query -l ../all.samples.F1.SNPs.good.10pct.filtered.recode.vcf | grep -v "Gar\|A2\|6177" > keep.A1
bcftools query -l ../all.samples.F1.SNPs.good.10pct.filtered.recode.vcf | grep -v "Ghe\|A1\|6177" > keep.A2

vcftools --vcf ../all.samples.F1.SNPs.good.10pct.filtered.recode.vcf --out A1 --keep keep.A1 --recode
vcftools --vcf ../all.samples.F1.SNPs.good.10pct.filtered.recode.vcf --out A2 --keep keep.A2 --recode
vcftools --vcf ../all.samples.F1.SNPs.good.10pct.filtered.recode.vcf --out both --keep keep.both --recode

vcftools --vcf A1.recode.vcf --out A1.SNP.filtered --max-non-ref-af 0.99 --recode
vcftools --vcf A2.recode.vcf --out A2.SNP.filtered --max-non-ref-af 0.99 --recode
vcftools --vcf both.recode.vcf --out both.SNP.filtered --max-non-ref-af 0.99 --recode

snpEff ann -dataDir . -stats A1.html -c snpEff.config longicalyx A1.SNP.filtered.recode.vcf > A1.snpeff.vcf
snpEff ann -dataDir . -stats A2.html -c snpEff.config longicalyx A2.SNP.filtered.recode.vcf > A2.snpeff.vcf
snpEff ann -dataDir . -stats both.html -c snpEff.config longicalyx both.SNP.filtered.recode.vcf > both.snpeff.vcf

cat both.snpeff.vcf | SnpSift caseControl "+++++++++++++++++++-------------------------------------------------------------------------------------+-+-------------" > both.snpeff.cc.vcf
cat both.snpeff.cc.vcf | SnpSift filter "(Cases[0] = 21) & (Controls[0] = 0)" > both.filtered.A1.fixed.vcf
cat both.snpeff.cc.vcf | SnpSift filter "(Cases[0] = 0) & (Controls[0] = 99)" > both.filtered.A2.fixed.vcf

cat both.snpeff.cc.vcf | SnpSift filter "(Cases[0] = 0) & ((Controls[0] +  Controls[1]) > 10)" | wc -l
cat both.snpeff.cc.vcf | SnpSift filter "(Controls[0] = 0) & ((Cases[0] +  Cases[1]) > 2)" | wc -l



# Approximately 12M SNP sites for all A accessions (A1 and A2); of these, 
# •	Individual samples are homozygous for reference (F1) allele at  5.6 – 7.9 M sites
# •	Individual samples are homozygous for the derived allele at 1.6 – 2 M sites
# •	Individual samples are heterozygous for the derived allele at 1.7 – 4.8 M sites
# 
# In this general background, we see
# •	352,691 sites where the derived allele is fixed in A1 and absent in A2
# 	o	4,042,689 sites where the derived allele is in at least 10% of the A1 population (het, or hom derived, 2 out of 21 samples), but absent in A2
# •	99,261 sites where the derived allele is fixed in A2 and absent in A1
# 	o	7,649,713 sites where the derived allele is in at least 10% of the A2 population (het, or hom derived, 10 out of 99 samples), but absent in A1
# 




# get each sample in an individual SNP vcf
cat keep.A* | while read line; do vcftools --vcf ../all.samples.F1.SNPs.good.10pct.filtered.recode.vcf --out $line.SNP --indv $line --max-non-ref-af 0.99 --recode; done


##### Indels
# extract each species 
vcftools --vcf ../all.samples.F1.indels.good.10pct.recode.vcf --out A1.indels --keep keep.A1 --recode
vcftools --vcf ../all.samples.F1.indels.good.10pct.recode.vcf --out A2.indels --keep keep.A2 --recode

# remove lines that are reference only
vcftools --vcf A1.indels.recode.vcf --out A1.indels.filtered --max-non-ref-af 0.99 --recode
vcftools --vcf A2.indels.recode.vcf --out A2.indels.filtered --max-non-ref-af 0.99 --recode

# annotate as ins or deletion
SnpSift varType ../all.samples.F1.indels.good.10pct.recode.vcf > both.indel.insVdel.vcf
SnpSift varType A1.indels.filtered.recode.vcf > A1.indel.insVdel.vcf
SnpSift varType A2.indels.filtered.recode.vcf > A2.indel.insVdel.vcf

# generate tables to get general numbers
grep -v '#' both.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep INS > both.indel.INS.tbl
grep -v '#' both.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep DEL > both.indel.DEL.tbl
grep -v '#' both.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep INS | grep -v DEL > both.indel.INSOnly.tbl
grep -v '#' both.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep DEL | grep -v INS > both.indel.DELOnly.tbl
grep -v '#' both.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep INS | grep DEL > both.indel.INSDEL.tbl

grep -v '#' A1.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep INS > A1.indel.INS.tbl
grep -v '#' A1.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep DEL > A1.indel.DEL.tbl
grep -v '#' A1.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep INS | grep -v DEL > A1.indel.INSOnly.tbl
grep -v '#' A1.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep DEL | grep -v INS > A1.indel.DELOnly.tbl
grep -v '#' A1.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep INS | grep DEL > A1.indel.INSDEL.tbl

grep -v '#' A2.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep INS > A2.indel.INS.tbl
grep -v '#' A2.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep DEL > A2.indel.DEL.tbl
grep -v '#' A2.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep INS | grep -v DEL > A2.indel.INSOnly.tbl
grep -v '#' A2.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep DEL | grep -v INS > A2.indel.DELOnly.tbl
grep -v '#' A2.indel.insVdel.vcf | cut -f1,2,4,5,8 | grep INS | grep DEL > A2.indel.INSDEL.tbl

# only get the simple insertion and deletions from the file
grep -v "INS,\|DEL,\|,DEL\|,INS\|SNP,\|,SNP\|,INTERVAL\|INTERVAL," A1.indel.insVdel.vcf > A1.indel.insVdelSimple.vcf
grep -v "INS,\|DEL,\|,DEL\|,INS\|SNP,\|,SNP\|,INTERVAL\|INTERVAL," A2.indel.insVdel.vcf > A2.indel.insVdelSimple.vcf
grep -v "INS,\|DEL,\|,DEL\|,INS\|SNP,\|,SNP\|,INTERVAL\|INTERVAL," both.indel.insVdel.vcf > both.indel.insVdelSimple.vcf

# annotate effects of simple insertion and deletion sites
snpEff ann -dataDir . -stats A1.indel.html -c snpEff.config longicalyx A1.indel.insVdelSimple.vcf > A1.indel.snpeff.vcf
snpEff ann -dataDir . -stats A2.indel.html -c snpEff.config longicalyx A2.indel.insVdelSimple.vcf > A2.indel.snpeff.vcf
snpEff ann -dataDir . -stats both.indel.html -c snpEff.config longicalyx both.indel.insVdelSimple.vcf > both.indel.snpeff.vcf

grep "frameshift" both.indel.snpeff.vcf | awk ' {OFS="\t"} {print $1,$2,$8} ' | sort | uniq > frameshift.indels.bed


# get individual samples 
cat keep.A* | while read line; do vcftools --vcf ../all.samples.F1.indels.good.10pct.recode.vcf --out $line.indels --indv $line --recode; done
