#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --cpus-per-task=7 
#SBATCH --mem=20G 
#SBATCH --time=4-02:30:00
#SBATCH --output=job.F.%J.out 
#SBATCH --error=job.F.%J.err
#SBATCH --job-name "F.prepare" 
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu 

DIR=/work/LAS/jfw-lab/corrinne/Adom/F.reference/genome
file=longicalyx.fasta

mkdir $DIR
cd $DIR

ml curl 

curl -O ftp://ftp.bioinfo.wsu.edu/species/Gossypium_longicalyx/F1_NSF/assembly/longicalyx.fasta.gz
curl -O ftp://ftp.bioinfo.wsu.edu/species/Gossypium_longicalyx/F1_NSF/genes/longicalyx.gff3.gz
curl -O ftp://ftp.bioinfo.wsu.edu/species/Gossypium_longicalyx/F1_NSF/genes/longicalyx.proteins.fa.gz
curl -O ftp://ftp.bioinfo.wsu.edu/species/Gossypium_longicalyx/F1_NSF/genes/longicalyx.transcripts.fa.gz
gunzip longicalyx.fasta.gz

ml samtools
ml bwa

bwa index $file

samtools faidx $file

module purge
module load gcc
wget https://cran.rstudio.com/src/base/R-4/R-4.0.2.tar.gz

tar -xvf R-4.0.2.tar.gz
cd R-4.0.2/

ml readline
ml bzip2
ml lzma
ml xz
ml pcre2
ml curl
ml jdk


./configure --prefix=/work/LAS/jfw-lab/corrinne/bin/R --with-x=no

make && make install
ml 



