#!/bin/bash

#SBATCH --job-name=array
#SBATCH --array=1-87%25
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=5
#SBATCH --mem=30G 
#SBATCH --time=4-02:30:00
#SBATCH --output=job.REsample.%A_%a.out 
#SBATCH --error=job.REsample.%A_%a.err
#SBATCH -J "REsample" 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu 

# prior to running the array, run this
# ml git gcc
# git clone --recursive https://github.com/IGBB/reservoir_sample.git
# gcc -lz -o reservoir_sample -Ireadfq reservoir_sample.c
# mv reservoir_sample ../

DIR=/work/LAS/jfw-lab/corrinne/Adom/reads
thr=5
name=$(cat samples.to.keep | sed -n ${SLURM_ARRAY_TASK_ID}p)
fastq=$name.1.fq.gz
trimDir=/work/LAS/jfw-lab/corrinne/Adom/trimmed_reads/
species=$(echo $name | cut -f1 -d '_')
cutLen=90

A1reads=185222
A2reads=190000
Freads=145667

# GSE
# 	Mbp	1% GSE	0.1% GSE	# 85 bp reads	# 90 bp reads	# 95 bp reads
# A1	1667	16.67	1.667	196,118	185,222	175,474
# A2	1710	17.1	1.71	201,176	190,000	180,000
# F	1311	13.11	1.311	154,235	145,667	138,000



if [ "$species" = "A1" ]
    then count=$A1reads
elif [[ "$species" == "A2" || "$species" == "SRR44570"* ]]
    then count=$A2reads
elif [ "$species" = "SRR617704" ] # F genome
    then count =$Freads
else count=0
fi


module load trimmomatic/0.36-lkktrba
trimmomatic SE -threads $thr $DIR/$fastq $name.trimmed.fq ILLUMINACLIP:Adapters.fa:2:30:15 LEADING:28 TRAILING:28 SLIDINGWINDOW:8:28 SLIDINGWINDOW:1:10 MINLEN:92 TOPHRED33
gzip $name.trimmed.fq

./reservoir_sample -l $cutLen -n $name -s $count /work/LAS/jfw-lab/corrinne/Adom/trimmed_reads/$name.trimmed.fq.gz > $name.sampled.fa

mkdir -p /work/LAS/jfw-lab/corrinne/Adom/trimmed_reads
mv $name.trimmed.fq.gz /work/LAS/jfw-lab/corrinne/Adom/trimmed_reads

