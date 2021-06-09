#!/bin/bash

#SBATCH --job-name=array
#SBATCH --array=1-121%25
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=5
#SBATCH --mem=30G 
#SBATCH --time=4-02:30:00
#SBATCH --output=job.REsample.%A_%a.out 
#SBATCH --error=job.REsample.%A_%a.err
#SBATCH -J "REsample" 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu 

fastq=$(ls reads/* | sed -n ${SLURM_ARRAY_TASK_ID}p)
name=$(basename $fastq | cut -f1 -d '.')
cutLen=90
count=400000

./reservoir_sample -l $cutLen -n $name -s $count $fastq > $name.sampled.2pct.fa


