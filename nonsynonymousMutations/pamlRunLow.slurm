#!/bin/bash

#SBATCH --job-name=array
#SBATCH --array=0-561
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH -A cotton_genomics
#SBATCH --time=6-22:30:00
#SBATCH -J "ksarray"
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu
#SBATCH --error=Agenome.lowFilter.pamlRun.%A_%a.err
#SBATCH --output=Agenome.lowFilter.pamlRun.%A_%a.out

module load parallel

echo -e 'Window' ' \t' 'Tree' ' \t' 'F' ' \t' 'F_haplotype' ' \t' 'A1' ' \t' 'A1_haplotype' ' \t' 'A2' ' \t' 'A2_haplotype' ' \t' 'A1_A2_dS'

tmpdir=/dev/shm/$SLURM_JOB_ID.$SLURM_ARRAY_TASK_ID.tmp

mkdir $tmpdir
cp /project/cotton_genomics/Adiversity/Agenome.trees $tmpdir
cp /project/cotton_genomics/Adiversity/paml4.9j/bin/codeml $tmpdir
cp /project/cotton_genomics/Adiversity/pamlRun.py $tmpdir
cp /project/cotton_genomics/Adiversity/lowFilter/lowFilter.fofn $tmpdir

file=$(sed -n $[SLURM_ARRAY_TASK_ID + 1]p /project/cotton_genomics/Adiversity/lowFilter/lowFilter.fofn)
mkdir -p $tmpdir/$(dirname $file)
cp /project/cotton_genomics/Adiversity/lowFilter/$file $tmpdir/$file


parallel -j 48 --timeout 600 --wd $tmpdir "python pamlRun.py lowFilter {} $SLURM_ARRAY_TASK_ID; rm $SLURM_ARRAY_TASK_ID.{}.paml $SLURM_ARRAY_TASK_ID.{}.ctl $SLURM_ARRAY_TASK_ID.{}.out; echo {} >&2" ::: {0..112831}
rm -r $tmpdir
