#!/bin/bash

#SBATCH --job-name=array
#SBATCH --array=114,190,228,422
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH -A cotton_genomics
#SBATCH --time=6-22:30:00
#SBATCH -J "REDOmedksarray"
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu
#SBATCH --error=Agenome.medFilter.pamlRun.%A_%a.err
#SBATCH --output=Agenome.medFilter.pamlRun.%A_%a.out

# missed files are found with
# for a in out/*.out; do comm -13 <(cut -f1 $a | cut -f2 -d '.' | sed 's/Window//g' | sort -n) iterations > $(basename $a).missed; done

# to incorporate new ones
# for a in out/*.out; do name=$(cut -f 4 -d '.' <(echo $a | sed 's/257955/256003/g')); cat Agenome.medFilter.pamlRun.$name.out >> $a; done


module load parallel

echo -e 'Window' ' \t' 'Tree' ' \t' 'F' ' \t' 'F_haplotype' ' \t' 'A1' ' \t' 'A1_haplotype' ' \t' 'A2' ' \t' 'A2_haplotype' ' \t' 'A1_A2_dS'

tmpdir=/dev/shm/$SLURM_JOB_ID.$SLURM_ARRAY_TASK_ID.tmp

mkdir $tmpdir
cp /project/cotton_genomics/Adiversity/Agenome.trees $tmpdir
cp /project/cotton_genomics/Adiversity/paml4.9j/bin/codeml $tmpdir
cp /project/cotton_genomics/Adiversity/pamlRun.py $tmpdir
cp /project/cotton_genomics/Adiversity/mediumFilter/medFilter.fofn $tmpdir
cp /project/cotton_genomics/Adiversity/mediumFilter/Agenome.medFilter.pamlRun.257955_$SLURM_ARRAY_TASK_ID.out.missed $tmpdir


file=$(sed -n $[SLURM_ARRAY_TASK_ID + 1]p /project/cotton_genomics/Adiversity/mediumFilter/medFilter.fofn)
mkdir -p $tmpdir/$(dirname $file)
cp /project/cotton_genomics/Adiversity/mediumFilter/$file $tmpdir/$file


parallel -j 48 --timeout 600 --wd $tmpdir "python pamlRun.py medFilter {} $SLURM_ARRAY_TASK_ID; rm $SLURM_ARRAY_TASK_ID.{}.paml $SLURM_ARRAY_TASK_ID.{}.ctl $SLURM_ARRAY_TASK_ID.{}.out; echo {} >&2" :::: Agenome.medFilter.pamlRun.257955_$SLURM_ARRAY_TASK_ID.out.missed
rm -r $tmpdir
