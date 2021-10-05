#!/bin/bash

#SBATCH --job-name=array
#SBATCH --array=0-383
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH -A cotton_genomics
#SBATCH --time=6-22:30:00
#SBATCH -J "REDOhiksarray"
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu
#SBATCH --error=Agenome.hiFilter.pamlRun.%A_%a.err
#SBATCH --output=Agenome.hiFilter.pamlRun.%A_%a.out

# missed files are found with
# for a in out/*.out; do comm -13 <(cut -f1 $a | cut -f2 -d '.' | sed 's/Window//g' | sort -n) iterations > $(basename $a).missed; done

# missed files are incorporated with
# for a in out/*.out; do name=$(cut -f 4 -d '.' <(echo $a | sed 's/260283/256003/g')); cat Agenome.hiFilter.pamlRun.$name.out >> $a; done

module load parallel

echo -e 'Window' ' \t' 'Tree' ' \t' 'F' ' \t' 'F_haplotype' ' \t' 'A1' ' \t' 'A1_haplotype' ' \t' 'A2' ' \t' 'A2_haplotype' ' \t' 'A1_A2_dS'

tmpdir=/dev/shm/$SLURM_JOB_ID.$SLURM_ARRAY_TASK_ID.tmp

mkdir $tmpdir
cp /project/cotton_genomics/Adiversity/Agenome.trees $tmpdir
cp /project/cotton_genomics/Adiversity/paml4.9j/bin/codeml $tmpdir
cp /project/cotton_genomics/Adiversity/pamlRun.py $tmpdir
cp /project/cotton_genomics/Adiversity/hiFilter/hiFilter.fofn $tmpdir
cp /project/cotton_genomics/Adiversity/hiFilter/Agenome.hiFilter.pamlRun.260283_$SLURM_ARRAY_TASK_ID.out.missed $tmpdir


file=$(sed -n $[SLURM_ARRAY_TASK_ID + 1]p /project/cotton_genomics/Adiversity/hiFilter/hiFilter.fofn)
mkdir -p $tmpdir/$(dirname $file)
cp /project/cotton_genomics/Adiversity/hiFilter/$file $tmpdir/$file


parallel -j 48 --timeout 300 --wd $tmpdir "python pamlRun.py hiFilter {} $SLURM_ARRAY_TASK_ID; rm $SLURM_ARRAY_TASK_ID.{}.paml $SLURM_ARRAY_TASK_ID.{}.ctl $SLURM_ARRAY_TASK_ID.{}.out; echo {} >&2" :::: Agenome.hiFilter.pamlRun.260283_$SLURM_ARRAY_TASK_ID.out.missed
rm -r $tmpdir
