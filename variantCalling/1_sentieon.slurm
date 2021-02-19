#!/bin/bash

#SBATCH --job-name=array
#SBATCH --array=1-330%20
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=70 
#SBATCH --mem=200G 
#SBATCH --time=7-02:30:00
#SBATCH --output=job.F.%A_%a.out 
#SBATCH --error=job.F.%A_%a.err
#SBATCH -J "F.sentieonArray" 
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu 

module load bwa
module load samtools

DIR=/work/LAS/jfw-lab/corrinne/Adom/reads
ref=longicalyx.fasta
ext=F1
thr=70

file1=$(ls -1 $DIR/*.1.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
file2=$(ls -1 $DIR/*.2.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)

name=$(basename $file1 .1.fq.gz)


module load sentieon-genomics

time ( bwa mem -M -K 10000000 -R "@RG\tID:$name \tSM:$name \tPL:ILLUMINA" -t $thr $ref $file1 $file2 > $name.$ext.sam  ) &> $name.$ext.bwa.timelog
time ( sentieon util sort -o $name.$ext.sort.bam -t $thr --sam2bam -i $name.$ext.sam  ) &> $name.$ext.sam2bam.timelog

#time ( sentieon driver -t $thr -r $ref -i $name.$ext.sort.bam --algo GCBias --summary $name.$ext.GC.summary $name.$ext.GC.metric --algo MeanQualityByCycle $name.$ext.MQ.metric --algo QualDistribution $name.$ext.QD.metric --algo InsertSizeMetricAlgo $name.$ext.IS.metric --algo AlignmentStat $name.$ext.ALN.metric ) &> $name.$ext.metric.timelog
#time ( sentieon plot metrics -o $name.$ext.metric.pdf gc=$name.$ext.GC.metric mq=$name.$ext.MQ.metric qd=$name.$ext.QD.metric isize=$name.$ext.IS.metric ) &> $name.$ext.metricPlot.timelog

time ( sentieon driver -t $thr -i $name.$ext.sort.bam --algo LocusCollector --fun score_info $name.$ext.score ) &> $name.$ext.locusCollect.timelog
time ( sentieon driver -t $thr -i $name.$ext.sort.bam --algo Dedup --rmdup --score_info $name.$ext.score --metrics $name.$ext.dedup.metric $name.$ext.dedup.bam ) &> $name.$ext.rmDup.timelog

time ( sentieon driver -t $thr -r $ref -i $name.$ext.dedup.bam --algo Realigner $name.$ext.realign.bam ) &> $name.$ext.realign.timelog

time ( sentieon driver -t $thr -r $ref -i $name.$ext.realign.bam  --algo Haplotyper $name.$ext.gVCF --emit_mode gvcf ) &> $name.$ext.gvcf.timelog


ml 
