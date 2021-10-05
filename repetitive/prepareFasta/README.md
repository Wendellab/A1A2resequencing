# A1A2resequencing
Repetitive analyses for Dual domestication, diversity, and differential introgression in Old World cotton diploids 

Adapters.fa: Adapters to screen for when quality trimming
RepeatExplorerPrep.slurm: gathers files, runs trimmomatic, and runs reservoir_sample
reservoir_sample: script to sample each fastq for a specified number of reads
samples.to.keep: specifies a list of files for RepeatExplorerPrep.slurm and clustering