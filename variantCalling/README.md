# A1A2resequencing
Variant analyses for Dual domestication, diversity, and differential introgression in Old World cotton diploids 

**snpeff**: directory containing snpeff related analyses

**0_prepareFiles.slurm**: download/index the reference genome and install R

**1_sentieon.slurm**: SNP-calling on individual samples

**2_sentieon.GVCFtyper.slurm**: joint genotyping

**3_filter.slurm**: filtering the results of the joint-genotyping 

**keep.samples**: samples passing filters

**PCAforSNP.r**: R script to generate PCA

**remove.F1**: file to remove G. longicalyx (cotton genome designation "F1")

**remove.hybrid**: file to remove putative hybrids

**species.table**: table specifying species information for each sample
