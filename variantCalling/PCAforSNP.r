




# You can either install these once and then comment them out, or you need to add the 'repos' argument
# If you do not specify where your repos is downloaded from, an interactive prompt comes up and will 
# probably stall the script

install.packages("vcfR", repos='http://cran.us.r-project.org')
install.packages("ggrepel", repos='http://cran.us.r-project.org')
install.packages("gdsfmt", repos='http://cran.us.r-project.org')
install.packages("BiocManager", repos='http://cran.us.r-project.org')
install.packages("stringr", repos='http://cran.us.r-project.org')
install.packages("ggplot2", repos='http://cran.us.r-project.org')
install.packages("dplyr", librepos='http://cran.us.r-project.org')
install.packages("devtools", repos='http://cran.us.r-project.org')
BiocManager::install("SNPRelate")


library(devtools)
library(gdsfmt)
library(SNPRelate)
library(vcfR)
library(stringr)
library(ggplot2)
library(ggrepel)
library(dplyr)



args = commandArgs(trailingOnly=TRUE)

#load in VCF file
vcf.fn <- "all.samples.F1.SNPs.noF1.4PCA.recode.vcf"

#Reformat
snpgdsVCF2GDS(vcf.fn, "vcf.gds", method="biallelic.only")
snpgdsSummary("vcf.gds")

#Open the GDS file
genofile <- snpgdsOpen("vcf.gds")

#Get create vectors from genofile
samplenames <- read.gdsn(index.gdsn(genofile, path = "sample.id" ))


############# Turning the sample.id list into a population code list. 
############# Probably an easier way but this works for now
pop_code <- read.csv("species.table", header=T,sep="\t")


#LD-based SNP pruning
set.seed(1000)

# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, autosome.only=FALSE)

snpset.id <- unlist(snpset)

#Principal Component Analysis(PCA)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=12, autosome.only=FALSE)
pc.percent <- pca$varprop*100
# head(round(pc.percent,2)) 
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

PCAtab <- merge(tab,pop_code, by="sample.id", all.x=T)

save.image(file = "all.samples.F1.10pct.gt10.PCA.RData")

ggplot(PCAtab, aes(x=EV2, y=EV1, colour=species, shape=origin, label=sample.id))+ geom_point() +geom_text_repel(aes(label=sample.id),size=2)+scale_colour_manual(values=c("blue","darkgreen"))  
ggsave("all.samples.F1.10pct.gt10.PCA.jpg")

