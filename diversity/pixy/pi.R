library(data.table)
library(tidyverse)
library(purrr)
library(ggplot2)

setDTthreads(10)
options(datatable.fread.datatable=FALSE)

pixyPi <- fread("pixy_pi.txt",header=T, sep="\t")

piA1 <- pixyPi %>% filter(pop=="A1") %>% select(pop,chromosome,avg_pi)
piA2 <- pixyPi %>% filter(pop=="A2") %>% select(pop,chromosome,avg_pi)

piA1Summary <- piA1 %>%
	split(.$chromosome) %>%
	map(summary)

piA2Summary <- piA2 %>%
	split(.$chromosome) %>%
	map(summary)



piFlip <- pixyPi %>%
	select(pop,chromosome,window_pos_1, avg_pi) %>%
	filter(!pop=="F1") %>%
	rename(Species=pop, Chrom=chromosome, Position=window_pos_1, Pi=avg_pi) %>% 
	transform(.,PiChart=ifelse(Species=="A1",Pi,-1*Pi))
	
panelA <- ggplot(piFlip, aes(x=Position,y=PiChart,color=Species)) + geom_line() + facet_grid(Chrom ~ .) + scale_color_manual(labels=c("G. herbaceum","G. arboreum"), values=c("green4","blue4"))

ggsave(panelA, filename="PixyPiDiversity.jpg", width=25, height=15, units="in")

pixyPiExon <- fread("all.exon_pi.txt",header=T, sep="\t")
pixyPiIntron <- fread("all.intron_pi.txt",header=T, sep="\t")
pixyPiIntergenic <- fread("all.intergenic_pi.txt",header=T, sep="\t")


piAE <- pixyPiExon %>% select(pop,chromosome,avg_pi) %>% filter(!pop=="F1") %>% filter(!is.na(avg_pi)) %>% mutate(gene="exon")

piAI <- pixyPiIntron %>% select(pop,chromosome,avg_pi) %>% filter(!pop=="F1") %>% filter(!is.na(avg_pi)) %>% mutate(gene="intron")

piAIG <- pixyPiIntergenic %>% select(pop,chromosome,avg_pi) %>% filter(!pop=="F1") %>% filter(!is.na(avg_pi)) %>% mutate(gene="intergenic")

piAC <- bind_rows(piAE, piAI, piAIG)


ggplot(piAC, aes(x=avg_pi,color=pop,linetype=gene)) +
	geom_density() +
	xlim(0,0.0021) +
	scale_y_log10() +
	scale_color_manual(labels=c("G. herbaceum","G. arboreum"), values=c("darkgreen","blue"))


