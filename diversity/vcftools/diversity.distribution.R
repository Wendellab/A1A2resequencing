library(fBasics)
library(ggplot2)

b <- read.table("both.piDiversity.windowed.pi", header=T)
a1 <- read.table("A1.piDiversity.windowed.pi", header=T)
a2 <- read.table("A2.piDiversity.windowed.pi", header=T)

b <- b[,c(1:2,5)]
a1 <- a1[,c(1:2,5)]
a2 <- a2[,c(1:2,5)]

names(b) <- c("chrom","start","Pi_both")
names(a1) <- c("chrom","start","Pi_A1")
names(a2) <- c("chrom","start","Pi_A2")

ba1 <- merge(b,a1,by=c("chrom","start"))
all <- merge(ba1,a2,by=c("chrom","start"))


allStats <- cbind(
    (basicStats(all[,3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(all[all$chrom=="F01",3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(all[all$chrom=="F02",3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(all[all$chrom=="F03",3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(all[all$chrom=="F04",3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(all[all$chrom=="F05",3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(all[all$chrom=="F06",3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(all[all$chrom=="F07",3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(all[all$chrom=="F08",3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(all[all$chrom=="F09",3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(all[all$chrom=="F10",3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(all[all$chrom=="F11",3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(all[all$chrom=="F12",3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(all[all$chrom=="F13",3:5])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]))

ggA1 <- a1
names(ggA1) <- c("chrom","species","Pi")
ggA1$species <- "A1"

ggA2 <- a2
names(ggA2) <- c("chrom","species","Pi")
ggA2$species <- "A2"

ggDiversity <- rbind(ggA1,ggA2)

ggplot(ggDiversity, aes(x=Pi)) + geom_line(aes(linetype=species,color=species),stat="density",size=1) + scale_color_manual(labels=c("G. herbaceum","G. arboreum"), values=c("green4","blue4"))
ggsave("Figure_diversity.distribution.jpg",width=8,height=5,units="in")


ggplot(all, aes(x=Pi_A1,y=Pi_A2)) + geom_point() + geom_smooth(method=lm,se=T)+ geom_abline(intercept = 0, slope = 1,color="red")+stat_regline_equation(label.y = 0.015)+stat_cor(label.y = 0.010) 
ggsave("Figure_diversity.distribution.stat.jpg",width=8,height=5,units="in")


ggplot(ggDiversity, aes(x=Pi)) + geom_line(aes(linetype=species),stat="density",size=1) + facet_wrap(~chrom)
ggsave("Figure_diversity.distribution.byChrom.jpg",width=8,height=5,units="in")

ggplot(all, aes(x=Pi_A1,y=Pi_A2)) + geom_point() + geom_smooth(method=lm,se=T)+ geom_abline(intercept = 0, slope = 1,color="red")+stat_regline_equation(label.y = 0.015)+stat_cor(label.y = 0.010)+ facet_wrap(~chrom)
ggsave("Figure_diversity.distribution.stat.byChrom.jpg",width=9,height=6,units="in")
 



bExon <- read.table("both.exon.piDiversity.windowed.pi", header=T)
bExon <- bExon[,c(1:2,5)]
names(bExon) <- c("chrom","start","Pi_both")
bExon$type <- "exon"

bIntron <- read.table("both.intron.piDiversity.windowed.pi", header=T)
bIntron <- bIntron[,c(1:2,5)]
names(bIntron) <- c("chrom","start","Pi_both")
bIntron$type <- "intron"

bIntergenic <- read.table("both.intergenic.piDiversity.windowed.pi", header=T)
bIntergenic <- bIntergenic[,c(1:2,5)]
names(bIntergenic) <- c("chrom","start","Pi_both")
bIntergenic$type <- "intergenic"

bothType <- rbind(bExon,bIntron,bIntergenic)

ggplot(bothType, aes(x=Pi_both)) + geom_line(aes(linetype=type,color=type),stat="density",size=1) + scale_x_log10()


bothMelt <- melt(bothType[,c(3:4)])


A1Exon <- read.table("A1.exon.piDiversity.windowed.pi", header=T)
A1Exon <- A1Exon[,c(1:2,5)]
names(A1Exon) <- c("chrom","start","Pi_A1")
A1Exon$type <- "exon"

A1Intron <- read.table("A1.intron.piDiversity.windowed.pi", header=T)
A1Intron <- A1Intron[,c(1:2,5)]
names(A1Intron) <- c("chrom","start","Pi_A1")
A1Intron$type <- "intron"

A1Intergenic <- read.table("A1.intergenic.piDiversity.windowed.pi", header=T)
A1Intergenic <- A1Intergenic[,c(1:2,5)]
names(A1Intergenic) <- c("chrom","start","Pi_A1")
A1Intergenic$type <- "intergenic"

A1Type <- rbind(A1Exon,A1Intron,A1Intergenic)

ggplot(A1Type, aes(x=Pi_A1)) + geom_line(aes(linetype=type,color=type),stat="density",size=1) + scale_x_log10()

A1Melt <- melt(A1Type[,c(3:4)])

A2Exon <- read.table("A2.exon.piDiversity.windowed.pi", header=T)
A2Exon <- A2Exon[,c(1:2,5)]
names(A2Exon) <- c("chrom","start","Pi_A2")
A2Exon$type <- "exon"

A2Intron <- read.table("A2.intron.piDiversity.windowed.pi", header=T)
A2Intron <- A2Intron[,c(1:2,5)]
names(A2Intron) <- c("chrom","start","Pi_A2")
A2Intron$type <- "intron"

A2Intergenic <- read.table("A2.intergenic.piDiversity.windowed.pi", header=T)
A2Intergenic <- A2Intergenic[,c(1:2,5)]
names(A2Intergenic) <- c("chrom","start","Pi_A2")
A2Intergenic$type <- "intergenic"

A2Type <- rbind(A2Exon,A2Intron,A2Intergenic)

ggplot(A2Type, aes(x=Pi_A2)) + geom_line(aes(linetype=type,color=type),stat="density",size=1) + scale_x_log10()


A2Melt <- melt(A2Type[,c(3:4)])

alllllMelt <- rbind(bothMelt, A1Melt, A2Melt)

ggplot(alllllMelt, aes(x=value)) + geom_line(aes(linetype=type,color=variable),stat="density",size=1) + scale_x_log10()
ggsave("PibyType.jpg", width=9, height=6, units="in",dpi=300)

ggplot(alllllMelt,aes(x=type,y=value,fill=variable))+geom_violin()+scale_y_log10()



intronStats <- cbind(
    (basicStats(bothType[bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F01" & bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F02" & bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F03" & bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F04" & bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F05" & bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F06" & bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F07" & bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F08" & bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F09" & bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F10" & bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F11" & bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F12" & bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F13" & bothType$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F01" & A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F02" & A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F03" & A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F04" & A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F05" & A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F06" & A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F07" & A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F08" & A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F09" & A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F10" & A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F11" & A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F12" & A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F13" & A1Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F01" & A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F02" & A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F03" & A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F04" & A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F05" & A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F06" & A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F07" & A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F08" & A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F09" & A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F10" & A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F11" & A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F12" & A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F13" & A2Type$type=="intron",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]))

intronStats <- intronStats[,seq(2,ncol(intronStats),2)]


exonStats <- cbind(
    (basicStats(bothType[bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F01" & bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F02" & bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F03" & bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F04" & bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F05" & bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F06" & bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F07" & bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F08" & bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F09" & bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F10" & bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F11" & bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F12" & bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F13" & bothType$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F01" & A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F02" & A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F03" & A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F04" & A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F05" & A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F06" & A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F07" & A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F08" & A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F09" & A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F10" & A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F11" & A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F12" & A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F13" & A1Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F01" & A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F02" & A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F03" & A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F04" & A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F05" & A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F06" & A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F07" & A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F08" & A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F09" & A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F10" & A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F11" & A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F12" & A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F13" & A2Type$type=="exon",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]))

exonStats <- exonStats[,seq(2,ncol(exonStats),2)]

intergenicStats <- cbind(
    (basicStats(bothType[bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F01" & bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F02" & bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F03" & bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F04" & bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F05" & bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F06" & bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F07" & bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F08" & bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F09" & bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F10" & bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F11" & bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F12" & bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(bothType[bothType$chrom=="F13" & bothType$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F01" & A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F02" & A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F03" & A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F04" & A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F05" & A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F06" & A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F07" & A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F08" & A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F09" & A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F10" & A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F11" & A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F12" & A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A1Type[A1Type$chrom=="F13" & A1Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F01" & A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F02" & A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F03" & A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F04" & A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F05" & A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F06" & A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F07" & A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F08" & A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F09" & A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F10" & A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F11" & A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F12" & A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(A2Type[A2Type$chrom=="F13" & A2Type$type=="intergenic",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]))

intergenicStats <- intergenicStats[,seq(2,ncol(intergenicStats),2)]

write.table(intronStats,file="intronDiversity.tbl", quote=F, sep='\t')
write.table(exonStats,file="exonDiversity.tbl", quote=F, sep='\t')
write.table(intergenicStats,file="intergenicDiversity.tbl", quote=F, sep='\t')

