library(ggplot2)
library(fBasics)


Fst <- read.table("A1vA2.fst.windowed.weir.fst", header=T)
Fst <- Fst[,c(1,2,5)]
names(Fst) <- c("Chrom","Pos","fst")

ggplot(Fst, aes(x=Pos,y=fst)) + geom_line() + facet_grid(Chrom ~ .) + geom_hline(yintercept=0.5, color="red3")+ geom_hline(yintercept=c(0.25,0.75), linetype="dashed", color="red3")
ggsave(filename="fst_chromosome.jpg", width=20, height=20, units="in")



allStats <- cbind(
    (basicStats(Fst[,2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(Fst[Fst$Chrom=="F01",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(Fst[Fst$Chrom=="F02",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(Fst[Fst$Chrom=="F03",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(Fst[Fst$Chrom=="F04",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(Fst[Fst$Chrom=="F05",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(Fst[Fst$Chrom=="F06",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(Fst[Fst$Chrom=="F07",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(Fst[Fst$Chrom=="F08",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(Fst[Fst$Chrom=="F09",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(Fst[Fst$Chrom=="F10",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(Fst[Fst$Chrom=="F11",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(Fst[Fst$Chrom=="F12",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(Fst[Fst$Chrom=="F13",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]))

allStats <- allStats[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28)]

ggplot(Fst, aes(x=fst)) + geom_line(stat="density",size=1) 
ggsave("Figure_fst.distribution.jpg",width=8,height=5,units="in")

ggplot(Fst, aes(x=fst)) + geom_line(stat="density",size=1) + facet_wrap(~Chrom)
ggsave("Figure_diversity.distribution.byChrom.jpg",width=8,height=5,units="in")

########
lrgFst <- read.table("A1vA2.fst.large.windowed.weir.fst", header=T)

lrgFst <- lrgFst[,c(1,2,5)]
names(lrgFst) <- c("Chrom","Pos","fst")

ggplot(lrgFst, aes(x=Pos,y=fst)) + geom_line(color="red3") + facet_grid(Chrom ~ .) + 
    geom_hline(yintercept=0.5, color="black") + geom_hline(yintercept=c(0.25,0.75), color="grey") + 
    theme(panel.background = element_rect(fill = "grey94"))

ggsave("Figure_Fst.1Mbwindow100kStep.byChr.pretty.jpg",width=15,height=12,units="in")


lrgStats <- cbind(
    (basicStats(lrgFst[,2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgFst[lrgFst$Chrom=="F01",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgFst[lrgFst$Chrom=="F02",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgFst[lrgFst$Chrom=="F03",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgFst[lrgFst$Chrom=="F04",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgFst[lrgFst$Chrom=="F05",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgFst[lrgFst$Chrom=="F06",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgFst[lrgFst$Chrom=="F07",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgFst[lrgFst$Chrom=="F08",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgFst[lrgFst$Chrom=="F09",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgFst[lrgFst$Chrom=="F10",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgFst[lrgFst$Chrom=="F11",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgFst[lrgFst$Chrom=="F12",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgFst[lrgFst$Chrom=="F13",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]))

lrgStats <- lrgStats[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28)]


lrgFst$hline <- ifelse(lrgFst$Chrom=="F01",0.43,
	ifelse(lrgFst$Chrom=="F02",0.43,
	ifelse(lrgFst$Chrom=="F03",0.42,
	ifelse(lrgFst$Chrom=="F04",0.41,
	ifelse(lrgFst$Chrom=="F05",0.48,
	ifelse(lrgFst$Chrom=="F06",0.44,
	ifelse(lrgFst$Chrom=="F07",0.40,
	ifelse(lrgFst$Chrom=="F08",0.41,
	ifelse(lrgFst$Chrom=="F09",0.44,
	ifelse(lrgFst$Chrom=="F10",0.39,
	ifelse(lrgFst$Chrom=="F11",0.44,
	ifelse(lrgFst$Chrom=="F12",0.43,
	ifelse(lrgFst$Chrom=="F13",0.44,0)))))))))))))

ggplot(lrgFst, aes(x=Pos,y=fst)) + geom_line(color="red4") + facet_grid(Chrom ~ .) + 
    geom_hline(aes(yintercept=hline), color="black") + 
    theme(panel.background = element_rect(fill = "grey94"),panel.grid.minor.y = element_blank())


ggsave("Figure1B_Fst.1Mbwindow100kStep.byChr.pretty.jpg",width=15,height=12,units="in")
















########

lrgrFst <- read.table("A1vA2.fst.larger.windowed.weir.fst", header=T)

lrgrFst <- lrgrFst[,c(1,2,5)]
names(lrgrFst) <- c("Chrom","Pos","fst")

ggplot(lrgrFst, aes(x=Pos,y=fst)) + geom_line() + facet_grid(Chrom ~ .) + geom_hline(yintercept=0.5, color="red3")+ geom_hline(yintercept=c(0.25,0.75), linetype="dashed", color="red3")


lrgrStats <- cbind(
    (basicStats(lrgrFst[,2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgrFst[lrgrFst$Chrom=="F01",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgrFst[lrgrFst$Chrom=="F02",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgrFst[lrgrFst$Chrom=="F03",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgrFst[lrgrFst$Chrom=="F04",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgrFst[lrgrFst$Chrom=="F05",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgrFst[lrgrFst$Chrom=="F06",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgrFst[lrgrFst$Chrom=="F07",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgrFst[lrgrFst$Chrom=="F08",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgrFst[lrgrFst$Chrom=="F09",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgrFst[lrgrFst$Chrom=="F10",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgrFst[lrgrFst$Chrom=="F11",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgrFst[lrgrFst$Chrom=="F12",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgrFst[lrgrFst$Chrom=="F13",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]))

lrgrStats <- lrgrStats[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28)]

########

lrgishFst <- read.table("A1vA2.fst.largish.windowed.weir.fst", header=T)

lrgishFst <- lrgishFst[,c(1,2,5)]
names(lrgishFst) <- c("Chrom","Pos","fst")

ggplot(lrgishFst, aes(x=Pos,y=fst)) + geom_line() + facet_grid(Chrom ~ .) + geom_hline(yintercept=0.4430, color="red3")+ geom_hline(yintercept=c(0.2215,0.6645), linetype="dashed", color="red3")
ggsave("Fst.5MbwindowHalfMbStep.byChr.jpg", height=8, width=4, dpi=600)

lrgishStats <- cbind(
    (basicStats(lrgishFst[,2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgishFst[lrgishFst$Chrom=="F01",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgishFst[lrgishFst$Chrom=="F02",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgishFst[lrgishFst$Chrom=="F03",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgishFst[lrgishFst$Chrom=="F04",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgishFst[lrgishFst$Chrom=="F05",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgishFst[lrgishFst$Chrom=="F06",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgishFst[lrgishFst$Chrom=="F07",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgishFst[lrgishFst$Chrom=="F08",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgishFst[lrgishFst$Chrom=="F09",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgishFst[lrgishFst$Chrom=="F10",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgishFst[lrgishFst$Chrom=="F11",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgishFst[lrgishFst$Chrom=="F12",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]),
    (basicStats(lrgishFst[lrgishFst$Chrom=="F13",2:3])[c("Mean", "Stdev", "Median", "Minimum", "Maximum","1. Quartile","3. Quartile"),]))

lrgishStats <- lrgishStats[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28)]
