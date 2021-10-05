if (!requireNamespace('devtools', quietly = TRUE))
       install.packages('devtools')
devtools::install_github('kevinblighe/PCAtools')

install.packages(c("BiocManager","devtools","ggplot2","tidyr","scales","cowplot","ggfortify","corrplot","tidyverse","ggrepel","ggdendro","caret","MASS","reshape2","ROCR"))

library(PCAtools)
library(scales)
library(cowplot)
library(ggfortify)
library(corrplot)
library(tidyverse)
library(ggrepel)
library(ggdendro)
library(caret)
library(MASS)
library(reshape2)
library(ROCR)
#library(gridExtra)
#library(testit)
#library(data.table)
#library(factoextra)

sessionInfo()
setwd('Y:/corrinne/Adom/RepeatExplorer')

# import file with all counts
# evaluate if 0.01% cutoff is reasonable (here, cluster 274)

# Columns added in order to rename for code development 
data <- read.table("all.counts", row.names =1, header = T, sep="\t")
data$annotation <- data$annotation %>% replace_na("none")

contam <- data[(data$annotation=="organelle/plastid" | data$annotation=="rDNA" | data$annotation=="organelle/mitochondria" |  data$annotation=="contamination"),]
contam$cluster <- NULL
contam$F1_lon_ <- NULL
total <- as.data.frame(colSums(contam[,-1]))
total$GS <- c(rep(185222,19),rep(190000,99))
names(total)[1] <- "reads"
total$pct <- total$reads*100/total$GS
#ltSix <- row.names(total[total$pct<6,])
write.table(contam,file="possible.contamination.by.cluster.tbl",quote=F, sep="\t", row.names=T, col.names=T)
write.table(total,file="possible.contamination.tbl",quote=F, sep="\t", row.names=T, col.names=T)

organ <- data[(data$annotation=="organelle/plastid" | data$annotation=="organelle/mitochondria" ),]
organ$cluster <- NULL
organ$F1_lon_ <- NULL
organSum <- aggregate(.~annotation, data=organ, FUN=sum)

totalO <- as.data.frame(colSums(organ[,-1]))
totalO$GS <- c(rep(185222,19),rep(190000,99))
names(totalO)[1] <- "reads"
totalO$pct <- totalO$reads*100/totalO$GS


# Add size/percent columns
data$size <- as.numeric(rowSums(data[,-(1:2)]))
data$percent <- cumsum(data$size)/sum(data$size)
which(data$size>sum(data$size)*0.0001)

# 229 clusters contain >= 0.01% of reads
png("cotton.cutoff.png", 5000, 5000, pointsize=12, res=600)
ggplot(data, aes(x=cluster, y=percent)) + geom_line(size=1) + geom_vline(xintercept=229, color='yellow3', size=1) + scale_x_log10(labels=comma) + scale_y_log10() + geom_vline(xintercept=0, color="grey")+ geom_hline(yintercept=0, color="grey")
dev.off()



##########################################################
################### PCA of scaled data ###################
##########################################################

annotClust <- data[c(1:500),c(2:(ncol(data)-3))]

# find rows where the rowMeans are >10 for either sample
row.names(annotClust[(rowMeans(annotClust[,c(2:20)])>10 | rowMeans(annotClust[,c(21:119)])>10),])
# returns most clusters up to CL0240

clusters <- data[c(1:240),c(2:(ncol(data)-3))]
names(clusters) <- gsub("_$","",names(clusters))


# make table for PCA
counts <- as.data.frame(t(clusters [,(2:(ncol(clusters)))]))
counts$species <- as.factor(gsub("_.*","",row.names(counts)))

# just a few stats, in case we want to see them
table(counts$species)
summary(counts)
ncol(counts)


##################### PCA with prcomp
pr.counts <- prcomp(counts[,-241], scale=TRUE)
summary(pr.counts)

### how many PC capture the information

pveC <- 100*pr.counts$sdev^2/sum(pr.counts$sdev^2)
par(mfrow=c(1,2))
plot(pveC, type ="o", ylab="PVE", xlab="Principal Component", col ="blue")
plot(cumsum(pveC), type="o", ylab ="Cumulative PVE", xlab="Principal Component", col ="brown3")

countpca <- data.frame('species' = counts$species, pr.counts$x[,1:8]) # first eight components
ggplot(data = countpca, aes(x = PC1, y = PC2, col = species)) + 
    geom_point() + scale_color_manual(values=c("A1"="green3", "A2"="darkslateblue")) 


countpca4fig <- countpca[,1:3]
countpca4fig <- countpca4fig %>% mutate(cluster=case_when(PC1 > 15 | PC1 < -20 | PC2 > 17 | PC2 < -15 ~ row.names(countpca4fig), FALSE ~ "")) 
countpca4fig[is.na(countpca4fig)] <- ""

countpcaFig <- ggplot(data = countpca4fig, aes(x = PC1, y = PC2, col = species)) + 
    geom_point() + scale_color_manual(values=c("A1"="green3", "A2"="darkslateblue")) +
        geom_text_repel(box.padding = 0.3,aes(label=cluster)) +
    labs(title = "", x="PC1 (28.1%)", y="PC2 (18.6%)")  



ggsave(filename="Figure_nicePCA.jpg",plot=countpcaFig,device="jpeg",width=8,height=8,units="in",dpi=400)

countpca <- countpca %>% mutate(cluster12=case_when(PC1 > 15 | PC1 < -20 | PC2 > 17 | PC2 < -15 ~ row.names(countpca), FALSE ~ "")) 
countpca <- countpca %>% mutate(cluster23=case_when(PC2 > 17 | PC2 < -15 | PC3 > 15 | PC3 < -15 ~ row.names(countpca), FALSE ~ "")) 
countpca <- countpca %>% mutate(cluster34=case_when(PC3 > 10 | PC3 < -15 | PC4 > 15 | PC4 < -15 ~ row.names(countpca), FALSE ~ "")) 
countpca <- countpca %>% mutate(cluster45=case_when(PC4 > 10 | PC4 < -15 | PC5 > 7 | PC5 < -15 ~ row.names(countpca), FALSE ~ "")) 
countpca <- countpca %>% mutate(cluster13=case_when(PC1 > 15 | PC1 < -20 | PC3 > 15 | PC3 < -15 ~ row.names(countpca), FALSE ~ "")) 
countpca[is.na(countpca)] <- ""

PC1v2 <- ggplot(data = countpca, aes(x = PC1, y = PC2, col = species)) + 
    geom_point() + scale_color_manual(values=c("A1"="green3", "A2"="darkslateblue"), labels=c("G.herbaceum","G.arboreum")) + 
    labs(title = "", x="PC1 (28.1%)", y="PC2 (18.6%)") +
    theme(legend.text=element_text(face="italic"),legend.title=element_blank(),
        legend.position="none",legend.background=element_blank())+
        geom_text_repel(box.padding = 0.3,aes(label=cluster12)) 
  

PC2v3 <- ggplot(data = countpca, aes(x = PC2, y = PC3, col = species)) + 
    geom_point() + scale_color_manual(values=c("A1"="green3", "A2"="darkslateblue")) + 
    labs(title = "", x="PC2 (18.6%)", y="PC3 (8.0%)") + 
    theme(legend.position = "none")+
        geom_text_repel(box.padding = 0.3,aes(label=cluster23)) 


PC3v4 <- ggplot(data = countpca, aes(x = PC3, y = PC4, col = species)) + 
    geom_point() + scale_color_manual(values=c("A1"="green3", "A2"="darkslateblue")) + 
    labs(title = "", x="PC3 (8.0%)", y="PC4 (5.4%)") + 
    theme(legend.position = "none")+
        geom_text_repel(box.padding = 0.3,aes(label=cluster34)) 


PC4v5 <- ggplot(data = countpca, aes(x = PC4, y = PC5, col = species)) + 
    geom_point() + scale_color_manual(values=c("A1"="green3", "A2"="darkslateblue")) + 
    labs(title = "", x="PC4 (5.4%)", y="PC5 (4.6%)") + 
    theme(legend.position = "none")+
        geom_text_repel(box.padding = 0.3,aes(label=cluster45)) 


PC1v3 <- ggplot(data = countpca, aes(x = PC1, y = PC3, col = species)) + 
    geom_point() + scale_color_manual(values=c("A1"="green3", "A2"="darkslateblue"), labels=c("G.herbaceum","G.arboreum")) + 
    labs(title = "", x="PC1 (28.1%)", y="PC3 (8.0%)") +
    theme(legend.text=element_text(face="italic"),legend.title=element_blank(),
        legend.position=c(0.13,0.94),legend.background=element_blank())+
        geom_text_repel(box.padding = 0.3,aes(label=cluster12)) 




pcaGrid <- plot_grid(PC1v2, PC2v3, PC3v4, PC4v5)
ggsave(filename="Figure_PCAgrid.jpg",plot=pcaGrid,device="jpeg",width=8,height=8,units="in",dpi=400)

pcaRow <- plot_grid(PC1v2, PC2v3, PC1v3)


##### remove A2_004, A2_101, A2_076, A2_045
removeThese <- c("A2_004", "A2_101", "A2_045")
filterCounts <- counts[!row.names(counts) %in% removeThese, ]
pr.filter <- prcomp(filterCounts[,-241], scale=TRUE)
summary(pr.filter)


 
filterpca <- data.frame('species' = filterCounts$species, pr.filter$x[,1:8]) # first eight components

filterpca <- filterpca %>% mutate(cluster12=case_when(PC1 > 15 | PC1 < -20 | PC2 > 17 | PC2 < -15 ~ row.names(filterpca), FALSE ~ "")) 
filterpca <- filterpca %>% mutate(cluster23=case_when(PC2 > 17 | PC2 < -15 | PC3 > 15 | PC3 < -15 ~ row.names(filterpca), FALSE ~ "")) 
filterpca <- filterpca %>% mutate(cluster34=case_when(PC3 > 10 | PC3 < -15 | PC4 > 15 | PC4 < -15 ~ row.names(filterpca), FALSE ~ "")) 
filterpca <- filterpca %>% mutate(cluster45=case_when(PC4 > 10 | PC4 < -15 | PC5 > 7 | PC5 < -15 ~ row.names(filterpca), FALSE ~ "")) 
filterpca$cluster <- ""
filterpca["A1_073","cluster"] <- "A1_073"
filterpca["A1_155","cluster"] <- "A1_155"
filterpca[is.na(filterpca)] <- ""

PC1v2f <- ggplot(data = filterpca, aes(x = PC1, y = PC2, col = species)) + 
    geom_point() + scale_color_manual(values=c("A1"="green3", "A2"="darkslateblue"), labels=c("G.herbaceum","G.arboreum")) + 
    labs(title = "", x="PC1 (24.1%)", y="PC2 (18.7%)") +
    theme(legend.text=element_text(face="italic"),legend.title=element_blank(),
        legend.position="none",legend.background=element_blank())+
        geom_text_repel(box.padding = 0.3,aes(label=cluster)) 
  

PC2v3f <- ggplot(data = filterpca, aes(x = PC2, y = PC3, col = species)) + 
    geom_point() + scale_color_manual(values=c("A1"="green3", "A2"="darkslateblue")) + 
    labs(title = "", x="PC2 (18.7%)", y="PC3 (6.9%)") + 
    theme(legend.position = "none")+
        geom_text_repel(box.padding = 0.3,aes(label=cluster)) 


PC3v4f <- ggplot(data = filterpca, aes(x = PC3, y = PC4, col = species)) + 
    geom_point() + scale_color_manual(values=c("A1"="green3", "A2"="darkslateblue")) + 
    labs(title = "", x="PC3 (6.9%)", y="PC4 (4.8%)") + 
    theme(legend.position = "none")+
        geom_text_repel(box.padding = 0.3,aes(label=cluster)) 


PC4v5f <- ggplot(data = filterpca, aes(x = PC4, y = PC5, col = species)) + 
    geom_point() + scale_color_manual(values=c("A1"="green3", "A2"="darkslateblue")) + 
    labs(title = "", x="PC4 (4.8%)", y="PC5 (3.2%)") + 
    theme(legend.position = "none")+
        geom_text_repel(box.padding = 0.3,aes(label=cluster)) 



pcaGridf <- plot_grid(PC1v2f, PC2v3f, PC3v4f, PC4v5f)
ggsave(filename="Figure_PCAgrid_filter.jpg",plot=pcaGridf,device="jpeg",width=8,height=8,units="in",dpi=400)

pcaGridNoFilFil <- plot_grid(PC1v2, PC2v3, PC1v2f, PC2v3f)

ggsave(filename="Figure3_PCAgrid_andfilter.jpg",plot=pcaGridNoFilFil,device="jpeg",width=8,height=8,units="in",dpi=400)

##################### PCA with PCAtools, for visualization
options(ggrepel.max.overlaps = Inf)

### counts

clusMeta <- data.frame(row.names = names(clusters[,2:119]))
clusMeta$species <- gsub("_.*","",row.names(clusMeta))

annotMat <- clusters[,2:119]

p <- pca(annotMat, metadata=clusMeta , scale=T, removeVar=0.1)

byplot <- biplot(p, 
    showLoadings = T, 
    lab = NULL,
    colby='species',
    colkey=c('A1'='green3','A2'='darkslateblue'),
    ellipse = TRUE,
    ellipseConf = 0.95,
    ellipseFill = TRUE,
    ellipseAlpha = 1/4,
    ellipseLineSize = 0,
    pointSize=1,
    legendPosition = 'none', legendLabSize = 12, legendIconSize = 3.0)
ggsave(filename="biplot.jpg",plot=byplot,device="jpeg",width=9,height=6,units="in",dpi=400)

elbow <- findElbowPoint(p$variance)
screeplot(p, components = getComponents(p, 1:20), vline = elbow) + geom_label(aes(x = elbow + 1, y = 50, label = 'Elbow method', vjust = -1, size = 8))

prplot <- pairsplot(p,
    components = getComponents(p, c(1:5)),
    triangle = TRUE, trianglelabSize = 12,
    hline = 0, vline = 0,
    pointSize = 0.4,
    gridlines.major = FALSE, gridlines.minor = FALSE,
    colby = 'species',
    colkey=c('A1'='green3','A2'='darkslateblue'),
    title = '', plotaxes = FALSE,
    margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))
ggsave(filename="pairsPlot.jpg",plot=prplot,device="jpeg",width=10,height=10,units="in",dpi=400)

plplot <- plotloadings(p,
    rangeRetain = 0.01,
    labSize = 4.0,
    labhjust = 1,
    # title = 'Loadings plot',
    # subtitle = 'PC1, PC2, PC3, PC4, PC5',
    caption = 'Top 1% variables',
    shape = 24, shapeSizeRange=c(4,4),
    col = c('lightblue', 'blue', 'darkblue'),
    drawConnectors = TRUE)
# -- variables retained:
# CL000016, CL000076, CL000015, CL000046, CL000057, CL000066, CL000130, CL000233, CL000031, CL000047, CL000053, CL000050, CL000162, CL000149, CL000207
ggsave(filename="loadingsPlot.jpg",plot=plplot,device="jpeg",width=10,height=15,units="in",dpi=400)



##### do samples cluster by species?
sd.data=scale(counts[,-241])

data.dist=dist(sd.data)
hc.out =hclust(dist(sd.data))

dend <- as.dendrogram(hc.out)
dendata <- dendro_data(dend, type = "rectangle")

denplot <- ggplot(dendata$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dendata$labels, aes(x, y, label = label),
            hjust = 0, size = 3)+
  ylim(-3, 80) + coord_flip() + scale_y_reverse()

ggsave(filename="dendroClusterPlot.jpg",plot=denplot,device="jpeg",width=15,height=15,units="in",dpi=400)

weird <- c("A2_045", "A1_155", "A1_073", "A2_004", "A2_101", "A1_133","A1_132")
countpcaHL <- countpca %>% mutate(highlight=case_when(row.names(countpca) %in% weird ~ gsub("_$","",row.names(countpca)), FALSE ~ ""))

weirdplot <- ggplot(data = countpcaHL, aes(x = PC1, y = PC2, col = species)) + 
    geom_point() + scale_color_manual(values=c("A1"="green3", "A2"="darkslateblue")) + 
    geom_text_repel(aes(label=highlight)) + scale_x_reverse()

ggsave(filename="weirdSamplePlot.jpg",plot=weirdplot,device="jpeg",width=10,height=10,units="in",dpi=400)

A2out <- c("A2_087","A2_085","A2_100","A2_123","A2_124","A2_096","A2_065","A2_064","A2_141","A2_119","A2_041","A2_066","A2_076")
countpcaA2 <- countpca %>% mutate(highlight=case_when(row.names(countpca) %in% A2out ~ gsub("_$","",row.names(countpca)), FALSE ~ ""))

A2plot <- ggplot(data = countpcaA2, aes(x = PC1, y = PC2, col = species)) + 
    geom_point() + scale_color_manual(values=c("A1"="green3", "A2"="darkslateblue")) + 
    geom_text_repel(box.padding = 1.5,aes(label=highlight),min.segment.length = 0,xlim=c(-10,-40)) #+ scale_x_reverse()

ggsave(filename="A2SamplePlot.jpg",plot=A2plot,device="jpeg",width=10,height=10,units="in",dpi=400)


########### are cluster abundances correlated
correlations <- cor(counts[,1:100])
res1 <- cor.mtest(counts[,1:100], conf.level = .99)
png("correlation.plot.png", 10000, 10000, pointsize=12, res=600)
corrplot(correlations, p.mat = res1$p, insig = "blank",type="upper",tl.pos = "td", tl.cex = 0.5,diag=F)
dev.off()


########### build a logistic regression model
### spoiler alert, fails
### probably due to high collinearity and low sampling

smpsize <- floor(0.75 * nrow(counts))
train <- sample(nrow(counts), size = smpsize)
trainct <- as.data.frame(counts[train, ])
testct <- as.data.frame(counts[-train, ])

countglm <- glm(species ~.,family=binomial(link='logit'),data=trainct)
summary(countglm)

# way too many variables are correlated
colin <- attributes(alias(countglm)$Complete)$dimnames[[1]]
#  [1] "CL000088" "CL000089" "CL000090" "CL000091" "CL000092" "CL000093" "CL000094" "CL000095" "CL000096"
# [10] "CL000097" "CL000098" "CL000099" "CL000100" "CL000101" "CL000102" "CL000103" "CL000104" "CL000105"
# [19] "CL000106" "CL000107" "CL000108" "CL000109" "CL000110" "CL000111" "CL000112" "CL000113" "CL000114"
# [28] "CL000115" "CL000116" "CL000117" "CL000118" "CL000119" "CL000120" "CL000121" "CL000122" "CL000123"
# [37] "CL000124" "CL000125" "CL000126" "CL000127" "CL000128" "CL000129" "CL000130" "CL000131" "CL000132"
# [46] "CL000133" "CL000134" "CL000135" "CL000136" "CL000137" "CL000138" "CL000139" "CL000140" "CL000141"
# [55] "CL000142" "CL000143" "CL000144" "CL000145" "CL000146" "CL000147" "CL000148" "CL000149" "CL000150"
# [64] "CL000151" "CL000152" "CL000153" "CL000154" "CL000155" "CL000156" "CL000157" "CL000158" "CL000159"
# [73] "CL000160" "CL000161" "CL000162" "CL000163" "CL000164" "CL000165" "CL000166" "CL000167" "CL000168"
# [82] "CL000169" "CL000170" "CL000171" "CL000172" "CL000173" "CL000174" "CL000175" "CL000176" "CL000177"
# [91] "CL000178" "CL000179" "CL000180" "CL000181" "CL000182" "CL000183" "CL000184" "CL000185" "CL000186"
#[100] "CL000187" "CL000188" "CL000189" "CL000190" "CL000191" "CL000192" "CL000193" "CL000194" "CL000195"
#[109] "CL000196" "CL000197" "CL000198" "CL000199" "CL000200" "CL000201" "CL000202" "CL000203" "CL000204"
#[118] "CL000205" "CL000206" "CL000207" "CL000208" "CL000209" "CL000210" "CL000211" "CL000212" "CL000213"
#[127] "CL000214" "CL000215" "CL000216" "CL000217" "CL000218" "CL000219" "CL000220" "CL000221" "CL000222"
#[136] "CL000223" "CL000224" "CL000225" "CL000226" "CL000227" "CL000228" "CL000229" "CL000230" "CL000231"
#[145] "CL000232" "CL000233" "CL000234" "CL000235" "CL000236" "CL000237" "CL000238" "CL000239" "CL000240"

redcounts <- counts[,!(names(counts) %in% colin)]

redsmpsize <- floor(0.75 * nrow(redcounts))
redtrain <- sample(nrow(redcounts), size = redsmpsize)
redtrainct <- as.data.frame(redcounts[redtrain, ])
redtestct <- as.data.frame(redcounts[-redtrain, ])

redcountglm <- glm(species ~.,family=binomial(link='logit'),data=redtrainct)
summary(redcountglm)


########### samples are small and variables are correlated
######## try a linear discriminant analysis on the PCA reduced data
##### elbow suggests that 9 PC capture most of the variance

ldapca <- countpca
ldapca$species <- gsub("A","",ldapca$species)


### make train/test dataset

smpsize <- floor(0.75 * nrow(ldapca))
trainmat <- sample(nrow(ldapca), size = smpsize)
traindf <- as.data.frame(ldapca[trainmat, ])
testdf <- as.data.frame(ldapca[-trainmat, ])

pca.lda <- lda(species ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8, data = traindf)
pca.lda.predict <- predict(pca.lda, newdata = testdf)


### ROC-AUC plot

# Get the posteriors as a dataframe.
ldapost <- as.data.frame(pca.lda.predict$posterior)

# Evaluate the model for ability to predict A1 and plot
# spoiler alert, it's terrible
predA1 <- prediction(ldapost[,1], testdf$species)
rocA1 = performance(predA1, measure = "tpr", x.measure = "fpr")
aucA1 <- performance(predA1, measure = "auc")
aucA1 <- aucA1@y.values

plot(rocA1)
abline(a=0, b= 1)
text(x = .25, y = .65 ,paste("AUC = ", round(aucA1[[1]],3), sep = ""))


# Evaluate the model for ability to predict A2 and plot
# I think it is predicting everything as A2 because of sampling bias
predA2 <- prediction(ldapost[,2], testdf$species)
rocA2 = performance(predA2, measure = "tpr", x.measure = "fpr")
aucA2 <- performance(predA2, measure = "auc")
aucA2 <- aucA2@y.values

plot(rocA2)
abline(a=0, b= 1)
text(x = .25, y = .65 ,paste("AUC = ", round(aucA2[[1]],3), sep = ""))

##### so, I don't think this will work
##### if we get wild, we can add over/under sampling
##### this will probably just compound collinearity problem


#################################################################
#################################################################
#################################################################
#########what if i used DESeq????################################
#################################################################
#################################################################
#################################################################

DEcountdata <- data[rowSums(data[,c(3:(ncol(data)-3))])>5,c(3:(ncol(data)-3))]
# nrow(DEcountdata)
# [1] 42560


#################################################################
#################################################################
#################################################################
############ let's do that later ################################
#################################################################
#################################################################
#################################################################

########### characterize composition ###########

# 9 multiplier represents # reads (x) * 90nt/read * 1 kb/1000nt * 100% = # reads * 9 = # Kb in entire genome for that cluster 
Kb <- data.frame(clusters[1], apply(clusters[2:ncol(clusters)], 2, function (x) x*9))
Kbsum <- aggregate(. ~annotation, data=Kb, FUN=sum)

Kbsum$A1 <- rowMeans(Kbsum[,2:20])
Kbsum$A2 <- rowMeans(Kbsum[,21:119])

Kbsum$A1min <- apply(Kbsum[,2:20], 1, min)
Kbsum$A2min <- apply(Kbsum[,21:119], 1, min)

Kbsum$A1max <- apply(Kbsum[,2:20], 1, max)
Kbsum$A2max <- apply(Kbsum[,21:119], 1, max)

Kbsum <- Kbsum[,-(2:119)]
Kbm <- melt(Kbsum[,-(4:7)])
min <- c(Kbsum$A1min, Kbsum$A2min)
max <- c(Kbsum$A1max, Kbsum$A2max)
Kbm$min <- min
Kbm$max <- max

limits <- aes(ymax=Kbm$max, ymin=Kbm$min)
dodge <- position_dodge(width=0.9)

TEamount <- ggplot(Kbm, aes(x=annotation, y=value, fill = variable)) + 
    geom_bar(stat = "identity",position = dodge) + 
    scale_y_continuous(labels=comma) + 
    scale_fill_manual(breaks=c("A1", "A2"), values=c("green3", "darkslateblue"), labels=c("G.herbaceum","G.arboreum")) + 
    geom_errorbar(limits, position = dodge) + 
    labs(title = "A", x="Broad element category", y="Aggregate amount (mean) in kilobases") + 
    theme_set(theme_grey(base_size=12)) + 
    theme(axis.text = element_text(size = rel(1.5)), 
        plot.title=element_text(face="bold", hjust=0.5),
        axis.title.x = element_text(face="bold", hjust=0.5, vjust=9), 
        axis.title.y = element_text(face="bold", vjust=2),
        axis.text.x = element_text(size = 12, angle = 45, hjust=0.7),
        axis.text.y = element_text(size = 12, angle = 45, vjust=0.1)) + 
    theme(legend.text=element_text(face="italic"),legend.title=element_blank(),
        legend.position=c(0.9,0.94),legend.background=element_blank())

    

ggsave(filename="TEamounts.jpg",plot=TEamount,device="jpeg",width=10,height=10,units="in",dpi=400)


sum(Kbsum$A1)/1000 # convert to Mb
# [1] 1138.838
sum(Kbsum$A2)/1000 # convert to Mb
# [1] 1196.713


########### compare A1 vs A2 ###########

### plot a 1:1 line and display cluster comparisons relative to 1:1 expection
# remember, these have nearly equivalent genome sizes

TEkeep <- c("LTR/Gypsy","LTR","LTR/Copia","DNA/MULE-MuDR","DNA/TcMar-Mariner")

# mean
A1A2 <- clusters[,c(1:119)]
A1A2<- A1A2[A1A2$annotation %in% TEkeep,]
A1A2$A1 <- rowMeans(A1A2[,2:20])
A1A2$A2 <- rowMeans(A1A2[,21:119])
A1A2$signK <- {ifelse((A1A2$A1 > A1A2$A2), "positive", "negative")}
A1A2$perc <- A1A2$A1/A1A2$A2

A1A2 <- A1A2 %>% mutate(highlight=case_when(perc > 1.25 | perc < 0.75 ~ annotation, FALSE ~ ""))
A1A2 <- A1A2 %>% mutate(cluster=case_when(perc > 1.25 | perc < 0.75 ~ row.names(A1A2), FALSE ~ ""))
A1A2[is.na(A1A2)] <- ""
A1A2$label <- paste0(A1A2$highlight,A1A2$cluster)
A1A2$label <- gsub("LTR/","",A1A2$label)
A1A2$label <- gsub("CL000","_CL",A1A2$label)


A1vA2 <- ggplot(A1A2, aes(x=A2, y=A1, shape=signK, color=signK)) + 
    geom_point(size=2) + geom_abline(intercept=0, slope=1) + 
    scale_color_manual(values=c("positive"="green3", "negative"="darkslateblue")) +  
    scale_x_continuous(expand = c(0, 0), limits=c(0,750)) + 
    scale_y_continuous(expand = c(0, 0), limits=c(0,750)) + 
    theme(legend.position="none", axis.title.x = element_text(face="italic", hjust=0.5), 
        axis.title.y = element_text(face="italic", vjust=0.5)) + 
    geom_text_repel(box.padding = 0.5,aes(label=highlight)) + 
    labs(title = "All clusters", x="G.arboreum", y="G.herbaceum")

ggsave(filename="TEcomparisonMean.jpg",plot=A1vA2,device="jpeg",width=10,height=10,units="in",dpi=400)


A1vA2c <- ggplot(A1A2, aes(x=A2, y=A1, shape=signK, color=signK)) + 
    geom_point(size=2) + geom_abline(intercept=0, slope=1) + 
    scale_color_manual(values=c("positive"="green3", "negative"="darkslateblue")) +  
    scale_x_continuous(expand = expansion(mult=c(0.05,0)), limits=c(0,3000)) + 
    scale_y_continuous(expand = expansion(mult=c(0.05,0)), limits=c(0,3000)) + 
    theme(legend.position="none", 
        axis.title.x = element_text(face="italic", hjust=0.53, size=14, vjust=5),
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.3), "cm"), 
        axis.title.y = element_text(face="italic", vjust=-1, size=14, hjust=0.53)) + 
    geom_text_repel(box.padding = 0.5,aes(label=label)) + 
    labs(title = "B", x="G.arboreum", y="G.herbaceum")


ggsave(filename="TEcomparisonMeanLabel.jpg",plot=A1vA2c,device="jpeg",width=10,height=10,units="in",dpi=400)


FigPlot <- plot_grid(TEamount,A1vA2c)

ggsave(filename="Figure_TEtwoPanel.jpg",plot=FigPlot,device="jpeg",width=15,height=9,units="in",dpi=400)


# let's also plot median
A1A2m <- clusters[,c(1:119)]
A1A2m<- A1A2m[A1A2m$annotation %in% TEkeep,]
A1A2m$A1 <- apply(A1A2m[,2:20], 1, median, na.rm = T)
A1A2m$A2 <- apply(A1A2m[,21:119], 1, median, na.rm = T)
A1A2m$signK <- {ifelse((A1A2m$A1 > A1A2m$A2), "positive", "negative")}
A1A2m$perc <- A1A2m$A1/A1A2m$A2

A1A2m <- A1A2m %>% mutate(highlight=case_when((perc > 1.25 | perc < 0.8) & (A1 > 500 | A2 > 500)  ~ annotation, FALSE ~ ""))
A1A2m[is.na(A1A2m)] <- ""

A1A2p <- ggplot(A1A2m, aes(x=A2, y=A1, shape=signK, color=signK)) + 
    geom_point(size=2) + geom_abline(intercept=0, slope=1) + 
    scale_color_manual(values=c("positive"="green3", "negative"="darkslateblue")) +  
    scale_x_continuous(expand = c(0, 0), limits=c(0,3100)) + 
    scale_y_continuous(expand = c(0, 0), limits=c(0,3100)) + 
    theme(legend.position="none", 
        axis.title.x = element_text(face="italic", hjust=0.5,vjust=5), 
        axis.title.y = element_text(face="italic", vjust=0.5,hjust=5)) + 
    geom_text_repel(box.padding = 0.5,aes(label=highlight)) + 
    labs(title = "All clusters", x="G.arboreum", y="G.herbaceum")

ggsave(filename="TEcomparisonMedian.jpg",plot=A1A2p,device="jpeg",width=10,height=10,units="in",dpi=400)



##### t.test between A1 and A2
A1A2 <- clusters[,c(1:119)]
A1A2<- A1A2[A1A2$annotation %in% TEkeep,]

A1A2$p.value <- "calculating..."

for (i in 1:nrow(A1A2))
{
  # Run t.test and then pull $p.value and assign 
  A1A2$p.value[i] <-  t.test(A1A2[i,c(2:20)], A1A2[i,c(21:119)])$p.value
}

A1A2$padj <- p.adjust(A1A2$p.value, method="BH")

nrow(A1A2)
# 175

nrow(A1A2[A1A2$padj< 0.01,])
# 64

row.names(A1A2[A1A2$padj< 0.01,])
# [1] "CL000002" "CL000003" "CL000004" "CL000005" "CL000008" "CL000009"
# [7] "CL000010" "CL000011" "CL000012" "CL000013" "CL000014" "CL000015"
#[13] "CL000016" "CL000017" "CL000020" "CL000021" "CL000023" "CL000024"
#[19] "CL000025" "CL000026" "CL000027" "CL000030" "CL000032" "CL000033"
#[25] "CL000034" "CL000037" "CL000038" "CL000040" "CL000046" "CL000047"
#[31] "CL000050" "CL000052" "CL000053" "CL000057" "CL000058" "CL000065"
#[37] "CL000066" "CL000067" "CL000072" "CL000081" "CL000093" "CL000099"
#[43] "CL000100" "CL000101" "CL000113" "CL000118" "CL000128" "CL000133"
#[49] "CL000137" "CL000146" "CL000148" "CL000149" "CL000152" "CL000155"
#[55] "CL000168" "CL000170" "CL000174" "CL000176" "CL000184" "CL000196"
#[61] "CL000198" "CL000210" "CL000219" "CL000238"

sigClus <- (A1A2[A1A2$padj< 0.01,])
table(sigClus$annotation)
write.table(A1A2,file="SuppTable8.distinguishingclusters.tbl",quote=F, sep="\t", row.names=T, col.names=T)


Allcorrelations <- cor(counts[,1:240])
allres1 <- cor.mtest(counts[,1:240], conf.level = .99)
pmat <- allres1$p

sigMat <- as.data.frame(replace(Allcorrelations, pmat >0.01 , ""))
sigMat[lower.tri(sigMat,diag=TRUE)]=""
write.table(sigMat,file="SuppTable9.maybe.tbl",quote=F, sep="\t", row.names=T, col.names=T)

ggcorrplot(Allcorrelations,type="lower", p.mat=allres1$p, insig="blank",outline.col = "white")

library(corrr)
library(igraph)
library(ggraph)

corMat <- clusters[,c(1:119)]
corMat <- corMat[corMat$annotation %in% TEkeep,-1]
#corMat <- corMat[row.names(corMat) %in% row.names(A1A2[A1A2$padj< 0.01,]),]


tidy_cors <- t(corMat) %>% correlate() %>% stretch()
graph_cors <- tidy_cors %>% 
  filter(abs(r) > 0.7 ) %>% 
  graph_from_data_frame(directed = FALSE)

passed <- row.names(A1A2[A1A2$padj< 0.01,])
nopass <- setdiff(row.names(corMat),passed)

V(graph_cors)$color <- NA
V(graph_cors)$color[V(graph_cors)$name %in% passed] <- "A1 versus A2"
V(graph_cors)$color[V(graph_cors)$name %in% nopass] <- "both As"



gg1 <- ggraph(graph_cors) +
  geom_edge_link(edge_color="grey") +
  geom_node_point(aes(color=color)) +
  geom_node_text(aes(label = name, color=color), repel = TRUE) +
  theme_graph() + 
  scale_color_manual(values=c("A1 versus A2"="black", "both As"="gray50")) 


ggsave(filename="SuppFigure5_corrNet.jpg",plot=gg1,device="jpeg",width=10,height=10,units="in",dpi=400)



# table(A1A2$annotation)
# DNA/MULE-MuDR           LTR     LTR/Copia     LTR/Gypsy 
#             2            15            17           141 

# table(sigClus$annotation)
# LTR LTR/Copia LTR/Gypsy 
#  3         7        54 
