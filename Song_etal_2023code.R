#Song et al. 2023
#Nature Plants

####################

library("phyloseq")
library ("ggthemes")
library("ggplot2")
library("phyloseq")
library("readr")
library("dplyr")
library("magrittr")
library("DESeq2")
library("BiocParallel")

setwd("/Users/zaydamoralesmoreira/Documents/HaneyLab/Song_etal_NaturePlants_2023")

#Alpha Diversity# 
################

mutants<-read.csv("Reads_psk_sequencing.csv",header=T)
row.names(mutants) <- mutants[,1]       
mutants <- mutants[,-1]

sam<-read.csv("Sample_info_psk_sequencing.csv",header=T)
row.names(sam) <- sam[,1]       
sam <- sam[,-1]

mutants_otu2<-phyloseq(otu_table(mutants, taxa_are_rows = TRUE), sample_data(sam))

p <-plot_richness(mutants_otu2, x="Mutant", measures=c("Observed","Chao1", "Shannon", "InvSimpson"))
p1<-p+theme_bw()+ geom_point(size = 3)+ theme(axis.text.x = element_text(angle = 65))
p2<-p1 + geom_boxplot(data=p1$mutants_otu2, aes(x=Mutant, y=value, fill=Mutant), alpha=0.3)
p3<-p2+scale_fill_brewer(palette="Dark2")
p3
ggsave("Alpha_Diversity_psk.svg", width=8, height=5, units="in", dpi=300)


#Beta Diversity# 
###############

#Hellinger transformation 
mutants<-read.csv("Reads_psk_sequencing_for_Hellinger.csv",header=T)
row.names(mutants) <- mutants[,1]       
mutants <- mutants[,-1]
hell_data <- sqrt(mutants / apply(mutants,1,sum))
write_csv(hell_data, path = "Hellinger_Transformed_reads_psk.csv")

#PCoA
mutants<-read.csv("Transposed_Hellinger_Transformed_reads_psk.csv",header=T)

row.names(mutants) <- mutants[,1]       
mutants <- mutants[,-1]
mutants_otu<-phyloseq(otu_table(mutants, taxa_are_rows = TRUE))

sam<-read.csv("Sample_info_psk_sequencing.csv",header=T)
row.names(sam) <- sam[,1]       
sam <- sam[,-1]

mutants_otu2<-phyloseq(otu_table(mutants, taxa_are_rows = TRUE), sample_data(sam))

ordu = ordinate(mutants_otu2, "PCoA", "bray")
plot_ordination(mutants_otu2, ordu, color="Mutant", axes = c(1, 2))
plot_ordination(mutants_otu2, ordu, color="Mutant", axes = c(1, 3))
plot_ordination(mutants_otu2, ordu, color="Mutant", axes = c(2, 3))

v <- data.frame(ordu$vectors, max.print = 200)

write_csv(v, path = "vectors_PCoA_psk.csv")

#PCoA scale colors
mutants<-read.csv("Coordinates_pco1_pco3_psk.csv",header=T)

mutants$Mutant <- factor(mutants$Mutant ,levels=c("Bulk.Soil", "Col", "hsm7", "pskr13") )
p6 <- ggplot(mutants, aes(x = x, y = y, color=Mutant)) +
  geom_point (size=2, stroke=2) +
  scale_shape_manual(values=c(16, 16, 16, 16,16))+
  ggtitle("") +
  labs(x = "PCO1 (26% of total variation)", y = "PCO3 (6% of total variation)", size=20)


p6 <- p6 +  scale_color_brewer(palette = "Dark2")

p6 + theme_bw()+theme(axis.text=element_text(size=12, color ="black"),
                      axis.title=element_text(size=12)) +theme(legend.text=element_text(size=12))+theme(legend.title =element_text(size=12,face="bold"))


ggsave("PCoA_psk.svg", width=8, height=5, units="in", dpi=300)

#Relative Abundance Family Level#
################################

Tax<- read.csv("Abundance_family_level_psk.csv", header=TRUE)

Tax$Mutant <- factor(Tax$Mutant, levels=unique(as.character(Tax$Mutant)) )

Tax$Family <-factor(Tax$Family, 
                    levels = c( "Others", "Methylophilaceae", "Pirellulaceae", "Acetobacteraceae", "Solirubrobacteraceae","Gemmataceae","Devosiaceae","Rhizobiaceae", "Solimonadaceae","Opitutaceae","Haliangiaceae","Chthoniobacteraceae",  "Microscillaceae", "Solibacteraceae", "Caulobacteraceae", "Pedosphaeraceae","Flavobacteriaceae", "Chitinophagaceae","Pseudomonadaceae","Nitrosomonadaceae","Gemmatimonadaceae","Xanthobacteraceae", "Burkholderiaceae"))

pal <- colorRampPalette(c("gray81", "plum3", "tan3", "darkolivegreen2", "lightgoldenrod2", "cyan", "orange", "deepskyblue","cornsilk", "olivedrab1", "lightcoral" , " royalblue","darkseagreen3", "turquoise", "mediumorchid", "thistle1", "lightpink","lightgoldenrod1","red3","lightblue1","burlywood1","palegreen","slateblue"))


Familybymutant<-ggplot(data=Tax, aes(x=Mutant, y=Abundance, fill=Family)) +
  geom_bar(color="black", stat="identity")+
  scale_fill_manual(values=pal(23)) +
  labs( y = "Relative Abundance(%)")  +theme_gray()+
  
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey80"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 20),
        strip.text = element_text(size = 12),
        axis.title.y = element_text(vjust= 3, size=12),
        axis.title.x = element_text(vjust= 3, size=12),
        panel.border = element_rect(colour="black"),
        axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1)) 
Familybymutant

ggsave("Relative_Abundance_Family_psk.svg", width=7, height=6, units="in", dpi=300)

#Pseudomonas ASVs Abundance# 
###########################

Tax<- read.csv("Pseudomonas_Abundance_ASV_level.csv", header=TRUE)
Tax$Mutant <- factor(Tax$Mutant, levels=unique(as.character(Tax$Mutant)) )

Tax$ASV <-factor(Tax$ASV, 
                 levels = c( "798fa1d63e7f3b60776b936b77cff5b5", "31abef2f0d2bf06a79c49cb18ecdd25b", "f93f9dd6eb5e90b86d7e8d14702a9be6", "3292efc3bafa0e0a73ace3fbc8e74014", "41265245f749afec0a34e50836f18c08","b5401d1c2c282dd45f4e922eca782f49","3df3f4cd23050b7d81be3aa2b44d8824","41b5b61cf5de00b2dddd6e3674c46ba6", "63661ab825ac3031d7e9e7404caf7450","ef7185f833fdb7bd46030c365fa9ee73","3e97b20551faf728055e6517bcbcabbb","3e7caeb1d4364dc0a95468107ca2c41d", "689f30c7d9a4c1ec733a78a5dbb43b76", "8f09bac714c02dc9415770271349dc99", "38bfbd3a2cf83da8fe03cea0f6b4f5f3", "7b5afc21508102a6c81750cd6ebfd8f2","a969769a23c7151eb73cd0e49edd0a10", "65222874c6abba1ac46976534b026449","402e5913597695a16d7cad415ffff02f","b2169de7b5980a96680c0a9cff5fbe7a","c3b9e9a07d3524d9447dfcd38e89a000","bf5f4f4da286a3c774d31edb1b6cdf32", "dc5030eb18fe395a50b8235f2e583dbd"))

pal <- colorRampPalette(c("#B8DE29FF", "#EDD9A3", "#FDE725FF", "#EFCC98", "#F3B584", "#73D055FF", "#F6A97A", "#29AF7FFF","#238A8DFF", "#F79C79", "#2D708EFF" , "#39568CFF","#F89078", "#FA7876", "#F66D7A", "#F2637F", "#EA4F88","#DF488D","#CA3C97","#B1339E","#952EA0","#783B9D","#4B2991"))

library(ggplot2)


ASVbymutant<-ggplot(data=Tax, aes(x=Mutant, y=Abundance, fill=ASV)) +
  geom_bar(color="black", stat="identity")+
  scale_fill_manual(values=pal(23)) +
  labs( y = "Relative abundance(%)")  +theme_gray()+
  theme_bw()  +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey80"),
        plot.title = element_text(size = rel(1.5), 
                                  face = "bold", vjust = 1.5, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 20),
        strip.text = element_text(size = 12),
        axis.title.y = element_text(vjust= 3, size=12),
        axis.title.x = element_text(vjust= 3, size=12), 
        panel.border = element_rect(colour="black"),
        axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 1)) 
ASVbymutant
ggsave("Pseudomonas_Abundance_ASV_psk.svg", width=4, height=6, units="in", dpi=300)

#Deseq Family level# 
###################

countData <- read.csv("Family_Deseq_Revisions.csv", header = TRUE, sep = ",")
head(countData)
metaData <- read.csv("Metadata_Deseq_family.csv", header = TRUE, sep = ",")
metaData
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~Treatment, tidy = TRUE)
dds

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$Treatment <- relevel(dds$Treatment, ref = "Col")

dds <- DESeq(dds)
countsdds <- counts(dds, normalized = TRUE)
write.csv(countsdds, file="dds_norm_countscolpsk.csv")
summary(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res)

res <- results(dds, name="Treatment_psk_vs_Col")
res
resultsNames(dds)

#Double check with contrast
res1 <- results( dds, contrast = c("Treatment", "psk", "Col") )
res1

resLFC <- lfcShrink(dds, coef="Treatment_psk_vs_Col", type="apeglm")
resLFC

register(SnowParam(4))
resOrdered <- res[order(res$pvalue),]
summary(res)

write.csv(as.data.frame(resOrdered), 
          file="Treatment_col_vs_psk_results.csv")

resSig1 <- subset(resOrdered, padj <= 0.05)
resSig1
write.csv(resSig1, file="psk0.05list.csv")

resSig2 <- subset(resOrdered, padj <= 0.1)
resSig2
write.csv(resSig2, file="psk0.1list.csv")

resSig3 <- subset(resSig1, log2FoldChange >=1 | log2FoldChange <=-1)
resSig3
write.csv(resSig3, file="psk0.05fold2list.csv")

psk<-read.csv("psk0.1list.csv")

psk$Family <-factor(psk$Family, 
                    levels = c("Pseudomonadaceae","Caulobacteraceae"))


da_taxa_plot <- ggplot(psk, aes(x = log2FoldChange, y = Family, color = Family,  size = -log10(padj)))+
  geom_point()+
  geom_vline(xintercept = 0, colour="blue", linetype = "dashed")+ 
  theme(axis.title.y=element_blank())+
  theme_bw()

da_taxa_plot
ggsave("Bubble_Graph_Family_Deseq_revisions.svg", width=7, height=6, units="in", dpi=300)                    

#Deseq RNA Seq#
##############

##Col##

countData <- read.csv("RNASeq_Col365_vs_buffer.csv", header = TRUE, sep = ",")
head(countData)
metaData <- read.csv("Metadata_Col365_vs_buffer.csv", header = TRUE, sep = ",")
metaData
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~Treatment, tidy = TRUE)
dds

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$Treatment <- relevel(dds$Treatment, ref = "Buffer")

dds <- DESeq(dds)
countsdds <- counts(dds, normalized = TRUE)
write.csv(countsdds, file="dds_norm_countscol.csv")
summary(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res)

res <- results(dds, name="Treatment_Bacteria_vs_Buffer")
res
resultsNames(dds)

#Double check with contrast
res1 <- results( dds, contrast = c("Treatment", "Bacteria", "Buffer") )
res1

resLFC <- lfcShrink(dds, coef="Treatment_Bacteria_vs_Buffer", type="apeglm")
resLFC

register(SnowParam(4))
resOrdered <- res[order(res$pvalue),]
summary(res)

write.csv(as.data.frame(resOrdered), 
          file="Treatment_Bacteria_vs_Buffer_resultsCol.csv")

resSig1 <- subset(resOrdered, padj <= 0.05)
resSig1
write.csv(resSig1, file="col0.05list.csv")

resSig2 <- subset(resOrdered, padj <= 0.1)
resSig2
write.csv(resSig2, file="col0.1list.csv")

resSig3 <- subset(resSig1, log2FoldChange >=1 | log2FoldChange <=-1)
resSig3
write.csv(resSig3, file="col0.05fold2list.csv")

resSig4 <- subset(resSig2, log2FoldChange >=1 | log2FoldChange <=-1)
resSig4
write.csv(resSig4, file="col0.1fold2list.csv")


##psk##

countData <- read.csv("RNASeq_psk365_vs_buffer.csv", header = TRUE, sep = ",")
head(countData)
metaData <- read.csv("Metadata_psk365_vs_buffer.csv", header = TRUE, sep = ",")
metaData
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~Treatment, tidy = TRUE)
dds

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$Treatment <- relevel(dds$Treatment, ref = "Buffer")

dds <- DESeq(dds)
countsdds <- counts(dds, normalized = TRUE)
write.csv(countsdds, file="dds_norm_countspsk.csv")
summary(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res)

res <- results(dds, name="Treatment_Bacteria_vs_Buffer")
res
resultsNames(dds)

#Double check with contrast
res1 <- results( dds, contrast = c("Treatment", "Bacteria", "Buffer") )
res1

resLFC <- lfcShrink(dds, coef="Treatment_Bacteria_vs_Buffer", type="apeglm")
resLFC

register(SnowParam(4))
resOrdered <- res[order(res$pvalue),]
summary(res)

write.csv(as.data.frame(resOrdered), 
          file="Treatment_Bacteria_vs_Buffer_resultspsk.csv")

resSig1 <- subset(resOrdered, padj <= 0.05)
resSig1
write.csv(resSig1, file="psk0.05list.csv")

resSig2 <- subset(resOrdered, padj <= 0.1)
resSig2
write.csv(resSig2, file="psk0.1list.csv")

resSig3 <- subset(resSig1, log2FoldChange >=1 | log2FoldChange <=-1)
resSig3
write.csv(resSig3, file="psk0.05fold2list.csv")

resSig4 <- subset(resSig2, log2FoldChange >=1 | log2FoldChange <=-1)
resSig4
write.csv(resSig4, file="psk0.1fold2list.csv")

