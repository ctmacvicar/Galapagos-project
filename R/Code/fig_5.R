# ============================================================
# Tutorial on drawing an NMDS plot using ggplot2 adopted by
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================

#loading necessary libraries
library(phyloseq)
library(ggplot2)
library(vegan)
library(tidyr)
library(tidyverse)
rm(list=ls())


#loading data
taxa=read.csv("data/Oct_Galapagos_BACT_L7.csv",row.names=1,check.names=FALSE);
#taxa=read.table("data/Galapagos16S_L7.txt", sep="\t", header=T);
metadata=read.table("data/GalapagosMeta.txt", sep="\t", header=T);


subsetmetadata <-subset(metadata, sample!="M5_C");
metadata <- subsetmetadata
#subset_taxa<-subset(taxa, index!="M5_C");
#taxa = subset_taxa 

otus <- as.data.frame(t(taxa))
otus <- select(otus, -M5_C)


metadata=metadata[order(metadata$sample),];

metadata$Site<-factor(metadata$Site, levels=c("Mirador", "Cerro Alto","el junco"));
metadata$Native <- factor(metadata$Native, labels = c("no", "yes"));
#metadata$Sample_Depth <- factor(metadata$Sample_Depth, labels = c("rhizo", "5cm", "15cm"));


#Just a check to ensure that the samples in meta_table are in the same order as in abund_table
#meta_table<-metadata[rownames(abund_table),]
head(metadata)
head(otus)

Site <-metadata$Site
Native <-metadata$Native
Sample_Depth <- metadata$Sample_Depth
Site_Plant <- metadata$Site_Plant


#transposing columns 'cause vegan likes it that way
t_otus <- as.data.frame(t(otus))

#rarefying data - there's some controversy about that so we might choose to normalize it instead in the future.
min_depth = min(colSums(otus))
t_otus_rarefied <- as.data.frame(round(rrarefy(t_otus, min_depth)))

#choose which method to use
otus_dist = as.matrix((vegdist(t_otus, "bray")))

#perform NMDS
NMDS = metaMDS(otus_dist, distance = "bray", trymax = 100)
stressplot(NMDS) 

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Site <-metadata$Site, Sample_Depth <- metadata$Sample_Depth, Plant = metadata$Plant)

head(NMDS)

#plot using ggplot2.it’s an excellent package for plotting and comes with a ton of functionality.  
ggplot(NMDS, aes(x=MDS1, y=MDS2, fill = Site)) +
  geom_point(mapping=aes(shape = Sample_Depth, fill = Site, color = Site_Plant),
             size = 4,
  )+
  stat_ellipse() +
  theme_bw() +
  scale_shape_manual(values = c(23, 24, 22)) +
  scale_color_manual(values = c("Mirador_Guava" = "yellow", "Mirador_Scalesia" = "orange2", "Cerro Alto_Guava" = "seagreen1", "Cerro Alto_Vervain" = "#008E51", "el junco_guava" = "skyblue", "el junco_Miconia" = "#043FFF")) +
  scale_fill_manual(values = c("Mirador" = "orange2", "Cerro Alto" = "#008E51", "el junco" = "#043FFF")) +
  labs(title = "BACT NMDS Site Plot", x = "NMDS1", y="NMDS2")

ggsave(filename="02_galapagos16S_NMDS_sites.pdf",
       device="pdf",path="./images",
       width=6,
       height=5,
       units="in",
       dpi=500);
