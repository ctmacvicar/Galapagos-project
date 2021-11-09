## CODE FOR: generating NMDS plots based on ITS microbial commmunity dissimilarity matrices


source(file="scripts/00_background.R"); #load necessary packages and specifications
library(ggplot2)
library(vegan)
library(dplyr)
rm(list=ls())

#Loading Data
taxa=read.table("data/R_Oct_Galapagos_L6_ITS.txt", sep="\t", header=T);
meta_data=read.table("data/GalapagosMeta.txt", sep="\t", header=T);
meta_data=meta_data[order(meta_data$sample),];

#removing the samples in the metadata that aren't in the taxa data file
remove_outlier <- (meta_data[-c(131, 130,2,17),])
meta_data <- remove_outlier

#taxa dataset for permanovas
taxa_mirador <- taxa %>% dplyr::select (M1_A, M1_B, M1_C, M10_A, 
                                        M10_B, M10_C, M11_A, M11_B, M11_C, M12_A, M12_B, M12_C, M2_A, M2_B, 
                                        M2_C, M3_A, M3_B, M3_C, M4_A, M4_B, M4_C, M5_A, M5_B,
                                        M6_B, M6_C, M7_A, M7_B, M7_C, M8_A, M8_B, M8_C, M9_A, M9_B, M9_C)

taxa_mirador_t<-t(taxa_mirador)

taxa_cerroA <- taxa %>% dplyr::select(C1_A, C10_A, C10_B, 
                                      C10_C, C11_A, C11_B, C11_C, C12_A, C12_B, C12_C, C2_A, C2_B, C2_C, C3_A, 
                                      C3_B, C4_A, C4_B, C4_C, C5_A, C5_B, C5_C, C6_B, C6_C, C7_A, C7_B, 
                                      C7_C, C8_A, C8_B, C8_C, C9_A, C9_B, C9_C);

taxa_cerroA_t= t(taxa_cerroA)

taxa_eljunco <- taxa %>% dplyr::select(J1_A, J1_B, J1_C, J10_A, 
                                       J10_B, J10_C, J11_A, J11_B, J11_C, J12_A, J12_B, J12_C, J13_A, J13_B, 
                                       J13_C, J14_A, J14_B, J14_C, J15_A, J15_B, J15_C, J16_A, J16_B, J16_C, 
                                       J17_A, J17_B, J17_C, J18_A, J18_B, J18_C, J19_A, J19_B, J19_C, J2_A,
                                       J2_B, J2_C, J20_A, J20_B, J20_C, J21_A, J21_B, J21_C, J22_A, J22_B, J22_C,
                                       J23_A, J23_B, J23_C, J24_A, J24_B, J24_C, J3_A, J3_B, J3_C, J4_A, J4_B, J4_C, J5_A, J5_B, J5_C, J6_A, J6_B, J6_C, J7_A, J7_B, J7_C, J8_A, J8_B, 
                                       J8_C, J9_A, J9_B, J9_C)

taxa_eljunco_t<- t(taxa_eljunco)

# Filtering Meta Data
meta_data$Plant_Invasiveness <- factor(meta_data$Plant_Invasiveness, labels = c("Native", "Invasive"));

meta_mirador_initial <- meta_data %>% filter(Site == "Mirador");

meta_cerro_A <- meta_data %>% filter(Site == "Cerro Alto");

meta_eljunco <- meta_data %>% filter(Site == "el junco");


################################################################################
#                               Mirador NMDS                                       #
################################################################################

head(meta_mirador)
head(taxa_mirador)


Plant_Invasiveness <-meta_data$Plant_Invasiveness
Sample_Depth <- meta_data$Sample_Depth
Plant <- meta_data$Plant


#transposing columns 'cause vegan likes it that way
t_otus <- as.data.frame(t(taxa_mirador))
t_otus[rowSums(t_otus[])>0,]

#rarefying data - there's some controversy about that so we might choose to normalize it instead in the future.
min_depth = min(colSums(taxa_mirador))
t_otus_rarefied <- as.data.frame(round(rrarefy(taxa_mirador, min_depth)))

#choose which method to use
otus_dist = as.matrix((vegdist(taxa_mirador_t, "bray")))
otus_dist[rowSums(otus_dist[])>0,]

#perform NMDS-- this is how far I got
NMDS = metaMDS(otus_dist, distance = "bray", trymax = 100)
stressplot(NMDS) 


#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, Plant <-meta_mirador_initial$Plant)
NMDS1 = data.frame(MDS1=MDS1, Sample_Depth <- meta_mirador_initial$Sample_Depth)

head(NMDS1)
head(NMDS)

NMDS2 = inner_join(NMDS1, NMDS)

#plot using ggplot2.it’s an excellent package for plotting and comes with a ton of functionality.  

box_col=c("#FF9200","#FF9200","#FF9200" )


mirador_NMDS <- ggplot(NMDS2, aes(x=MDS1, y=MDS2, fill=Plant, group= c(Plant))) +
  geom_point(mapping=aes(shape=as.factor(Sample_Depth)),
             size = 4, color="#FF9200")+
  scale_shape_manual(values = c(23, 24, 22)) +
  stat_ellipse() +
  theme_bw() +
  scale_fill_manual(
    values = c("white", "#FF9200"),
    guide = guide_legend(override.aes = list(shape = 21)))+
  scale_color_manual(values = c("#FF9200")) +
  labs(title = "NMDS Site Plot", shape= "Sample Depth", fill= "Plant")
plot(mirador_NMDS)

#save NMDS plot
ggsave(filename="galapagosITS_alldepths_mirador_NMDS_oct.pdf",
       device="pdf",path="./images",
       plot=mirador_NMDS,
       width=6,
       height=5,
       units="in",
       dpi=500);

################################################################################
#                         Cerro Alto NMDS                                        #
################################################################################
head(meta_cerro_A)
head(taxa_cerro_A)


Plant_Invasiveness <-meta_data$Plant_Invasiveness
Plant <- meta_data$Plant


#transposing columns 'cause vegan likes it that way
t_otus <- as.data.frame(t(taxa_cerroA))
t_otus[rowSums(t_otus[])>0,]

#rarefying data
min_depth = min(colSums(taxa_cerroA))
t_otus_rarefied <- as.data.frame(round(rrarefy(taxa_cerroA, min_depth)))
#tt_otus_rarefied <- as.data.frame(t(t_otus_rarefied))


#choose which method to use
otus_dist = as.matrix((vegdist(taxa_cerroA_t, "bray")))
otus_dist[rowSums(otus_dist[])>0,]

#perform NMDS
NMDS = metaMDS(otus_dist, distance = "bray", trymax = 100)
stressplot(NMDS) 


#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, Plant <-meta_cerro_A$Plant)
NMDS1 = data.frame(MDS1=MDS1, Sample_Depth <- meta_cerro_A$Sample_Depth)

head(NMDS1)
head(NMDS)

NMDS2 = inner_join(NMDS1, NMDS)


#plot using ggplot2.it’s an excellent package for plotting and comes with a ton of functionality.  

box_col=c("#008E51","#008E51")

cerro_A_NMDS <- ggplot(NMDS2, aes(x=MDS1, y=MDS2, fill=Plant, group= c(Plant))) +
  geom_point(mapping=aes(shape=as.factor(Sample_Depth)),
             size = 4, color="#008E51")+
  scale_shape_manual(values = c(23, 24, 22)) +
  stat_ellipse() +
  theme_bw() +
  scale_fill_manual(
    values = c("white", "#008E51"),
    guide = guide_legend(override.aes = list(shape = 21)))+
  scale_color_manual(values = c("#008E51")) +
  labs(title = "NMDS Site Plot", shape= "Sample Depth", fill= "Plant")

#view plot
plot(cerro_A_NMDS)

#save plot
ggsave(filename="galapagosITS_Cerro_A_alldepths_NMDS_october.pdf",
       device="pdf",path="./images",
       plot=cerro_A_NMDS,
       width=6,
       height=5,
       units="in",
       dpi=500);

################################################################################
#                           el junco NMDS                                     #
################################################################################

head(meta_mirador)
head(taxa_eljunco)


#transposing columns 'cause vegan likes it that way
t_otus <- as.data.frame(t(taxa_eljunco))
t_otus[rowSums(t_otus[])>0,]

#rarefying data - there's some controversy about that so we might choose to normalize it instead in the future.
min_depth = min(colSums(taxa_eljunco))
t_otus_rarefied <- as.data.frame(round(rrarefy(taxa_eljunco, min_depth)))
#tt_otus_rarefied <- as.data.frame(t(t_otus_rarefied))

#choose which method to use
otus_dist = as.matrix((vegdist(taxa_eljunco_t, "bray")))
otus_dist[rowSums(otus_dist[])>0,]

#perform NMDS
NMDS = metaMDS(otus_dist, distance = "bray", trymax = 100)
stressplot(NMDS) 


#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, Plant <-meta_eljunco$Plant)
NMDS1 = data.frame(MDS1=MDS1, Sample_Depth <- meta_eljunco$Sample_Depth)

head(NMDS1)
head(NMDS)

NMDS2 = inner_join(NMDS1, NMDS)


#plot using ggplot2.it’s an excellent package for plotting and comes with a ton of functionality.  

box_col=c("#043FFF","#043FFF")

eljunco_NMDS <- ggplot(NMDS2, aes(x=MDS1, y=MDS2, fill=Plant, group= c(Plant))) +
  geom_point(mapping=aes(shape=as.factor(Sample_Depth)),
             size = 4, color="#043FFF")+
  scale_shape_manual(values = c(23, 24, 22)) +
  stat_ellipse() +
  theme_bw() +
  scale_fill_manual(
    values = c("white", "#043FFF"),
    guide = guide_legend(override.aes = list(shape = 21)))+
  scale_color_manual(values = c("#043FFF")) +
  labs(title = "NMDS Site Plot", shape= "Sample Depth", fill= "Plant")

#view plot 
plot(eljunco_NMDS)

#save plot
ggsave(filename="galapagosITS_eljunco_alldepths_NMDS_oct.pdf",
       device="pdf",path="./images",
       plot=eljunco_NMDS,
       width=6,
       height=5,
       units="in",
       dpi=500);
