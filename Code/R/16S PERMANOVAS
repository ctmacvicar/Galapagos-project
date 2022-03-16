#16S October 2021 Permanovas for Galapagos data

source(file="scripts/00_background.R"); #load necessary packages and specifications
library(ggplot2)
library(vegan)
library(dplyr)
rm(list=ls())

#Loading Data
taxa=read.table("data/R_Oct_Galapagos_L6_16S.txt",sep="\t", header=T);
meta_data=read.table("data/GalapagosMeta.txt", sep="\t", header=T);
meta_data=meta_data[order(meta_data$sample),];


#remove taxa in meta data that aren't in taxa file
#remove_outlier <- (meta_data[-c(2,17,131),])
#meta_data <- remove_outlier

#16S taxa at all sample depths
taxa_mirador_perm <- taxa %>% dplyr::select (M1_A, M1_B, M1_C, M10_A, 
                                              M10_B, M10_C, M11_A, M11_B, M11_C, M12_A, M12_B, M12_C, M2_A, M2_B, 
                                              M2_C, M3_A, M3_B, M3_C, M4_A, M4_B, M4_C, M5_A, M5_B, M5_C,
                                              M6_B, M6_C, M7_A, M7_B, M7_C, M8_A, M8_B, M8_C, M9_A, M9_B, M9_C)
taxa_cerroA_perm <- taxa %>% dplyr::select(C1_A, C10_A, C10_B, 
                                            C10_C, C11_A, C11_B, C11_C, C12_A, C12_B, C12_C, C2_A, C2_B, C2_C, C3_A, 
                                            C3_B, C4_A, C4_B, C4_C, C5_A, C5_B, C5_C, C6_B, C6_C, C7_A, C7_B, 
                                            C7_C, C8_A, C8_B, C8_C, C9_A, C9_B, C9_C);
taxa_eljunco_perm <- taxa %>% dplyr::select(J1_A, J1_B, J1_C, J10_A, 
                                             J10_B, J10_C, J11_A, J11_B, J11_C, J12_A, J12_B, J12_C, J13_A, J13_B, 
                                             J13_C, J14_A, J14_B, J14_C, J15_A, J15_B, J15_C, J16_A, J16_B, J16_C, 
                                             J17_A, J17_B, J17_C, J18_A, J18_B, J18_C, J19_A, J19_B, J19_C, J2_A,
                                             J2_B, J2_C, J20_A, J20_B, J20_C, J21_A, J21_B, J21_C, J22_A, J22_B, J22_C,
                                             J23_A, J23_B, J23_C, J24_A, J24_B, J24_C, J3_A, J3_B, J3_C, J4_A, J4_B, J4_C, J5_A, J5_B, J5_C, J6_A, J6_B, J6_C, J7_A, J7_B, J7_C, J8_A, J8_B, 
                                             J8_C, J9_A, J9_B, J9_C)

taxa_all <- taxa  %>% dplyr::select (M1_A, M1_B, M1_C, M10_A, 
                                     M10_B, M10_C, M11_A, M11_B, M11_C, M12_A, M12_B, M12_C, M2_A, M2_B, 
                                     M2_C, M3_A, M3_B, M3_C, M4_A, M4_B, M4_C, M5_A, M5_B, M5_C,
                                     M6_B, M6_C, M7_A, M7_B, M7_C, M8_A, M8_B, M8_C, M9_A, M9_B, M9_C, C1_A, C10_A, C10_B, 
                                           C10_C, C11_A, C11_B, C11_C, C12_A, C12_B, C12_C, C2_A, C2_B, C2_C, C3_A, 
                                           C3_B, C4_A, C4_B, C4_C, C5_A, C5_B, C5_C, C6_B, C6_C, C7_A, C7_B, 
                                           C7_C, C8_A, C8_B, C8_C, C9_A, C9_B, C9_C, J1_A, J1_B, J1_C, J10_A, 
                                            J10_B, J10_C, J11_A, J11_B, J11_C, J12_A, J12_B, J12_C, J13_A, J13_B, 
                                            J13_C, J14_A, J14_B, J14_C, J15_A, J15_B, J15_C, J16_A, J16_B, J16_C, 
                                            J17_A, J17_B, J17_C, J18_A, J18_B, J18_C, J19_A, J19_B, J19_C, J2_A,
                                            J2_B, J2_C, J20_A, J20_B, J20_C, J21_A, J21_B, J21_C, J22_A, J22_B, J22_C,
                                            J23_A, J23_B, J23_C, J24_A, J24_B, J24_C, J3_A, J3_B, J3_C, J4_A, J4_B, J4_C, J5_A, J5_B, J5_C, J6_A, J6_B, J6_C, J7_A, J7_B, J7_C, J8_A, J8_B, 
                                            J8_C, J9_A, J9_B, J9_C)

# Filtering Meta Data by Site
meta_data$Plant_Invasiveness <- factor(meta_data$Plant_Invasiveness, labels = c("Native", "Invasive"));

meta_mirador <- meta_data %>% filter(Site == "Mirador");

meta_cerro_A <- meta_data %>% filter(Site == "Cerro Alto");

meta_eljunco <- meta_data %>% filter(Site == "el junco");


#Mirador
################################################################################
#             2. Generate the 2 types of distance matrices for permanova
################################################################################
#transpose OTU table so OTUs are columns
taxadf=taxa_mirador_perm
taxadf=t(taxadf);
taxadf=taxadf[order(row.names(taxadf)),];

###BRAY-CURTIS distance
bray<-apply(taxadf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
#print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distance
jac=(taxadf>0)*1;
#print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");

#Mirador
################################################################################
#             3. Conduct PERMANOVAs for Plant, Sample Depth, and interaction
################################################################################

#Does Plant type influence microbial communities within Mirador?
#Bray-Curtis
print(adonis(bray.dist~Plant,             
             data=meta_mirador, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Plant,             
             data=meta_mirador, method = "jaccard",         
             permutations = 999));

#Does Sampling Depth influence microbial communities within Mirador?
#Bray-Curtis
print(adonis(bray.dist~Sample_Depth,             
             data=meta_mirador, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Sample_Depth,             
             data=meta_mirador, method = "jaccard",         
             permutations = 999));

# Does interaction between plant and depth microbial communities within Mirador?
#Bray-Curtis
print(adonis(bray.dist~(Plant*Sample_Depth),             
             data=meta_mirador, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~(Plant*Sample_Depth),             
             data=meta_mirador, method = "jaccard",         
             permutations = 999));


#Cerro Alto
################################################################################
#             2. Generate the 2 types of distance matrices for permanova
################################################################################
#transpose OTU table so OTUs are columns
taxadf=taxa_cerroA_perm
taxadf=t(taxadf);
taxadf=taxadf[order(row.names(taxadf)),];

###BRAY-CURTIS distance
bray<-apply(taxadf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
#print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distance
jac=(taxadf>0)*1;
#print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");



#Cerro Alto
################################################################################
#             3. Conduct PERMANOVAs for Plant, Sample Depth, and interaction
################################################################################

#Does Plant type influence microbial communities in Cerro Alto?
#Bray-Curtis
print(adonis(bray.dist~Plant,             
             data=meta_cerro_A, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Plant,             
             data=meta_cerro_A, method = "jaccard",         
             permutations = 999));

#Does Sampling Depth influence microbial communities in Cerro Alto?
#Bray-Curtis
print(adonis(bray.dist~Sample_Depth,             
             data=meta_cerro_A, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Sample_Depth,             
             data=meta_cerro_A, method = "jaccard",         
             permutations = 999));

# Does interaction between plant and depth influence microbial communities in Cerro Alto?
#Bray-Curtis
print(adonis(bray.dist~(Plant*Sample_Depth),             
             data=meta_cerro_A, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~(Plant*Sample_Depth),             
             data=meta_cerro_A, method = "jaccard",         
             permutations = 999));


# El Junco
################################################################################
#             2. Generate the 2 types of distance matrices for permanova
################################################################################
#transpose OTU table so OTUs are columns
taxadf=taxa_eljunco_perm
taxadf=t(taxadf);
taxadf=taxadf[order(row.names(taxadf)),];

###BRAY-CURTIS distance
bray<-apply(taxadf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
#print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distance
jac=(taxadf>0)*1;
#print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");




#El Junco
################################################################################
#             3. Conduct PERMANOVAs for Plant, Sample Depth, and interaction
################################################################################

#Does Plant type influence microbial communities in El Junco?
#Bray-Curtis
print(adonis(bray.dist~Plant,             
             data=meta_eljunco, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Plant,             
             data=meta_eljunco, method = "jaccard",         
             permutations = 999));

#Does Sampling Depth influence microbial communities in El Junco?
#Bray-Curtis
print(adonis(bray.dist~Sample_Depth,             
             data=meta_eljunco, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Sample_Depth,             
             data=meta_eljunco, method = "jaccard",         
             permutations = 999));

# Does interaction between plant and depth influence microbial communities in El Junco?
#Bray-Curtis
print(adonis(bray.dist~(Plant*Sample_Depth),             
             data=meta_eljunco, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~(Plant*Sample_Depth),             
             data=meta_eljunco, method = "jaccard",         
             permutations = 999));


# All Sites Together
################################################################################
#             2. Generate the 2 types of distance matrices for permanova
################################################################################
#transpose OTU table so OTUs are columns
taxadf=taxa_all
taxadf=t(taxadf);
taxadf=taxadf[order(row.names(taxadf)),];

###BRAY-CURTIS distance
bray<-apply(taxadf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
#print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distance
jac=(taxadf>0)*1;
#print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");




#All Sites
################################################################################
#             3. Conduct PERMANOVAs for Plant, Sample Depth, and interaction
################################################################################

#Does Site influence microbial communities?
#Bray-Curtis
print(adonis(bray.dist~Site,             
             data=meta_data, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Site,             
             data=meta_data, method = "jaccard",         
             permutations = 999));





################################################################################
#Filtering Data by depth-- only 5 cm and rhizosphere-- in all sites

#rhizo and 5 cm sample depths only in meta data

sample_depth <- c("rhizo", "5cm")
meta_rhizo5_mirador <- meta_mirador %>% filter (Sample_Depth %in% sample_depth)

meta_cerro_A <- meta_data %>% filter(Site == "Cerro Alto");
meta_rhizo5_cerro_A <- meta_cerro_A %>% filter(Sample_Depth %in% sample_depth);

meta_eljunco <- meta_data %>% filter(Site == "el junco");
meta_rhizo5_eljunco <- meta_eljunco %>% filter(Sample_Depth %in% sample_depth);

meta_data_rhizo5 <- meta_data %>% filter(Sample_Depth %in% sample_depth)

#Mirador taxa (rhizo and 5 cm only)
taxa_rhizo5_mirador <- taxa %>% dplyr::select(M1_A, M1_B, M10_A, M10_B, M11_B, M12_B,
                                              M11_A, M12_A, M2_A, M2_B, M3_A, M3_B, M4_A, M4_B, M5_A, 
                                              M5_B, M6_B, M7_A, M7_B, M8_A, M8_B, M9_A, M9_B, M6_C);

#Cerro Alto taxa (rhizo and 5 cm only)
taxa_rhizo5_cerro_A <- taxa %>% dplyr::select(C1_A, C10_A, C10_B, C11_A, C11_B, C12_A, C12_B, C2_A, C2_B, C3_A, 
                                              C3_B, C4_B, C4_A, C5_A, C5_B, C6_B, C7_A, C7_B, C8_A, C8_B, C9_A, C9_B);

# el junco taxa (rhizo and 5 cm only)
taxa_rhizo5_eljunco <- taxa %>% dplyr::select(J1_A, J1_B, J10_A, J10_B, J11_A, J11_B, J12_A, J12_B, J13_A, J13_B, J14_A, J14_B,
                                              J15_A, J15_B, J16_A, J16_B, J17_B, J18_B, J19_B, J20_B, J21_B, J22_B, J23_B, J24_B,
                                              J17_A, J18_A, J19_A, J20_A, J21_A, J22_A, J2_A, J2_B,
                                              J23_A, J24_A, J3_A, J3_B, J4_B, J4_A,, J5_A, J5_B, J6_A, J6_B,  J7_A, J7_B, J8_A, J8_B, J9_A, J9_B);

taxa_all_rhizo5 <- taxa %>% dplyr::select(M1_A, M1_B, M10_A, M10_B, M11_B, M12_B,
                                          M11_A, M12_A, M2_A, M2_B, M3_A, M3_B, M4_A, M4_B, M5_A, 
                                          M5_B, M6_B, M7_A, M7_B, M8_A, M8_B, M9_A, M9_B, M6_C, C1_A, C10_A, C10_B, C11_A, C11_B, C12_A, C12_B, C2_A, C2_B, C3_A, 
                                              C3_B, C4_B, C4_A, C5_A, C5_B, C6_B, C7_A, C7_B, C8_A, C8_B, C9_A, C9_B, J1_A, J1_B, J10_A, J10_B, J11_A, J11_B, J12_A, J12_B, J13_A, J13_B, J14_A, J14_B,
                                              J15_A, J15_B, J16_A, J16_B, J17_B, J18_B, J19_B, J20_B, J21_B, J22_B, J23_B, J24_B,
                                              J17_A, J18_A, J19_A, J20_A, J21_A, J22_A, J2_A, J2_B,
                                              J23_A, J24_A, J3_A, J3_B, J4_B, J4_A,, J5_A, J5_B, J6_A, J6_B,  J7_A, J7_B, J8_A, J8_B, J9_A, J9_B);
  
################################################################################
#                               PERMANOVAS
################################################################################


#Mirador-- only 5 cm and rhizosphere
################################################################################
#             2. Generate the 2 types of distance matrices for permanova
################################################################################
#transpose OTU table so OTUs are columns
taxadf=taxa_rhizo5_mirador
taxadf=t(taxadf);
taxadf=taxadf[order(row.names(taxadf)),];

###BRAY-CURTIS distance
bray<-apply(taxadf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
#print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distance
jac=(taxadf>0)*1;
#print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");

#Mirador-- only 5 cm and rhizosphere
################################################################################
#             3. Conduct PERMANOVAs for Plant, Sample Depth, and interaction
################################################################################

#Does Plant type influence microbial communities in Mirador?
#Bray-Curtis
print(adonis(bray.dist~Plant,             
             data=meta_rhizo5_mirador, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Plant,             
             data=meta_rhizo5_mirador, method = "jaccard",         
             permutations = 999));

#Does Sampling Depth influence microbial communities in Mirador?
#Bray-Curtis
print(adonis(bray.dist~Sample_Depth,             
             data=meta_rhizo5_mirador, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Sample_Depth,             
             data=meta_rhizo5_mirador, method = "jaccard",         
             permutations = 999));

# Does interaction between plant and depth influence microbial communities in Mirador?
#Bray-Curtis
print(adonis(bray.dist~(Plant*Sample_Depth),             
             data=meta_rhizo5_mirador, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~(Plant*Sample_Depth),             
             data=meta_rhizo5_mirador, method = "jaccard",         
             permutations = 999));


#Cerro Alto-- only 5 cm and rhizosphere
################################################################################
#             2. Generate the 2 types of distance matrices for permanova
################################################################################
#transpose OTU table so OTUs are columns
taxadf=taxa_rhizo5_cerro_A
taxadf=t(taxadf);
taxadf=taxadf[order(row.names(taxadf)),];

###BRAY-CURTIS distance
bray<-apply(taxadf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
#print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distance
jac=(taxadf>0)*1;
#print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");



#Cerro Alto-- only 5 cm and rhizosphere
################################################################################
#             3. Conduct PERMANOVAs for Plant, Sample Depth, and interaction
################################################################################

#Does Plant type influence microbial communities in Cerro Alto?
#Bray-Curtis
print(adonis(bray.dist~Plant,             
             data=meta_rhizo5_cerro_A, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Plant,             
             data=meta_rhizo5_cerro_A, method = "jaccard",         
             permutations = 999));

#Does Sampling Depth influence microbial communities in Cerro Alto?
#Bray-Curtis
print(adonis(bray.dist~Sample_Depth,             
             data=meta_rhizo5_cerro_A, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Sample_Depth,             
             data=meta_rhizo5_cerro_A, method = "jaccard",         
             permutations = 999));

# Does interaction between plant and depth influence microbial communities in Cerro Alto?
#Bray-Curtis
print(adonis(bray.dist~(Plant*Sample_Depth),             
             data=meta_rhizo5_cerro_A, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~(Plant*Sample_Depth),             
             data=meta_rhizo5_cerro_A, method = "jaccard",         
             permutations = 999));


# El Junco-- only 5 cm and rhizosphere
################################################################################
#             2. Generate the 2 types of distance matrices for permanova
################################################################################
#transpose OTU table so OTUs are columns
taxadf=taxa_rhizo5_eljunco
taxadf=t(taxadf);
taxadf=taxadf[order(row.names(taxadf)),];

###BRAY-CURTIS distance
bray<-apply(taxadf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
#print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distance
jac=(taxadf>0)*1;
#print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");


#El Junco-- only 5 cm and rhizosphere
################################################################################
#             3. Conduct PERMANOVAs for Plant, Sample Depth, and interaction
################################################################################

#Does Plant type influence microbial communities in El Junco?
#Bray-Curtis
print(adonis(bray.dist~Plant,             
             data=meta_rhizo5_eljunco, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Plant,             
             data=meta_rhizo5_eljunco, method = "jaccard",         
             permutations = 999));

#Does Sampling Depth influence microbial communities in El Junco?
#Bray-Curtis
print(adonis(bray.dist~Sample_Depth,             
             data=meta_rhizo5_eljunco, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Sample_Depth,             
             data=meta_rhizo5_eljunco, method = "jaccard",         
             permutations = 999));

# Does interaction between plant and depth influence microbial communities in El Junco?
#Bray-Curtis
print(adonis(bray.dist~(Plant*Sample_Depth),             
             data=meta_rhizo5_eljunco, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~(Plant*Sample_Depth),             
             data=meta_rhizo5_eljunco, method = "jaccard",         
             permutations = 999));

# All Sites-- only 5 cm and rhizosphere
################################################################################
#             2. Generate the 2 types of distance matrices for permanova
################################################################################
#transpose OTU table so OTUs are columns
taxadf=taxa_all_rhizo5
taxadf=t(taxadf);
taxadf=taxadf[order(row.names(taxadf)),];

###BRAY-CURTIS distance
bray<-apply(taxadf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
#print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distance
jac=(taxadf>0)*1;
#print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");


#All Sites-- only 5 cm and rhizosphere
################################################################################
#             3. Conduct PERMANOVAs for Plant, Sample Depth, and interaction
################################################################################

#Does Site influence microbial communities?
#Bray-Curtis
print(adonis(bray.dist~Site,             
             data=meta_data_rhizo5, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Site,             
             data=meta_data_rhizo5, method = "jaccard",         
             permutations = 999));


#Comparing 15 cm to rhizo+5cm sample depth
################################################################################
#Reading and Filtering data

taxa <- taxa %>% dplyr::select (M1_A, M1_B, M1_C, M10_A, 
                                M10_B, M10_C, M11_A, M11_B, M11_C, M12_A, M12_B, M12_C, M2_A, M2_B, 
                                M2_C, M3_A, M3_B, M3_C, M4_A, M4_B, M4_C, M5_A, M5_B, 
                                M6_B, M6_C, M7_A, M7_B, M7_C, M8_A, M8_B, M8_C, M9_A, M9_B, M9_C, C1_A, C10_A, C10_B, 
                                C10_C, C11_A, C11_B, C11_C, C12_A, C12_B, C12_C, C2_A, C2_B, C2_C, C3_A, 
                                C3_B,C4_A, C4_B, C4_C, C5_A, C5_B, C5_C, C6_B, C6_C, C7_A, C7_B, 
                                C7_C, C8_A, C8_B, C8_C, C9_A, C9_B, C9_C, J1_A, J1_B, J1_C, J10_A, 
                                J10_B, J10_C, J11_A, J11_B, J11_C, J12_A, J12_B, J12_C, J13_A, J13_B, 
                                J13_C, J14_A, J14_B, J14_C, J15_A, J15_B, J15_C, J16_A, J16_B, J16_C, 
                                J17_A, J17_B, J17_C, J18_A, J18_B, J18_C, J19_A, J19_B, J19_C, J2_A,
                                J2_B, J2_C, J20_A, J20_B, J20_C, J21_A, J21_B, J21_C, J22_A, J22_B, J22_C,
                                J23_A, J23_B, J23_C, J24_A, J24_B, J24_C, J3_A, J3_B, J3_C, J4_A, J4_B, J4_C, J5_A, J5_B, J5_C, J6_A, J6_B, J6_C, J7_A, J7_B, J7_C, J8_A, J8_B, 
                                J8_C, J9_A, J9_B, J9_C)


#meta data for comparing 15 cm to rhizo+5cm through permanova
meta_data2 <- read.table("data/GalapagosMeta_15cm&r5.txt", sep="\t", header=T)

#remove taxa in meta data that aren't in taxa file
remove_outlier <- (meta_data2[-c(2,17,131, 130),])
meta_data2 <- remove_outlier


################################################################################
#             2. Generate the 2 types of distance matrices
################################################################################
taxadf=taxa
taxadf=t(taxadf);
taxadf=taxadf[order(row.names(taxadf)),];

###BRAY-CURTIS distance
bray<-apply(taxadf, 1, function(i) (i/sum(i)));
bray=as.data.frame(t(bray));
#print(rowSums(bray));
bray.dist=vegdist(bray, method="bray");

###JACCARD distance
jac=(taxadf>0)*1;
#print(rowSums(jac));
jac.dist=vegdist(jac, method="jaccard");


################################################################################
#             3. Conduct PERMANOVAs for Depth
################################################################################
#Does sample depth influence microbial communities all sites?
#Bray-Curtis
print(adonis(bray.dist~Sample_Depth,             
             data=meta_data2, method = "bray",         
             permutations = 999)); 
#Jaccard
print(adonis(jac.dist~Sample_Depth,             
             data=meta_data2, method = "jaccard",         
             permutations = 999));
