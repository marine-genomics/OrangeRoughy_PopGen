#### R script file to analyse vcf obtained with GBS for orange roughy ####

#### 1. create map of sampling sites ####

library("ggplot2")
library("cowplot")
library("ggthemes")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
library("ggOceanMaps")
library("readxl")
library("gdistance")
library("raster")


# make map
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)


africa <- ggplot(data=world) + theme_bw() + 
  geom_sf() + coord_sf(xlim=c(-6, 50), ylim=c(35,-33), datum=NA) +
  #   geom_rect(xmin = -19, xmax = 53.1, ymin = 38.5, ymax = -36, 
  #     fill = NA, colour = "black", size = 1.5) +
  geom_rect(xmin = 6.5, xmax = 19, ymin = -17, ymax = -36, 
            fill = NA, colour = "black", size = 0.75) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
africa

world_points<- st_centroid(world)
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))

seafo <- data.frame(longitude=c(5.8,6.47933333,9.04766667),latitude=c(-26.0771667,-24.6336667,-20.5755))
namibia_c <- data.frame(longitude=c(13.3832,10.1242,13.5988,12.7315,13.327),latitude=c(-24.537,-19.3753,-26.36,-22.4928,-24.6772))
south_africa <- data.frame(longitude=c(21.584167),latitude=c(-36.668333))

dt <- data.frame(lon = c(-15, 30), lat = c(-17, -38))

map <- basemap(data=dt, bathymetry=FALSE, land.col = "grey") +
  theme_bw() +
  geom_sf() 
annotation_scale(location = "bl", width_hint = 0.5) +
annotation_north_arrow(location = "bl", which_north = "true", 
pad_x = unit(0.5, "in"), pad_y = unit(0.25, "in"),
style = north_arrow_fancy_orienteering) +
xlab("Longitude") + ylab("Latitude") +
coord_sf(xlim = c(5.5, 30), ylim = c(-17, -37.5), expand = FALSE) + 
annotate(geom = "text", x = 16, y = -21, label = "Namibia", 
         fontface = "bold", color = "black", size = 6) +
(geom = "text", x = 23, y = -30, label = "South Africa", 
         fontface = "bold", color = "black", size = 6) +
(geom = "text", x = 19.55, y = -33.9, label = "Cape Town", 
           fontface = "plain", color = "black", size = 4) +
annotate(geom = "text", x = 16.2, y = -26.6, label = "Luderitz", 
         fontface = "plain", color = "black", size = 4) +
annotate(geom = "text", x = 18, y = -28.5, label = "Orange River", 
          fontface = "plain", color = "black", size = 4) +
annotate(geom = "text", x = 15.56, y = -22.95, label = "Walvis Bay", 
          fontface = "plain", color = "black", size = 4) +
annotate(geom = "text", x = 10, y = -36, label = "Atlantic Ocean", 
        fontface = "italic", color = "black", size = 6) +
annotate(geom = "text", x = 27, y = -36, label = "Indian Ocean", 
        fontface = "italic", color = "black", size = 6) + 
geom_point(data = seafo, aes(x = longitude, y = latitude), size = 5, 
           shape = 21, fill = "darkgrey") +
annotate(geom = "text", x = 27, y = -36, label = "SEAFO 26", 
           fontface = "italic", color = "black", size = 6) + 
geom_point(data = namibia_c, aes(x = longitude, y = latitude), size = 5, 
             shape = 21, fill ="darkgreen") +  
geom_point(data = south_africa, aes(x = longitude, y = latitude), size = 5, 
         shape = 21, fill = "darkgoldenrod2") +
annotation_custom(grob = ggplotGrob(africa),
                   xmin = 24,
                   xmax = 24 + (50-(-16))/10.5,
                   ymin = -22,
                   ymax = -22 + (35-(-33))/13.5)

ggsave("map.tiff")
raster::writeRaster(map, filename = "map.tiff", overwrite = TRUE)

col <- c("darkgrey","darkgreen","darkgoldenrod2")

#### 2. basic population analyses ####

library(vcfR)
library(dartR)
library(adegenet)
library(hierfstat)
library(readxl)
library(poppr)
library(tidyverse)
library(factoextra)
library(Relatedness)
library(StAMPP)
library(ape)
library(ggplot2)
library(vegan)
library(qqman)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(ggrepel)

#### import vcf and check ####
vcf <- file.choose(new=FALSE) #this opens an interactive window to select file
data_all_vcf <- read.vcfR(vcf,verbose=TRUE)
data_all <- vcfR2genlight(data_all_vcf)

data_all1 <- gl.compliance.check(data_all)

# recode individuals names
names_ind <- as.data.frame(data_all1$ind.names)
ind.recode <- read.csv("samples_names.csv", header=FALSE)
ind.recode2 <- cbind(names_ind,ind.recode)
write_csv(ind.recode2, "ind.recode2.csv") #don't forget to open it and delete the first row

data_all2 <- gl.recode.ind(data_all1, "ind.recode3.csv", recalc = FALSE, mono.rm = FALSE, verbose = NULL)

#### import metadata ####
metadata <- read_excel("metadata1.xlsx", sheet="popmap_3")


#### calculate stats ####

# overall
data_all2$pop <- as.factor(metadata$Overall)
as.data.frame(data_all2$ind.names)
unique(as.data.frame(data_all2$pop))
all_gi <- gl2gi(data_all2)

overall_stats<-basic.stats(all_gi)
ci_fis <- boot.ppfis(all_gi)
overall_ar <- allelic.richness(all_gi)

allele_freq_region <- pop.freq(all_gi,diploid=TRUE)

overall_HO <- as.data.frame(overall_stats$Ho)
overall_HE <- as.data.frame(overall_stats$Hs)
overall_FIS <- as.data.frame(overall_stats$Fis)
overall_AR <- as.data.frame(overall_ar)

write_csv(overall_HO, "overall_HO_new.csv")
write_csv(overall_HE, "overall_HE_new.csv")
write_csv(overall_FIS, "overall_FIS_new.csv")
write_csv(overall_AR, "overall_AR_new.csv")

data_pca_sc <- scaleGen(all_gi, NA.method="mean",scale=F)
data_pca <- dudi.pca(data_pca_sc, scale=F, nf= 10, scannf = F)

fviz_pca_ind(data_pca,label="none",habillage = all_gi$pop,axes=c(1,3),
             legend.title ="Region",mean.point=F,pointsize=4,
             addEllipses = FALSE, ellipse.level=0.5,ellipse.alpha = 0, title="PCA") + theme_classic()


# per region
data_all2$pop <- as.factor(metadata$Region)
as.data.frame(data_all2$ind.names)
unique(as.data.frame(data_all2$pop))
all_gi2 <- gl2gi(data_all2)

col_region <- c("limegreen","azure3","darkgoldenrod2","navyblue")

overall_stats2<-basic.stats(all_gi2)
ci_fis2 <- boot.ppfis(all_gi2)
overall_ar2 <- allelic.richness(all_gi2)

allele_freq_region <- pop.freq(all_gi2,diploid=TRUE)
combined_allele_freq_region2 <- do.call(rbind, lapply(allele_freq_region, as.vector))
combined_allele_freq_region3 <- as.data.frame(combined_allele_freq_region2)
colnames(combined_allele_freq_region3)[1] <- "MAA_Namibia"
colnames(combined_allele_freq_region3)[2] <- "MIA_Namibia"
colnames(combined_allele_freq_region3)[3] <- "MAA_SEAFO"
colnames(combined_allele_freq_region3)[4] <- "MIA_SEAFO"
colnames(combined_allele_freq_region3)[5] <- "MAA_South_Africa"
colnames(combined_allele_freq_region3)[6] <- "MIA_South_Africa"

combined_allele_freq_region3$diff_NAM_SEAFO <- combined_allele_freq_region3$MAA_Namibia - combined_allele_freq_region3$MAA_SEAFO
combined_allele_freq_region3$diff_NAM_SA <- combined_allele_freq_region3$MAA_Namibia - combined_allele_freq_region3$MAA_South_Africa
combined_allele_freq_region3$diff_SEAFO_SA <- combined_allele_freq_region3$MAA_SEAFO - combined_allele_freq_region3$MAA_South_Africa


# make heatmaps of allele frequencies
library(pheatmap)

combined_allele_freq_region6 <- cbind(combined_allele_freq_region3$MAA_SEAFO,combined_allele_freq_region3$MAA_Namibia,combined_allele_freq_region3$MAA_South_Africa)

pheatmap(combined_allele_freq_region6, cluster_rows=FALSE, cluster_cols=FALSE, 
         color = colorRampPalette(c("lightblue", "darkred"))(200),show_rownames = F,labels_col = c("SEAFO", "Namibia", "South Africa"),
         main = "Heatmap of Allele Frequencies: Region")

# calculate correlation tests between minor allele frequencies and locations
cor.test(combined_allele_freq_region3$MIA_Namibia,combined_allele_freq_region3$MIA_SEAFO, method="pearson")cor.test(MAA_South_Africacombined_allele_freq_region3$MIA_Namibia,combined_allele_freq_region3$MIA_SEAFO, method="pearson")

boxplot(combined_allele_freq_region3$MIA_Namibia,combined_allele_freq_region3$MIA_SEAFO, combined_allele_freq_region3$MIA_South_Africa)

# prepare the data to build graphics of diversity metrics per region
overall_HO2 <- as.data.frame(overall_stats2$Ho)
ho_all <- read.csv("overall_HO_new.csv", header = TRUE)
overall_HO2 <- cbind(overall_HO2,ho_all[,-2])
colnames(overall_HO2)[4] <- "SE_Atlantic"

mean(overall_HO2$SEAFO)
mean(overall_HO2$Namibia)
mean(overall_HO2$South_Africa)
mean(overall_HO2$SE_Atlantic)

overall_HO2_2 <- stack(overall_HO2)
colnames(overall_HO2_2)[1] <- "HO"
colnames(overall_HO2_2)[2] <- "Location"

overall_HE2 <- as.data.frame(overall_stats2$Hs)
he_all <- read.csv("overall_HE_new.csv", header = TRUE)
overall_HE2 <- cbind(overall_HE2,he_all[,-2])
colnames(overall_HE2)[4] <- "SE_Atlantic"

overall_HE2_2 <- stack(overall_HE2)
colnames(overall_HE2_2)[1] <- "HE"
colnames(overall_HE2_2)[2] <- "Location"

mean(overall_HE2$SEAFO)
mean(overall_HE2$Namibia)
mean(overall_HE2$South_Africa)
mean(overall_HE2$SE_Atlantic)

overall_FIS2 <- as.data.frame(overall_stats2$Fis)
fis_all <- read.csv("overall_FIS_new.csv", header = TRUE)
overall_FIS2 <- cbind(overall_FIS2,fis_all[,-2])
colnames(overall_FIS2)[4] <- "SE_Atlantic"
overall_FIS2_2 <- stack(overall_FIS2)
colnames(overall_FIS2_2)[1] <- "FIS"
colnames(overall_FIS2_2)[2] <- "Location"

mean(overall_FIS2$SEAFO, na.rm=TRUE)
mean(overall_FIS2$Namibia, na.rm=TRUE)
mean(overall_FIS2$South_Africa, na.rm=TRUE)
mean(overall_FIS2$SE_Atlantic, na.rm=TRUE)

overall_AR2 <- as.data.frame(overall_ar2)
overall_AR2 <- overall_AR2[,-1]
colnames(overall_AR2)[1] <- "Namibia"
colnames(overall_AR2)[2] <- "SEAFO"
colnames(overall_AR2)[3] <- "South_Africa"
ar_all <- read.csv("overall_AR_new.csv", header = TRUE)
overall_AR2 <- cbind(overall_AR2,ar_all[,2])
colnames(overall_AR2)[4] <- "SE_Atlantic"
overall_AR2_2 <- stack(overall_AR2)
colnames(overall_AR2_2)[1] <- "AR"
colnames(overall_AR2_2)[2] <- "Location"

mean(overall_AR2$SEAFO, na.rm=TRUE)
mean(overall_AR2$Namibia, na.rm=TRUE)
mean(overall_AR2$South_Africa, na.rm=TRUE)
mean(overall_AR2$SE_Atlantic, na.rm=TRUE)

overall_stats_region <- cbind(overall_HO2_2,overall_HE2_2$HE,overall_FIS2_2$FIS,overall_AR2_2$AR)
colnames(overall_stats_region)[3] <- "HE"
colnames(overall_stats_region)[4] <- "FIS"
colnames(overall_stats_region)[5] <- "AR"

str(overall_stats_region)

overall_stats_region$Location <- factor(overall_stats_region$Location, levels = c("SEAFO", "Namibia", "South_Africa","SE_Atlantic"))

measure <- gather(overall_stats_region, key = "measure", value = "values", c("HO", "HE", "FIS", "AR"))

facet.labels <- c("Observed Heterozygosity", "Expected Heterozygosity", "Inbreeding Coefficient", "Allelic Richness") #correct symbols
names(facet.labels) <- c("HO", "HE", "FIS", "AR")

facet <- ggplot(data = measure, aes(x = Location, y = values, fill = Location)) +
  geom_violin(scale = "area", draw_quantiles = c(0.5))+
  facet_wrap(measure~., scales = "free_y", nrow = 2, labeller = labeller(measure = facet.labels))+
  #labs(x = "Locations")+
  scale_fill_manual(values=c("azure3","limegreen","darkgoldenrod2","navyblue"))+
  theme_bw() + 
  theme(legend.position="none", #for no legend, otherwise write "right"/"bottom"/"top", etc 
        legend.key.size = unit(1, 'cm'), 
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10), 
        strip.text.x = element_text(size = 10), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), #if required element_text(size = 10)
        axis.text = element_text(size = 10)) #if required face = "italic"

#### calculate statistical differences between metrics and locations ####
anova_HO <- aov(HO~Location, data = overall_stats_region)
summary(anova_HO)
plot(anova_HO)

kw_HO <- kruskal.test(HO~Location, data = overall_stats_region)
wt_HO<- pairwise.wilcox.test(overall_stats_region$HO, overall_stats_region$Location,
                              p.adjust.method = "BH") #BH is equivalent to false discovery rate

anova_HE <- aov(HE~Location, data = overall_stats_region)
summary(anova_HE)
plot(anova_HE)

kw_HE <- kruskal.test(HE~Location, data = overall_stats_region)
wt_HE<- pairwise.wilcox.test(overall_stats_region$HE, overall_stats_region$Location,
                             p.adjust.method = "BH") #BH is equivalent to false discovery rate

anova_FIS <- aov(FIS~Location, data = overall_stats_region)
summary(anova_FIS)
plot(anova_HE)

kw_FIS <- kruskal.test(FIS~Location, data = overall_stats_region)
wt_FIS <- pairwise.wilcox.test(overall_stats_region$FIS, overall_stats_region$Location,
                             p.adjust.method = "BH") #BH is equivalent to false discovery rate


#### calculate population structure ####

# overall
data_pca_sc2 <- scaleGen(all_gi2, NA.method="mean",scale=F)
data_pca2 <- dudi.pca(data_pca_sc2, scale=F, nf= 10, scannf = F)

pca_region <- fviz_pca_ind(data_pca2,label="none",habillage = all_gi2$pop,axes=c(1,2),
             legend.title ="Region",mean.point=F,pointsize=4,palette = col_region,
             addEllipses = FALSE, ellipse.level=0.5,ellipse.alpha = 0, title="PCA per region") + theme_classic()

all_gl2 <- gi2gl(all_gi2)
region_fst <- stamppFst(all_gl2, nboots = 10000, percent = 95)

write.table(region_fst$Fsts, "region_fst.txt")
write.table(region_fst$Pvalues, "region_pvalues.txt")

# per location
data_all2$pop <- as.factor(metadata$Location)
check_ind <- as.data.frame(data_all2$ind.names)
check_pop <- as.data.frame(data_all2$pop)
check_ind_pop <- cbind(check_ind,check_pop)

gl2structure(data_all2, outfile ="all_indv.str", outpath ='.')
write_csv(data_all2_strc, "data_all2_structure.csv")

unique(as.data.frame(data_all2$pop))
all_gi3a <- gl2gi(data_all2)

locations <- seppop(all_gi3a)
all_gi3b <- repool(locations$SEAFO_26,locations$SEAFO_24,locations$SEAFO_20,locations$Hotspot,locations$Rix,locations$Flats,locations$Three,locations$Johnies,
                  locations$South_Africa)

col_location <- c("azure3","darkgrey","azure4","olivedrab1","limegreen","green1","green2","palegreen4","darkgoldenrod2")

overall_stats3<-basic.stats(all_gi3b)
ci_fis3 <- boot.ppfis(all_gi3b)
overall_ar3 <- allelic.richness(all_gi3b)

allele_freq_loc <- pop.freq(all_gi3b,diploid=TRUE)
combined_allele_freq_loc2 <- do.call(rbind, lapply(allele_freq_loc, as.vector))
combined_allele_freq_loc3 <- as.data.frame(combined_allele_freq_loc2)

library(readr)
allele_freq_locations <- read_csv("allele_freq_locations.csv")

combined_allele_freq_loc4 <- cbind(allele_freq_locations$MAA_S26,allele_freq_locations$MAA_S24,
                                   allele_freq_locations$MAA_S20, allele_freq_locations$MAA_Hot,
                                   allele_freq_locations$MAA_Rix, allele_freq_locations$MAA_Three,
                                   allele_freq_locations$MAA_Flats,allele_freq_locations$MAA_Johnies,
                                   allele_freq_locations$MAA_SA)
                                                                
pheatmap(combined_allele_freq_loc4, cluster_rows=FALSE, cluster_cols=FALSE, 
         color = colorRampPalette(c("lightblue", "darkred"))(200),show_rownames = F,labels_col = c("SEAFO_26", "SEAFO_24", "SEAFo_20","Hotspot","Rix","Three","Flats","Johnies","South Africa"),
         main = "Heatmap of Allele Frequencies: Catch Location")

# prep the data for violin plots
overall_HO3 <- read_csv("overall_HO3.csv")
overall_HO3 <- cbind(overall_HO3,ho_all[,-2])
colnames(overall_HO3)[10] <- "SE_Atlantic"
overall_HO3_2 <- stack(overall_HO3)
colnames(overall_HO3_2)[1] <- "HO"
colnames(overall_HO3_2)[2] <- "Location"

overall_HE3 <- read_csv("overall_HE3.csv")
overall_HE3 <- cbind(overall_HE3,he_all[,-2])
colnames(overall_HE3)[10] <- "SE_Atlantic"
overall_HE3_2 <- stack(overall_HE3)
colnames(overall_HE3_2)[1] <- "HE"
colnames(overall_HE3_2)[2] <- "Location"

overall_FIS3 <- read_csv("overall_FIS3.csv")
overall_FIS3 <- cbind(overall_FIS3,fis_all[,-2])
colnames(overall_FIS3)[10] <- "SE_Atlantic"
overall_FIS3_2 <- stack(overall_FIS3)
colnames(overall_FIS3_2)[1] <- "FIS"
colnames(overall_FIS3_2)[2] <- "Location"

overall_AR3 <- read_csv("overall_AR3.csv")
overall_AR3 <- overall_AR3[,-1]
colnames(overall_AR3)[1] <- "SEAFO_26"
colnames(overall_AR3)[2] <- "SEAFO_24"
colnames(overall_AR3)[3] <- "SEAFO_20"
colnames(overall_AR3)[4] <- "Hotspot"
colnames(overall_AR3)[5] <- "Rix"
colnames(overall_AR3)[6] <- "Flats"
colnames(overall_AR3)[7] <- "Three"
colnames(overall_AR3)[8] <- "Johnies"
colnames(overall_AR3)[9] <- "South_Africa"

ar_all <- read.csv("overall_AR_new.csv", header = TRUE)
overall_AR3 <- cbind(overall_AR3,ar_all[,2])
colnames(overall_AR3)[10] <- "SE_Atlantic"
overall_AR3_2 <- stack(overall_AR3)
colnames(overall_AR3_2)[1] <- "AR"
colnames(overall_AR3_2)[2] <- "Location"

overall_stats_location <- cbind(overall_HO3_2,overall_HE3_2$HE,overall_FIS3_2$FIS,overall_AR3_2$AR)
colnames(overall_stats_location)[3] <- "HE"
colnames(overall_stats_location)[4] <- "FIS"
colnames(overall_stats_location)[5] <- "AR"

str(overall_stats_location)

overall_stats_location$Location <- factor(overall_stats_location$Location, levels = c("SEAFO_26", "SEAFO_24", "SEAFO_20","Hotspot", "Rix","Flats","Three","Johnies","South_Africa", "SE_Atlantic"))

measure2 <- gather(overall_stats_location, key = "measure", value = "values", c("HO", "HE", "FIS", "AR"))

facet.labels <- c("Observed Heterozygosity", "Expected Heterozygosity", "Inbreeding Coefficient", "Allelic Richness") #correct symbols
names(facet.labels) <- c("HO", "HE", "FIS", "AR")

facet2 <- ggplot(data = measure2, aes(x = Location, y = values, fill = Location)) +
  geom_violin(scale = "area", draw_quantiles = c(0.5))+
  facet_wrap(measure~., scales = "free_y", nrow = 3, labeller = labeller(measure = facet.labels))+
  #labs(x = "Locations")+
  scale_fill_manual(values=c("azure3","darkgrey","azure4","olivedrab1","limegreen","green1","green2","palegreen4","darkgoldenrod2","navyblue"))+
  theme_bw() + 
  theme(legend.position="none", 
        legend.key.size = unit(1, 'cm'), 
        legend.title = element_text(size=10), 
        legend.text = element_text(size=10), 
        strip.text.x = element_text(size = 10), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10, angle=90, hjust = 1))

#### calculate statistics ####
anova_HO3 <- aov(HO~Location, data = overall_stats_location)
summary(anova_HO3)
plot(anova_HO)

kw_HO3 <- kruskal.test(HO~Location, data = overall_stats_location)
wt_HO3<- pairwise.wilcox.test(overall_stats_location$HO, overall_stats_location$Location,
                     p.adjust.method = "BH") #BH is equivalent to false discovery rate

anova_HE3 <- aov(HE~Location, data = overall_stats_location)
summary(anova_HE)
plot(anova_HE)

kw_HE3 <- kruskal.test(HE~Location, data = overall_stats_location)
wt_HE3<- pairwise.wilcox.test(overall_stats_location$HE, overall_stats_location$Location,
                              p.adjust.method = "BH") #BH is equivalent to false discovery rate

anova_FIS3 <- aov(FIS~Location, data = overall_stats_location)
summary(anova_FIS)
plot(anova_HE)

kw_FIS3 <- kruskal.test(FIS~Location, data = overall_stats_location)
wt_FIS3<- pairwise.wilcox.test(overall_stats_location$FIS, overall_stats_location$Location,
                              p.adjust.method = "BH") #BH is equivalent to false discovery rate


# combine diversity plots into one plot
ggarrange(facet,facet2, align = "v", ncol = 1, labels = c("A)","B)"))


#### perform PCA ####
data_pca_sc3 <- scaleGen(all_gi3b, NA.method="mean",scale=F)
data_pca3 <- dudi.pca(data_pca_sc3, scale=F, nf= 10, scannf = F)

fviz_pca_ind(data_pca3,label="none",habillage = all_gi3b$pop,
             mean.point=F,pointsize=4,palette = col_location,
             addEllipses = FALSE, ellipse.level=0.5,ellipse.alpha = 0) + 
  theme_classic() + scale_shape_manual(values = c(17,17,17,19,19,19,19,19,15))



loadingplot(data_pca3, at=NULL, threshold=quantile(x,0.75),
            axis=1, fac=NULL, byfac=FALSE,
            lab=NULL, cex.lab=0.7, cex.fac=1, lab.jitter=0,
            main="Loading plot", xlab="Variables", ylab="Loadings",
            srt = 0, adj = NULL, ...)

#### perform DAPC ####
dapc2 <- dapc(all_gi3b,var.contrib = TRUE, scale = FALSE, n.pca = 85, n.da=NULL) #data is a genind object  

# number of discriminant factors: n.da = nPop(dedup_gi2) - 1

scatter(dapc2, scree.da = FALSE, bg = "white", 
        posi.pca = "none", legend = TRUE, 
        cellipse = 1, cstar = 0,
        col = c(col_location),
        pch = 16,
        txt.leg = paste(levels(all_gi3b$pop)), 
        clabel = 1, cex = 2)
txt.leg = paste(levels(data2$pop))

# cross-validate dapc
dedup_gi2_xval <- xvalDapc(tab(all_gi3b, NA.method = "mean"), pop(all_gi3b))
dedup_gi2_xval2 <- xvalDapc(tab(all_gi3b, NA.method = "mean"), pop(all_gi3b),
                            #n.pca = 18:52, 
                            n.rep = 100,
                            parallel = "multicore", ncpus = 4L)

compoplot(dapc2, 
          #subset = 1:65, 
          col.pal = col_location,
          border = NA,
          #col = c("blue","olivegreen", "grey"),
          #lab = "", border = NA,
          legend = F,
          show.lab = F,
          posi = "top")
dev.off()

# kmeans clustering

nb_k_all <- find.clusters(all_gi3b)
kmeans_all_dapc <- dapc(all_gi3b, nb_k_all$grp)

scatter(kmeans_all_dapc, scree.da = TRUE, bg = "white", 
        posi.pca = "bottomleft", legend = TRUE, 
        cellipse = 1, cstar = 0,
        col = c(col_location),
        pch = 16,
        txt.leg = paste(levels(all_gi3b)), 
        clabel = 1, cex = 2)

compoplot(kmeans_all_dapc, 
          #subset = 1:65, 
          col.pal = col_location,
          border = NA,
          #col = c("blue","olivegreen", "grey"),
          #lab = "", border = NA,
          legend = F,
          show.lab = F,
          posi = "top")


#### calculate FST ####

all_gl3b <- gi2gl(all_gi3b)
location_fst <- stamppFst(all_gl3b, nboots = 10000, percent = 95)

write.table(location_fst$Fsts, "location_fst.txt")
write.table(location_fst$Pvalues, "location_pvalues.txt")

# adjusted p-value using BF

# first, turn make pvalues into a vector

pvalues <- c(0,0,0,0,0,0,0,0,0,0,0.7878,0,0,0.2006,0,0,0,0,0,0,0,0,0,0,0,0,0,0.001,0.5094,0,1.00E-04,0,0,0,0)
fdrs<-p.adjust(pvalues, method="BH")

print(fdrs)

#### Calculate Admixture ####

library(LEA)
gl2geno(all_gl3b, outfile = "gl_geno_primary2", outpath = "C:/Users/u05086605/Desktop/PROJECTS/SEAFO_orange_roughy_2022", verbose = NULL)

lea_pca <- pca("gl_geno_primary2.lfmm", scale=FALSE)
tw <- tracy.widom(lea_pca)
plot(tw$percentage)

project1=NULL
project1 = snmf("gl_geno_primary2.geno",
               K = 1:10,
               entropy = TRUE,
               repetitions = 100,project = "new")

# check cross-entropy
plot(project1, col = "blue", pch = 19, cex = 1.2)

# plot ancestry proportions
best = which.min(cross.entropy(project1, K = 2))
my.colors <- c("palegreen4","darkgoldenrod2","blue")

barplot(t(Q(project1, K = 2, run = best)), 
        col = my.colors,
        xlab = "Individuals",
        ylab = "Ancestry proportions")

barchart(project1, K = 2, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp

axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)

#### Perform migration analyses ####

# BayesAss was run in a high performance cluster, results were imported into R for visualization

# BayesAss results
install.packages("circlize")
install.packages("migest")

library(circlize)
library(migest)
library(dplyr)

bayesass <- tibble(Location = c("SEAFO_26","SEAFO_24","SEAFO_20","Hotspot","Rix","Flats","Three","Johnies","South_Africa"),
                   SEAFO_26=c(0.8386,0.0174,0.0172,0.0178,0.0181,0.0169,0.017,0.0172,0.0398),
                   SEAFO_24=c(0.0178,0.8522,0.0174,0.0174,0.0172,0.0176,0.0172,0.0175,0.0258),
                   SEAFO_20=c(0.0175,0.0175,0.8592,0.0175,0.0178,0.0179,0.0172,0.0177,0.0177),
                   Hotspot=c(0.011,0.0115,0.0115,0.9095,0.0115,0.0111,0.0113,0.0111,0.0114),
                   Rix=c(0.0128,0.0131,0.0123,0.0126,0.8972,0.0137,0.013,0.0127,0.0126),
                   Flats=c(0.019,0.0186,0.0191,0.0185,0.0186,0.8385,0.0212,0.0287,0.0178),
                   Three=c(0.0175,0.0174,0.0168,0.0179,0.0168,0.0358,0.8422,0.018,0.0176),
                   Johnies=c(0.0118,0.0119,0.0117,0.0117,0.0124,0.0121,0.012,0.9047,0.0118),
                   South_Africa=c(0.012,0.0119,0.0118,0.0119,0.0119,0.0113,0.0121,0.0123,0.9048))


bayesass.mat = as.matrix(bayesass[, 2:10])
bayesass.mat

dimnames(bayesass.mat) = list(source = bayesass$Location, sink = bayesass$Location)
bayesass.mat


circos.par(start.degree = 90, gap.degree = 6)
chordDiagram(x = bayesass.mat, grid.col = col_location, grid.border = "black", transparency = 0.25,
             order = bayesass$Location, directional = F, direction.type = "arrows",
             self.link = 1, preAllocateTracks = list(track.height = 0.1),
             annotationTrack = "grid", annotationTrackHeight = c(0.1, 0.1),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.15, link.arr.lty = 3, link.arr.col = "#252525", 
             link.largest.ontop = F)


circos.trackPlotRegion(track.index = 1,
                       bg.border = NA,
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         sector.index = get.cell.meta.data("sector.index")
                         # Text direction
                         theta = circlize(mean(xlim), 1)[1, 1] %% 360
                         dd = ifelse(theta < 180 || theta > 360, "bending.inside", "bending.outside")
                         circos.text(x = mean(xlim), y = 0.8, labels = sector.index, facing = dd,
                                     niceFacing = TRUE, cex = 1.5, font = 1)
                       })



#### 3. Isolation-by-Distance tests ####

#### calculate distance along bathymetry using marmap #### 

library(marmap)

bathymetry2 <- getNOAA.bathy(lon1 = 5, lon2 = 22, lat1 =-19, lat2 = -37, resolution = 10)
trans1 <- trans.mat(bathymetry2)
trans2 <- trans.mat(bathymetry2, min.depth=300, max.depth =-2000)

lon <- c(5.8,6.47933333,9.04766667,10.1242,12.7315,13.3832,13.327,13.5988,21.584167)
lat <- c(-26.0771667,-24.6336667,-20.5755,-19.3753,-22.4928,-24.537,-24.6772,-26.36,-36.668333)

sites <- as.data.frame(cbind(lon, lat))
rownames(sites)[1] <- "SEAFO_26"
rownames(sites)[2] <- "SEAFO_24"
rownames(sites)[3] <- "SEAFO_20"
rownames(sites)[4] <- "Hotspot"
rownames(sites)[5] <- "Rix"
rownames(sites)[6] <- "Flats"
rownames(sites)[7] <- "Three"
rownames(sites)[8] <- "Johnies"
rownames(sites)[9] <- "South_Africa"

out1 <- lc.dist(trans1, sites, res = "path") #--> least cost path avoiding landmasses
out2 <- lc.dist(trans2, sites, res = "path") #--> least cost path for a maximum depth of 2000m

dist1 <- lc.dist(trans1, sites, res = "dist")
dist2 <- lc.dist(trans2, sites, res = "dist") 

write_csv(as.data.frame(dist1), "distances_bathymetry.csv")
write.table(as.matrix(dist2), "distances_bathymetry2.txt")

str(dist1)
plot(bathymetry2, 
     deep = c(-5000,-2000, -1000, -500), shallow = c(-300, -150, 0), col = c("darkblue", "blue", "turquoise", "lightblue"), 
     step = c(500, 100,100, 1),lty = c(1,1, 1, 1), lwd = c(0.6,0.6, 0.6, 1.2), draw = c(FALSE,FALSE, FALSE, FALSE)) 
points(sites, pch = 21, col = "orange", bg = col2alpha("orange", .9), cex = 1.2) 
text(sites[,1], sites[,2], lab = rownames(sites), pos = c(3, 4, 1, 2), col = "blue") 
lapply(out1, lines, col = "orange", lwd = 5, lty = 1)-> dummy 
lapply(out2, lines, col = "black", lwd = 5, lty = 1)-> dummy

d <- dist2isobath(bathymetry2, lon, lat, isobath = -900)

plot(bathymetry2, image = TRUE, lwd = 0.1, land = TRUE, bpal = list(c(0, max(bathymetry2), "lightgrey"), c(min(bathymetry2), 0, "lightblue")))
plot(bathymetry2, deep = 0, shallow = 0, step = 0, lwd = 0.6, add = TRUE)
points(lon, lat, pch = 21, bg = "orange2", cex = 0.8)


#### calculate isolation by distance using linear models ####

library(readxl)
distance_model1 <- read_excel("distances_3models.xlsx", 
                              sheet = "model1_general_LCP")
distance1_vector <- distance_model1[upper.tri(distance_model1)]

distance_model2 <- read_excel("distances_3models.xlsx", 
                              sheet = "model2_bathymetry_2000")
distance2_vector <- distance_model2[upper.tri(distance_model2)]


distance_model3<- read_excel("distances_3models.xlsx", 
                             sheet = "model3_Stepping_Stone")
distance3_vector <- distance_model3[upper.tri(distance_model3)]

genetic_distance <- read_excel("pairwise_fst_locations.xlsx",sheet = "FST_fullmatrix")
genetic_distance_vector <- genetic_distance[upper.tri(genetic_distance)]


# fit the model
model1 <- lm(genetic_distance_vector ~ distance1_vector)
summary(model1)

model2 <- lm(genetic_distance_vector ~ distance2_vector)
summary(model2)

model3 <- lm(genetic_distance_vector ~ distance3_vector)
summary(model3)
plot(model3)

# test which model is the most likely
AIC(model1)
AIC(model2)
AIC(model3)

BIC(model1)
BIC(model2)
BIC(model3)

#calculate correlation values

corr1 <- cor.test(genetic_distance_vector,distance1_vector,method="pearson")
corr2 <- cor.test(genetic_distance_vector,distance2_vector,method="spearman")
corr3 <- cor.test(genetic_distance_vector,distance3_vector,method="spearman")

genetic_distance2 <- read_excel("pairwise_fst_locations.xlsx",sheet = "FST(1-FST)_fullmatrix")
genetic_distance_vector2 <- genetic_distance[upper.tri(genetic_distance)]

corr1_2 <- cor.test(genetic_distance_vector2,distance1_vector,method="spearman")
corr2_2 <- cor.test(genetic_distance_vector2,distance2_vector,method="spearman")
corr3_2 <- cor.test(genetic_distance_vector2,distance3_vector,method="spearman")

# Create a data frame for plotting
plot_data <- data.frame(
  Genetic_Distance = genetic_distance_vector,
  Geographic_Distance = distance1_vector
)

# Create scatter plot
ggplot(plot_data, aes(x = Geographic_Distance, y = Genetic_Distance)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(x = "Geographic Distance", y = "Genetic Distance", title = "Isolation by Distance") +
  theme_minimal()


#### calculate IBD using mantel tests ####
library(vegan)

dist_model1 <- read_excel("distances_3models.xlsx", 
                          sheet = "model1_halfmatrix")

dist_model2 <- read_excel("distances_3models.xlsx", 
                          sheet = "model2_halfmatrix")

dist_model3<- read_excel("distances_3models.xlsx", 
                         sheet = "model3_halfmatrix")
gene_distance <- read_excel("pairwise_fst_locations.xlsx",sheet = "FST_halfmatrix")


mantel1 <- mantel(dist_model1,gene_distance, permutations=9999)
mantel2 <- mantel(dist_model2,gene_distance)
mantel3 <- mantel(dist_model3,gene_distance)

gene_distance_2  <- read_excel("pairwise_fst_locations.xlsx",sheet = "FST(1-FST)")

mantel1_2 <- mantel(dist_model1,gene_distance_2)
mantel2_2 <- mantel(dist_model2,gene_distance_2)
mantel3_2  <- mantel(dist_model3,gene_distance_2)


#### 4. Re-do analyses without SA ####

vcf2 <- file.choose(new=FALSE) 
data_subset_vcf <- read.vcfR(vcf2,verbose=TRUE)
data_subset <- vcfR2genlight(data_subset_vcf)

data_subset1 <- gl.compliance.check(data_subset)

# recode individuals names
names_ind <- as.data.frame(data_subset1$ind.names)
ind.recode <- read.csv("samples_names_withouSA.csv", header=FALSE)
ind.recode2 <- cbind(names_ind,ind.recode)
write_csv(ind.recode2, "ind.recode_withouSA.csv") #don't forget to open it and delete the first row

data_subset2 <- gl.recode.ind(data_subset1, "ind.recode_withouSA.csv", recalc = TRUE, mono.rm = TRUE, verbose = NULL)


# import metadata
metadata2 <- read_excel("metadata1.xlsx", sheet="popmap_2")

# overall
data_subset2$pop <- as.factor(metadata2$Overall)
as.data.frame(data_subset2$ind.names)
unique(as.data.frame(data_subset2$pop))
subset_gi <- gl2gi(data_subset2)

subset_pca_sc <- scaleGen(subset_gi, NA.method="mean",scale=F)
subset_pca <- dudi.pca(subset_pca_sc, scale=F, nf= 10, scannf = F)

fviz_pca_ind(subset_pca,label="none",habillage = subset_gi$pop,axes=c(1,2),
             legend.title ="Region",mean.point=F,pointsize=4,
             addEllipses = FALSE, ellipse.level=0.5,ellipse.alpha = 0, title="PCA") + theme_classic()

# region
data_subset2$pop <- as.factor(metadata2$Region)
as.data.frame(data_subset2$ind.names)
unique(as.data.frame(data_subset2$pop))
subset_gi2 <- gl2gi(data_subset2)

subset_pca_sc2 <- scaleGen(subset_gi2, NA.method="mean",scale=F)
subset_pca2 <- dudi.pca(subset_pca_sc2, scale=F, nf= 10, scannf = F)

col_region2 <- c("limegreen","azure3")
fviz_pca_ind(subset_pca2,label="none",habillage = subset_gi2$pop,axes=c(1,2),
             legend.title ="Region",mean.point=F,pointsize=4,palette = col_region2,
             addEllipses = FALSE, ellipse.level=0.5,ellipse.alpha = 0, title="PCA per region") + theme_classic()

# location
data_subset2$pop <- as.factor(metadata2$Location)
as.data.frame(data_subset2$ind.names)
unique(as.data.frame(data_subset2$pop))
subset_gi3 <- gl2gi(data_subset2)

locations2 <- seppop(subset_gi3)
subset_gi3b <- repool(locations2$SEAFO_26,locations2$SEAFO_24,locations2$SEAFO_20,locations2$Hotspot,locations2$Rix,locations2$Flats,locations2$Three,locations2$Johnies)

subset_pca_sc3 <- scaleGen(subset_gi3b, NA.method="mean",scale=F)
subset_pca3 <- dudi.pca(subset_pca_sc3, scale=F, nf= 10, scannf = F)

col_location2 <- c("azure3","darkgrey","azure4","olivedrab1","limegreen","green1","green2","palegreen4")

fviz_pca_ind(subset_pca3,label="none",habillage = subset_gi3b$pop,axes=c(1,2),
             mean.point=F,pointsize=4,palette = col_location2,
             addEllipses = FALSE, ellipse.level=0.5,ellipse.alpha = 0, title="Orange roughy") + 
  theme_classic() + scale_shape_manual(values = c(17,17,17,19,19,19,19,19))

grp <- find.clusters(subset_gi3b, max.n.clust=10) 
table(pop(subset_gi3b), grp$grp)
dapc3 <- dapc(subset_gi3b, grp$grp) #data is a genind object  
scatter(dapc3, scree.da = FALSE, bg = "white", 
        posi.pca = "bottomleft", legend = FALSE, 
        cellipse = 1, cstar = 0,
        col = c(col),
        pch = 16,
        txt.leg = paste(levels(subset_gi3$pop)), 
        clabel = 1, cex = 2)

compoplot(dapc3, 
          #subset = 1:65, 
          col.pal = col,
          border = NA,
          #col = c("blue","olivegreen", "grey"),
          #lab = "", border = NA,
          legend = F,
          show.lab = TRUE,
          posi = "topright")

subset_gl3b <- gi2gl(subset_gi3b)
subset_fst <- stamppFst(subset_gl3b, nboots = 10000, percent = 95)

write.table(subset_fst$Fsts, "subset_location_fst.txt")
write.table(subset_fst$Pvalues, "subset_location_pvalues.txt")

# IBD

distance_model1_1 <- read_excel("distances_3models.xlsx", 
                              sheet = "model1_general_LCP_noSA")
distance1_vector1_1 <- distance_model1_1[upper.tri(distance_model1_1)]

distance_model2_2 <- read_excel("distances_3models.xlsx", 
                              sheet = "model2_bathymetry_2000_noSA")
distance2_vector2_2 <- distance_model2_2[upper.tri(distance_model2_2)]

distance_model3_3<- read_excel("distances_3models.xlsx", 
                             sheet = "model3_Stepping_Stone_noSA")
distance3_vector3_3 <- distance_model3_3[upper.tri(distance_model3_3)]

genetic_distance2 <- read_excel("pairwise_fst_locations.xlsx",sheet = "FST_fullmatrix_noSA")
genetic_distance_vector2 <- genetic_distance2[upper.tri(genetic_distance2)]


# fit the model
model1_1 <- lm(genetic_distance_vector2 ~ distance1_vector1_1)
summary(model1_1)

model2_2 <- lm(genetic_distance_vector2 ~ distance2_vector2_2)
summary(model2_2)

model3_3 <- lm(genetic_distance_vector2 ~ distance3_vector3_3)
summary(model3_3)
plot(model3)

# test which model is the most likely
AIC(model1)
AIC(model2)
AIC(model3)

BIC(model1)
BIC(model2)
BIC(model3)

#calculate correlation values

corr1_1 <- cor.test(genetic_distance_vector2,distance1_vector1_1,method="spearman")
corr2_2 <- cor.test(genetic_distance_vector2,distance2_vector2_2,method="spearman")
corr3_3 <- cor.test(genetic_distance_vector2,distance3_vector3_3,method="spearman")

genetic_distance2 <- read_excel("pairwise_fst_locations.xlsx",sheet = "FST(1-FST)_fullmatrix")
genetic_distance_vector2 <- genetic_distance[upper.tri(genetic_distance)]

corr1_2 <- cor.test(genetic_distance_vector2,distance1_vector,method="spearman")
corr2_2 <- cor.test(genetic_distance_vector2,distance2_vector,method="spearman")
corr3_2 <- cor.test(genetic_distance_vector2,distance3_vector,method="spearman")

# Create a data frame for plotting
plot_data <- data.frame(
  Genetic_Distance = genetic_distance_vector2,
  Geographic_Distance = distance1_vector1_1
)

# Create scatter plot
ggplot(plot_data, aes(x = Geographic_Distance, y = Genetic_Distance)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(x = "Geographic Distance", y = "Genetic Distance", title = "Isolation by Distance") +
  theme_minimal()

# conduct a mantel test
library(vegan)

dist_model1_1 <- read_excel("distances_3models.xlsx", 
                          sheet = "model1_halfmatrix_noSA")

dist_model2_2 <- read_excel("distances_3models.xlsx", 
                          sheet = "model2_halfmatrix_noSA")

dist_model3_3<- read_excel("distances_3models.xlsx", 
                         sheet = "model3_halfmatrix_noSA")
gene_distance2 <- read_excel("pairwise_fst_locations.xlsx",sheet = "FST_halfmatrix_noSA")

mantel1_1 <- mantel(dist_model1_1,gene_distance2, permutations=9999)
mantel2_2 <- mantel(dist_model2_2,gene_distance2)
mantel3_3 <- mantel(dist_model3_3,gene_distance2)

gene_distance_2  <- read_excel("pairwise_fst_locations.xlsx",sheet = "FST(1-FST)")

mantel1_2 <- mantel(dist_model1,gene_distance_2)
mantel2_2 <- mantel(dist_model2,gene_distance_2)
mantel3_2  <- mantel(dist_model3,gene_distance_2)


#### 5. Outlier Analyses ####

#### bayescan ####
# bayescan --> runs in the CHPC
gl2bayescan(all_gl3b, outfile = "all_sites_bayescan_primary2.txt", outpath = "C:/Users/u05086605/Desktop/PROJECTS/SEAFO_orange_roughy_2022", verbose = NULL)

source("plot_R.r")

bayescan <- read.table("all_sites_bayescan_primary2_fst.txt")
loci <- as.data.frame(data_all2$loc.names)
write_csv(loci, "SNPs_names_primary2_2.csv")

bayescan <- cbind(loci,bayescan)
colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST") 

attach(bayescan)
class(bayescan$Q_VALUE)  
bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE) 
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 

bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4)) 
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4)) 
bayescan$ALPHA <- (round(bayescan$ALPHA, 4)) 
bayescan$FST <- (round(bayescan$FST, 6))

bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.05,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing")) 
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION) 

positive <- bayescan[bayescan$SELECTION=="diversifying",] #--> there are no loci under positive selection
neutral <- bayescan[bayescan$SELECTION=="neutral",] 
balancing <- bayescan[bayescan$SELECTION=="balancing",]  

xtabs(data=bayescan, ~SELECTION) 

range(bayescan$Q_VALUE) 
bayescan$LOG10_Q <- -log10(bayescan$Q_VALUE) 

# there are no loci under positive selection, only balancing

#positive <- bayescan[bayescan$SELECTION=="diversifying",] 
neutral <- bayescan[bayescan$SELECTION=="neutral",] 
balancing <- bayescan[bayescan$SELECTION=="balancing",]  

xtabs(data=bayescan, ~SELECTION) 

range(bayescan$Q_VALUE) 
bayescan$LOG10_Q <- -log10(bayescan$Q_VALUE) 

x_title="Log(q-value)" 
y_title="Fst" 

graph_1<-ggplot(bayescan,aes(x=LOG10_Q,y=FST, label=SELECTION)) 
graph_1+geom_point(aes(fill=SELECTION), pch=21, size=2)+  
  scale_fill_manual(name="Selection",values=c("white","grey"))+ 
  labs(x=x_title)+ 
  labs(y=y_title)+ 
  theme(axis.title=element_text(size=12, family="Helvetica",face="bold"), legend.position="none")+ 
  theme(axis.text.x=element_text(colour="black"))+ 
  theme(axis.text.y=element_text(colour="black",size=12))+ 
  theme(axis.text.x=element_text(colour="black",size=12))+ 
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),  
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold")) +
  theme_classic()


plot_bayescan("all_sites_bayescan_fst.txt")



#### basic top FST loci for outlier detection ####
library(readr)
regions_fst <- as.data.frame(read_delim("fst_regions_62364loci.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE))
View(regions_fst)
regions_fst <- cbind(regions_fst,loci)

regions_fst$seafo_nam <- as.numeric(as.character(regions_fst$seafo_nam))
regions_fst$seafo_sa <- as.numeric(unlist(regions_fst$seafo_sa))
regions_fst$nam_sa <- as.numeric(unlist(regions_fst$nam_sa))

mean(regions_fst$seafo_nam) # 0.0004719614
mean(regions_fst$seafo_sa,na.rm=TRUE) # 0.001443721 --> had quite a few missing data
mean(regions_fst$nam_sa,na.rm=TRUE) # 0.001493804

quantile(regions_fst$seafo_nam,0.95) # 0.03411893
quantile(regions_fst$seafo_nam,0.99) # 0.0657623 

quantile(regions_fst$seafo_sa,na.rm=TRUE,0.95) # 0.0616312
quantile(regions_fst$seafo_sa,na.rm=TRUE,0.99) # 0.1126109 

quantile(regions_fst$nam_sa,na.rm=TRUE,0.95) # 0.0489152
quantile(regions_fst$nam_sa,na.rm=TRUE,0.99) # 0.0967133 

outliers_region_fst99 <- regions_fst[regions_fst$seafo_nam > 0.064 & regions_fst$seafo_sa > 0.113 & regions_fst$nam_sa > 0.096,]
outliers_seafo_nam_fst99 <- regions_fst[regions_fst$seafo_nam > 0.064,6]
outliers_seafo_sa_fst99 <- regions_fst[regions_fst$seafo_sa > 0.112,6]
outliers_nam_sa_fst99 <- regions_fst[regions_fst$nam_sa > 0.095,6]

# Example: create a factor column for color groups before plotting

seafo_nam <- ggplot(regions_fst, aes(x = loci, y = seafo_nam)) +
  geom_point(aes(colour = factor(case_when(
    seafo_nam > 0.064 ~ "High",
    seafo_nam > 0.033 & seafo_nam <= 0.064 ~ "Medium",
    TRUE ~ "Low"
  ))), show.legend = FALSE) + 
  scale_color_manual(values = c("Low" = "lightgrey", "Medium" = "lightpink", "High" = "maroon")) +
  geom_hline(yintercept = 0.033, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 0.45)) +
  theme_classic() + theme(axis.text.x = element_blank()) + scale_x_discrete() +
  labs(x = NULL, y = "SEAFO vs Namibia") 

seafo_sa <- ggplot(regions_fst, aes(x = loci, y = seafo_sa)) +
  geom_point(aes(colour = factor(case_when(
    seafo_sa > 0.113 ~ "High",
    seafo_sa > 0.061 & seafo_nam <= 0.113 ~ "Medium",
    TRUE ~ "Low"
  ))), show.legend = FALSE) + 
  scale_color_manual(values = c("Low" = "lightgrey", "Medium" = "lightpink", "High" = "maroon")) +
  geom_hline(yintercept = 0.0616, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 0.45)) +
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(x = NULL, y = "SEAFO vs South Africa") 

nam_sa <- ggplot(regions_fst, aes(x = loci, y = nam_sa)) +
  geom_point(aes(colour = factor(case_when(
    nam_sa > 0.096 ~ "High",
    nam_sa > 0.048 & seafo_nam <= 0.096 ~ "Medium",
    TRUE ~ "Low"
  ))), show.legend = FALSE) + 
  scale_color_manual(values = c("Low" = "lightgrey", "Medium" = "lightpink", "High" = "maroon")) +
  geom_hline(yintercept = 0.048, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 0.45))  +
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(x = NULL, y = "Namibia vs South Africa") 

region_fst_fig <- grid.arrange(seafo_nam,seafo_sa,nam_sa,ncol=1,nrow=3)

locations_fst <- as.data.frame(read_delim("locations_fst.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE))
View(locations_fst)
locations_fst$loci <- paste(locations_fst$CHROM,"_",locations_fst$POS)

#99%
mean(locations_fst$seafo26_20,na.rm=TRUE) # 0.222
mean(locations_fst$seafo26_hotspot,na.rm=TRUE) # 0.183673  
mean(locations_fst$seafo26_johnies,na.rm=TRUE) # 0.1974
mean(locations_fst$seafo26_sa,na.rm=TRUE) # 0.185953

mean(locations_fst$seafo20_hotspot,na.rm=TRUE) # 0.187839 
mean(locations_fst$seafo20_johnies,na.rm=TRUE) # 0.193869
mean(locations_fst$seafo20_sa,na.rm=TRUE) # 0.193863 

mean(locations_fst$hotspot_johnies,na.rm=TRUE) # 0.127799
mean(locations_fst$hotspot_sa,na.rm=TRUE) # 0.138007 

mean(locations_fst$johnies_sa,na.rm=TRUE) # 0.137286

#95%
quantile(locations_fst$seafo26_20,0.95,na.rm=TRUE) # 0.122807 
quantile(locations_fst$seafo26_hotspot,0.95,na.rm=TRUE) # 0.0997809 
quantile(locations_fst$seafo26_johnies,0.95,na.rm=TRUE) # 0.104588 
quantile(locations_fst$seafo26_sa,0.95,na.rm=TRUE) # 0.104588

quantile(locations_fst$seafo20_hotspot,0.95,na.rm=TRUE) # 0.0997809 
quantile(locations_fst$seafo20_johnies,0.95,na.rm=TRUE) # 0.105
quantile(locations_fst$seafo20_sa,0.95,na.rm=TRUE) # 0.107706

quantile(locations_fst$hotspot_johnies,0.95,na.rm=TRUE) # 0.0729592 
quantile(locations_fst$hotspot_sa,0.95,na.rm=TRUE) # 0.076901

quantile(locations_fst$johnies_sa,0.95,na.rm=TRUE) # 0.0774278

seafo26_20 <- ggplot(locations_fst, aes(x = loci, y = seafo26_20)) +
  geom_point(aes(colour = factor(case_when(
    seafo26_20 > 0.22 ~ "High",
    seafo26_20 > 0.12 & seafo26_20 <= 0.22 ~ "Medium",
    TRUE ~ "Low"
  ))), show.legend = FALSE) + 
  scale_color_manual(values = c("Low" = "lightgrey", "Medium" = "lightpink", "High" = "maroon")) +
  geom_hline(yintercept = 0.12, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(x = NULL, y = "SEAFO 26 vs SEAFO 20") 

seafo26_hotspot <- ggplot(locations_fst, aes(x = loci, y = seafo26_hotspot)) +
  geom_point(aes(colour = factor(case_when(
    seafo26_hotspot > 0.18 ~ "High",
    seafo26_hotspot > 0.09 & seafo26_hotspot <= 0.18 ~ "Medium",
    TRUE ~ "Low"
  ))), show.legend = FALSE) + 
  scale_color_manual(values = c("Low" = "lightgrey", "Medium" = "lightpink", "High" = "maroon")) +
  geom_hline(yintercept = 0.09, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(x = NULL, y = "SEAFO 26 vs Hotspot") 

seafo26_johnies <- ggplot(locations_fst, aes(x = loci, y = seafo26_johnies)) +
  geom_point(aes(colour = factor(case_when(
    seafo26_johnies > 0.19 ~ "High",
    seafo26_johnies > 0.10 & seafo26_johnies <= 0.19 ~ "Medium",
    TRUE ~ "Low"
  ))), show.legend = FALSE) + 
  scale_color_manual(values = c("Low" = "lightgrey", "Medium" = "lightpink", "High" = "maroon")) +
  geom_hline(yintercept = 0.10, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(x = NULL, y = "SEAFO 26 vs Johnies")

seafo26_sa <- ggplot(locations_fst, aes(x = loci, y = seafo26_sa)) +
  geom_point(aes(colour = factor(case_when(
    seafo26_sa > 0.18 ~ "High",
    seafo26_sa > 0.10 & seafo26_sa <= 0.18 ~ "Medium",
    TRUE ~ "Low"
  ))), show.legend = FALSE) + 
  scale_color_manual(values = c("Low" = "lightgrey", "Medium" = "lightpink", "High" = "maroon")) +
  geom_hline(yintercept = 0.10, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(x = NULL, y = "SEAFO 26 vs South Africa")

grid.arrange(seafo26_20,seafo26_hotspot,seafo26_johnies,seafo26_sa)

seafo20_hotspot <- ggplot(locations_fst, aes(x = loci, y = seafo20_hotspot)) +
  geom_point(aes(colour = factor(case_when(
    seafo20_hotspot > 0.18 ~ "High",
    seafo20_hotspot > 0.10 & seafo20_hotspot <= 0.18 ~ "Medium",
    TRUE ~ "Low"
  ))), show.legend = FALSE) + 
  scale_color_manual(values = c("Low" = "lightgrey", "Medium" = "lightpink", "High" = "maroon")) +
  geom_hline(yintercept = 0.10, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(x = NULL, y = "SEAFO 20 vs Hotspot") 

seafo20_johnies <- ggplot(locations_fst, aes(x = loci, y = seafo20_johnies)) +
  geom_point(aes(colour = factor(case_when(
    seafo20_johnies > 0.19 ~ "High",
    seafo20_johnies > 0.10 & seafo20_johnies <= 0.19 ~ "Medium",
    TRUE ~ "Low"
  ))), show.legend = FALSE) + 
  scale_color_manual(values = c("Low" = "lightgrey", "Medium" = "lightpink", "High" = "maroon")) +
  geom_hline(yintercept = 0.10, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(x = NULL, y = "SEAFO 20 vs Johnies")

seafo20_sa <- ggplot(locations_fst, aes(x = loci, y = seafo20_sa)) +
  geom_point(aes(colour = factor(case_when(
    seafo20_sa > 0.19 ~ "High",
    seafo20_sa > 0.10 & seafo20_sa <= 0.19 ~ "Medium",
    TRUE ~ "Low"
  ))), show.legend = FALSE) + 
  scale_color_manual(values = c("Low" = "lightgrey", "Medium" = "lightpink", "High" = "maroon")) +
  geom_hline(yintercept = 0.10, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 0.6)) +
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(x = NULL, y = "SEAFO 20 vs South Africa")

grid.arrange(seafo20_hotspot,seafo20_johnies,seafo20_sa)

hotspot_johnies <- ggplot(locations_fst, aes(x = loci, y = hotspot_johnies)) +
  geom_point(aes(colour = factor(case_when(
    hotspot_johnies > 0.19 ~ "High",
    hotspot_johnies > 0.10 & hotspot_johnies <= 0.19 ~ "Medium",
    TRUE ~ "Low"
  ))), show.legend = FALSE) + 
  scale_color_manual(values = c("Low" = "lightgrey", "Medium" = "lightpink", "High" = "maroon")) +
  geom_hline(yintercept = 0.10, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(x = NULL, y = " Hotspot vs Johnies")

hotspot_sa <- ggplot(locations_fst, aes(x = loci, y = hotspot_sa)) +
  geom_point(aes(colour = factor(case_when(
    hotspot_sa > 0.19 ~ "High",
    hotspot_sa > 0.10 & hotspot_sa <= 0.19 ~ "Medium",
    TRUE ~ "Low"
  ))), show.legend = FALSE) + 
  scale_color_manual(values = c("Low" = "lightgrey", "Medium" = "lightpink", "High" = "maroon")) +
  geom_hline(yintercept = 0.10, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(x = NULL, y = "Hotspot vs South Africa")

johnies_sa <- ggplot(locations_fst, aes(x = loci, y = johnies_sa)) +
  geom_point(aes(colour = factor(case_when(
    johnies_sa > 0.19 ~ "High",
    johnies_sa > 0.10 & johnies_sa <= 0.19 ~ "Medium",
    TRUE ~ "Low"
  ))), show.legend = FALSE) + 
  scale_color_manual(values = c("Low" = "lightgrey", "Medium" = "lightpink", "High" = "maroon")) +
  geom_hline(yintercept = 0.10, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_classic() + theme(axis.text.x = element_blank()) + 
  labs(x = NULL, y = " Johnies vs South Africa")

grid.arrange(hotspot_johnies,hotspot_sa,johnies_sa)

outliers_region_fst99 <- regions_fst[regions_fst$seafo_nam > 0.064 & regions_fst$seafo_sa > 0.113 & regions_fst$nam_sa > 0.096,]
outliers_seafo_nam_fst99 <- regions_fst[regions_fst$seafo_nam > 0.064,11]
outliers_seafo_sa_fst99 <- regions_fst[regions_fst$seafo_sa > 0.113,11]
outliers_nam_sa_fst99 <- regions_fst[regions_fst$nam_sa > 0.096,11]

#### use pcadapt for outlier detection ####
library(pcadapt)

all_snps <- read.pcadapt(vcf, type = "vcf")
k_nb <- pcadapt(input = all_snps, K = 20) 
plot(k_nb, option = "screeplot") # We recommend to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule)

# in this case, k_nb = 7  
# also going to test k_nb =5 

# compute z-scores based on PCA
test_stat <- pcadapt(all_snps, K = 5)
summary(test_stat)

# plot p-values 
plot(test_stat , option = "manhattan")
plot(test_stat, option = "qqplot")

hist(test_stat$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

plot(test_stat, option = "stat.distribution")

# choose the threshold to detect outliers

#qvalues
BiocManager::install("qvalue")
library(qvalue)

# bind loci names with the outputs
snps <- cbind(loci,test_stat$pvalues)
colnames(snps)=c("SNP","pvalues") 

qval <- qvalue(snps$pvalues)$qvalues
qval <- cbind(loci,qval)
alpha <- 0.05
outliers <- subset(qval, qval<alpha)
length(outliers)
write_csv(outliers, "pcadapt_outliers_k5.csv")

# BH
padj <- p.adjust(test_stat$pvalues,method="BH")
padj <- cbind(loci,padj)
outliers2 <- subset(padj, padj < alpha)
length(outliers2)

# associate PCs with outliers
snp_pc <- get.pc(test_stat, outliers)


#### outlier detection wiht outflank ####
outflank1 <- gl.outflank(data_all2)
View(outflank1$outflank$results)

outflank1$outflank$dfInferred
#[1] 8.681712 --> need to report it

outflank1$outflank$FSTbar
#[1] 0.001573408

outflank1$outflank$FSTNoCorrbar
#[1] 0.03705964

# there was a good fit of the model to the data, but we need to check if all assumptions were met

# check the quality of the results (http://rstudio-pubs-static.s3.amazonaws.com/305384_9aee1c1046394fb9bd8e449453d72847.html)
plot(outflank1$outflank$results$FST,outflank1$outflank$results$FSTNoCorr)
abline(0,1)
text(outflank1$outflank$results$FST>0.15, outflank1$outflank$results$FSTNoCorr>0.15, labels=outflank1$outflank$results$LocusName)
# should be correlated without outliers 
# looks like loci for which FSTNoCorr is above 0.15 are outliers and should be removed before analyses
# are re-ran (7 loci)

plot(outflank1$outflank$results$He,outflank1$outflank$results$FSTNoCorr)
#should be relatively evenly distributed

hist(outflank1$outflank$results$FSTNoCorr,breaks=seq(0,0.3, by=0.001))

write_csv(outflank1$outflank$results, "outlier_analyses_outflank1.csv")

# found 3 outliers, with this approach. Not going to use it further

#### find common outliers across the top 1%, pcadapt and bayescan ####

# pcadapt and fst 1% per regions 
outliers_seafo_nam99 <- as.data.frame(outliers_seafo_nam_fst99)
colnames(outliers_seafo_nam99)[1] <- "POS"

outliers_seafo_sa99 <- as.data.frame(outliers_seafo_sa_fst99)
colnames(outliers_seafo_sa99)[1] <- "POS"

outliers_nam_sa99 <- as.data.frame(outliers_nam_sa_fst99)
colnames(outliers_nam_sa99)[1] <- "POS"

outliers_pcadaptK7 <- as.data.frame(outliers$`data_all2$loc.names`)
colnames(outliers_pcadaptK7)[1] <- "POS"

common1 <- intersect(outliers_seafo_nam99,outliers_pcadaptK7) #23
common2 <- intersect(outliers_seafo_sa99,outliers_pcadaptK7) #23
common3 <- intersect(outliers_nam_sa99,outliers_pcadaptK7) #36

  common_all <- intersect(common2,common3) #19
  
# pcadapt and bayescan
balancing1 <- as.data.frame(balancing[,1])
colnames(balancing1)[1] <- "POS"

neutral1 <- as.data.frame(neutral[,1])
colnames(neutral1)[1] <- "POS"

common4 <- intersect(balancing1,outliers_pcadaptK7) # 205
common5 <- intersect(neutral1,outliers_pcadaptK7) # 271

write_csv(common4, "outliers_pcadapt_bayescan.csv")

# fst 1% and bayescan
common6 <- intersect(balancing1,outliers_seafo_nam99) # 235
common7 <- intersect(balancing1,outliers_seafo_sa99) # 159
common8 <- intersect(balancing1,outliers_nam_sa99) # 136
common9 <- intersect(common6,common7)#8
common10 <- intersect(common6,common8)#11
common11 <- intersect(common7,common8)#16

# fst 1% and pcadapt and bayescan
common12 <- intersect(common4,common6) #9
common13 <- intersect(common4,common7) #4
common14 <- intersect(common4,common8) #5

#### Recalculate measures of differentiation focusing only on the 19 common outlier SNPs ####
# only for fst and pcadapt due to the high number of false positives of bayescan

outliers_pcBS <- gl.keep.loc(all_gl3b, loc.list = c("JAXDDX010000007_1_9220081","JAXDDX010000011_1_4529079","JAXDDX010000021_1_3352009","JAXDDX010000044_1_1819047","JAXDDX010000131_1_543394","JAXDDX010000002_1_1172427","JAXDDX010000002_1_14125432","JAXDDX010000003_1_16804485","JAXDDX010000007_1_9220081","JAXDDX010000009_1_805903","JAXDDX010000009_1_6832665","JAXDDX010000012_1_1590355","JAXDDX010000023_1_1956043","JAXDDX010000035_1_1945829","JAXDDX010000045_1_1702739","JAXDDX010000047_1_2461939","JAXDDX010000095_1_1589864","JAXDDX010000139_1_916301","JAXDDX010000176_1_327736"))
outliers_pcBS_gi <- gl2gi(outliers_pcBS)

data_pca_sc5 <- scaleGen(outliers_pcBS_gi, NA.method="mean",scale=F)
data_pca5 <- dudi.pca(data_pca_sc5, scale=F, nf= 10, scannf = F)

fviz_pca_ind(data_pca5,label="none",habillage = outliers_pcBS_gi$pop,
             mean.point=F,pointsize=4,palette = col_location,
             addEllipses = FALSE, ellipse.level=0.5,ellipse.alpha = 0, title="Orange roughy") + 
  theme_classic() + scale_shape_manual(values = c(17,17,17,19,19,19,19,19,15))


outliers_fst <- stamppFst(outliers_pcBS, nboots = 10000, percent = 95)

write.table(outliers_fst$Fsts, "outliers_fst.txt")
write.table(outliers_fst$Pvalues, "outliers_pvalues.txt")


# dapc
dapc3 <- dapc(outliers_pcBS_gi,var.contrib = TRUE, scale = FALSE, n.pca = NULL, n.da=NULL) #data is a genind object  

# n.da = nPop(dedup_gi2) - 1

scatter(dapc3, scree.da = TRUE, bg = "white", 
        posi.pca = "bottomleft", legend = TRUE, 
        cellipse = 1, cstar = 0,
        col = c(col_location),
        pch = 16,
        txt.leg = paste(levels(outliers_pcadapt3b$pop)), 
        clabel = 1, cex = 2)


# cross-validate dapc
dedup_gi3_xval <- xvalDapc(tab(outliers_pcBS_gi, NA.method = "mean"), pop(outliers_pcBS_gi))

nb_k <- find.clusters(outliers_pcadapt3b_gl)
kmeans_dapc <- dapc(outliers_pcadapt3b_gl, nb_k$grp)

scatter(kmeans_dapc, scree.da = TRUE, bg = "white", 
        posi.pca = "bottomleft", legend = TRUE, 
        cellipse = 1, cstar = 0,
        col = c(col_location),
        pch = 16,
        txt.leg = paste(levels(outliers_pcadapt3b$pop)), 
        clabel = 1, cex = 2)

compoplot(kmeans_dapc, 
          #subset = 1:65, 
          col.pal = col_location,
          border = NA,
          #col = c("blue","olivegreen", "grey"),
          #lab = "", border = NA,
          legend = F,
          show.lab = F,
          posi = "top")

# admixture
gl2geno(outliers_pcBS, outfile = "gl_outlier_pcadapt_geno", outpath = "C:/Users/u05086605/Desktop/PROJECTS/SEAFO_orange_roughy_2022", verbose = NULL)

lea_pca2 <- pca("gl_outlier_pcadapt_geno.lfmm", scale=FALSE)
tw2 <- tracy.widom(lea_pca2)
plot(tw2$percentage)

project2=NULL
project2 = snmf("gl_outlier_pcadapt_geno.geno",
               K = 1:10,
               entropy = TRUE,
               repetitions = 100,project = "new")

# check cross-entropy
plot(project2, col = "blue", pch = 19, cex = 1.2)

# plot ancestry proportions
best2 = which.min(cross.entropy(project2, K = 3))
my.colors <- c("grey","palegreen4","darkgoldenrod2")
barplot(t(Q(project2, K = 3, run = best)), 
        col = my.colors,
        xlab = "Individuals",
        ylab = "Ancestry proportions")

res = G(project2, K = 3, run = best)
