# Packages ----
library(vegan)
library(ggplot2)
library(ggpubr)
#library(microViz)
library(AICcPermanova)
library(dplyr)
library(ggordiplots)
library(base)
library(glmulti)
library(lme4)
library(Matrix)
library(CCA)
library(reshape2)
library(indicspecies)
library(cowplot)
set.seed(123)

# import data ----
# this has functional group data (richness, diversity, RA)
metadata1 <- read.csv("CliMush-SoilLitter-metaWfunctions_13Aug25.csv", header=T, sep=",")
head(metadata1)
metadata1 <- metadata1[,-1]

metadata2 <- read.csv("final_data.csv", header=T, sep=",")
head(metadata2)
# change habitat to C, G, O
# change fire to B, UB
metadata2 <- metadata2 %>%
  mutate(habitat=recode(habitat,
                       "Conifer"="C",
                       "Oak"="O",
                       "Grassland"="G")) %>%
  mutate(fire=recode(fire,
                       "Burned"="B",
                       "Unburned"="UB"))

# add OTU_ID to metadata
metadata <- metadata2 %>%
  left_join(
    metadata1 %>% dplyr::select(domain, habitat, season, fire, year, subplot, compartment, OTU_ID),
    by = c("domain", "habitat", "season", "fire", "year", "compartment", "subplot"), relationship =
      "many-to-many"
  )

metadata$ecoregion <- as.factor(metadata$ecoregion)
metadata$habitat <- as.factor(metadata$habitat)
metadata$subplot <- as.factor(metadata$subplot)
metadata$compartment <- as.factor(metadata$compartment)
metadata$fire <- as.factor(metadata$fire)
metadata$season <- as.factor(metadata$season)
metadata$habitat <- as.factor(metadata$habitat)

#metadata1 %>%
#  count(domain, habitat, season, fire, year, compartment) %>%
#  filter(n > 1)

#md1_collapsed <- metadata1 %>%
#  group_by(domain, habitat, season, fire, year, compartment) %>%
#  summarise(OTU_ID = paste(unique(OTU_ID), collapse = "; "), .groups = "drop")

#metadata <- metadata2 %>%
#  left_join(md1_collapsed, by = c("domain", "habitat", "season", "fire", "year", "compartment"))


# import SL_FT
#SL_FT <- read.csv("CliMush_SL_otu_functions1.csv", head=T, sep=",")
#head(SL_FT)
#SL_FT <- SL_FT[,-1]
## subset for tests
#SL_FT <- SL_FT[-c(20000:821300),]
#write.csv(SL_FT, "SL_FT_subset.csv")
#SL_FT <- read.csv("SL_FT_subset.csv", head=T, sep=",")
#head(SL_FT)
#SL_FT <- SL_FT[,-1]

#SL_FT2 <- mutate_all(SL_FT, function(x) as.numeric(as.character(x)))
#SL_FT3 <- SL_FT2[,-c(1,1509:1543)]
#otu.table.filter <- SL_FT3[,colSums(SL_FT3)>10]

# pathogens ----
pathrds <- readRDS("SLpath_ATCH2.rds")
pathnames <- read.table("SLpath_rclr4_names2.txt", header=F)
# permanova
#SL_FT3_t <- t(SL_FT3)
pathmeta_perm <- metadata[(metadata$OTU_ID %in% pathnames$V1),]
head(pathmeta_perm)

# check
metadata$OTU_ID[!metadata$OTU_ID %in% pathnames$V1]
pathnames$V1[!pathnames$V1 %in% metadata$OTU_ID]

#SLpathmeta_sub <- SLpathmeta[common_samples, , drop = FALSE]

#adonis2(pathrds ~ ecoregion + habitat + fire + subplot + compartment + year + season + tmean_Annual + ppt_mean_Annual, data=pathmeta_perm, by="terms") #  in order to include elevation, pH, CN, and moisture, need to remove all samples where NA, then re-run vegdist on that subset

# from https://chrischizinski.github.io/rstats/vegan-ggplot2/
path.site <- read.csv("pathsite_scores2.csv", header=T)
head(path.site)

SLpathmeta <- metadata[((metadata$pprich1)>0),]
SLpathmeta_filt <- SLpathmeta[(SLpathmeta$OTU_ID %in% path.site$X),]
head(SLpathmeta_filt)
path.site$OTU_ID <- path.site$X

SLpathmeta_mds <- left_join(SLpathmeta_filt, path.site, by="OTU_ID")

# NMDS plots
path1 <- ggplot(data=SLpathmeta_mds,aes(x=NMDS2,y=NMDS3)) +  geom_point(data=SLpathmeta_mds,aes(x=NMDS2,y=NMDS3,colour=habitat),size=1.5,alpha=0.5) + 
  scale_colour_manual(values=c("C" = "darkgreen", "G" = "orange", "O"="brown")) +
  coord_equal() +
  theme_bw() + theme(legend.position = "none")

path2 <- ggplot(data=SLpathmeta_mds, aes(x=NMDS2,y=NMDS3)) + geom_point(data=SLpathmeta_mds,aes(x=NMDS2,y=NMDS3,colour=compartment),size=1.5,alpha=0.5) + 
  scale_colour_manual(values=c("soil" = "darkmagenta", "litter" = "darkorange")) +
  theme_bw() +
  labs(title = "") + theme(legend.position = "none")
#+ ylim(-0.3,0.4) + xlim(2.75,4)

path3 <- ggplot(data=SLpathmeta_mds, aes(x=NMDS2,y=NMDS3)) + geom_point(data=SLpathmeta_mds,aes(x=NMDS2,y=NMDS3, colour=fire),size=1.5,alpha=0.5) + 
  scale_colour_manual(values=c("B" = "darkred", "UB" = "orange")) +
  theme_bw() +
  labs(title = "Pathogen NMDS by fire")
#+  ylim(-0.3,0.4) + xlim(2.75,4)

path4 <- ggplot(data=SLpathmeta_mds, aes(x=NMDS2,y=NMDS3)) + geom_point(data=SLpathmeta_mds,aes(x=NMDS2,y=NMDS3, colour=ecoregion),size=1.5,alpha=0.5) +
  coord_equal() +
  theme_bw() + theme(legend.position = "none")
     #  +  ylim(-0.3,0.4) + xlim(2.75,4)

pathNMDS <- ggarrange(path1, path2, path3, path4, 
                     labels = c("A", "B", "C", "E"),
                     ncol = 2, nrow = 2) + bgcolor("white")
ggsave("pathNMDS23.png", path="figures/", width=8, height =6)

# saprotrophs ----
sap.site <- read.csv("sapsite_scores.csv", header=T)
head(sap.site)

SLsapmeta <- metadata[((metadata$pprich1)>0),]
SLsapmeta_filt <- SLsapmeta[(SLsapmeta$OTU_ID %in% sap.site$X),]
sap.site$OTU_ID <- sap.site$X

SLsapmeta_mds <- left_join(SLsapmeta_filt, sap.site, by="OTU_ID")

# NMDS plots
sap1 <- ggplot(data=SLsapmeta_mds,aes(x=NMDS1,y=NMDS2)) + 
  geom_point(data=SLsapmeta_mds,aes(x=NMDS1,y=NMDS2,colour=habitat),size=1.5,alpha=0.5) + 
  scale_colour_manual(values=c("C" = "darkgreen", "G" = "orange", "O"="brown")) +
  coord_equal() +
  theme_bw() + theme(legend.position = "none") +
  labs(title = "                 Habitat") 

sap2 <- ggplot(data=SLsapmeta_mds, aes(x=NMDS1,y=NMDS2)) +
  geom_point(data=SLsapmeta_mds,aes(x=NMDS1,y=NMDS2,colour=compartment),size=1.5,alpha=0.5) + 
  scale_colour_manual(values=c("soil" = "darkmagenta", "litter" = "darkorange")) +
  coord_equal() +
  theme_bw() +
  labs(title = "Compartment") + theme(legend.position = "none")

sap3 <- ggplot(data=SLsapmeta_mds, aes(x=NMDS1,y=NMDS2)) +
  geom_point(data=SLsapmeta_mds,aes(x=NMDS1,y=NMDS2,colour=fire),size=1.5,alpha=0.5) + 
  scale_colour_manual(values=c("B" = "darkred", "UB" = "orange")) +
  coord_equal() +
  theme_bw() +
  labs(title = "Saprotroph NMDS by fire")

sap4 <- ggplot(data=SLsapmeta_mds, aes(x=NMDS1,y=NMDS2)) +
  geom_point(data=SLsapmeta_mds,aes(x=NMDS1,y=NMDS2,colour=ecoregion),size=1.5,alpha=0.5) + coord_equal() +
  theme_bw() + theme(legend.position = "none") +
  labs(title = "Ecoregion") 

sapNMDS <- ggarrange(sap1, sap2, sap3, sap4, 
          labels = c("A", "B", "C", "E"),
          ncol = 2, nrow = 2) + bgcolor("white")
ggsave("sapNMDS.png", path="figures/", width=8, height
=6)

# ectomycorrhizae ----
ectords <- readRDS("SLecto_ATCH2.rds")
ectonames <- read.table("SLecto_rclr4_names2.txt", header=F)
# permanova
ectometa_perm <- metadata[(metadata$OTU_ID %in% ectonames$V1),]
head(ectometa_perm)

# check
metadata$OTU_ID[!metadata$OTU_ID %in% ectonames$V1]
ectonames$V1[!ectonames$V1 %in% metadata$OTU_ID]

#SLpathmeta_sub <- SLpathmeta[common_samples, , drop = FALSE]

#adonis2(ectords ~ ecoregion + habitat + fire + subplot + compartment + year + season + tmean_Annual + ppt_mean_Annual, data=ectometa_perm, by="terms") #  in order to include elevation, pH, CN, and moisture, need to remove all samples where NA, then re-run vegdist on that subset

ecto.site <- read.csv("ectosite_scores.csv", header=T)
head(ecto.site)

ectoSLmeta <- metadata[((metadata$ectorich)>0),] # 1513
ectoSLmeta_filt <- ectoSLmeta[(ectoSLmeta$OTU_ID %in% ectonames$V1),]
ecto.site$OTU_ID <- ecto.site$X

SLectometa_mds <- left_join(ectoSLmeta_filt, ecto.site, by="OTU_ID")

# NMDS plots
ecto1 <- ggplot(data=SLectometa_mds,aes(x=NMDS1,y=NMDS2)) + geom_point(data=SLectometa_mds,aes(x=NMDS1,y=NMDS2,colour=habitat),size=1.5,alpha=0.5) + 
  scale_colour_manual(values=c("C" = "darkgreen", "G" = "orange", "O"="brown")) +
  coord_equal() +
  theme_bw() + theme(legend.position = "none")

ecto2 <- ggplot(data=SLectometa_mds, aes(x=NMDS1,y=NMDS2)) + geom_point(data=SLectometa_mds,aes(x=NMDS1,y=NMDS2,colour=compartment),size=1.5,alpha=0.5) + 
  scale_colour_manual(values=c("soil" = "darkmagenta", "litter" = "darkorange")) +
  coord_equal() +
  theme_bw() +
  labs(title = "") + 
  ylim(-0.25,0.25) + xlim(-0.25,0.25) + theme(legend.position = "none")

ecto3 <- ggplot(data=SLectometa_mds, aes(x=NMDS1,y=NMDS2)) + geom_point(data=SLectometa_mds,aes(x=NMDS1,y=NMDS2,colour=fire),size=1.5,alpha=0.5) + 
  scale_colour_manual(values=c("B" = "darkred", "UB" = "orange")) +
  coord_equal() +
  theme_bw() +
  labs(title = "Ectomycorrhizae NMDS by fire") + 
  ylim(-0.25,0.25) + xlim(-0.25,0.25)

ecto4 <- ggplot(data=SLectometa_mds, aes(x=NMDS1,y=NMDS2)) + geom_point(data=SLectometa_mds,aes(x=NMDS1,y=NMDS2,colour=ecoregion),size=1.5,alpha=0.5) + coord_equal() +
  theme_bw() + theme(legend.position = "none")

ectoNMDS <- ggarrange(ecto1, ecto2, ecto3, ecto4, 
                     labels = c("A", "B", "C", "E"),
                     ncol = 2, nrow = 2) + bgcolor("white")
ggsave("ectoNMDS_changelim.png", path="figures/", width=8, height=6)

# pyrofiles ----
pyro.site <- read.csv("pyrosite_scores.csv", header=T)
head(pyro.site)
pyro.site$OTU_ID <- pyro.site$X
pyronames <- read.table("SLpyro_rclr4_names2.txt", header=F)

SLpyrometa <- metadata[((metadata$pyrorich)>0),]
SLpyrometa_filt <- SLpyrometa[(SLpyrometa$OTU_ID %in% pyronames$V1),]
SLpyrometa_mds <- left_join(SLpyrometa_filt, pyro.site, by="OTU_ID")

# NMDS plots
pyro1 <- ggplot(data=SLpyrometa_mds,aes(x=NMDS1,y=NMDS2)) + geom_point(data=SLpyrometa_mds,aes(x=NMDS1,y=NMDS2,colour=habitat),size=1.5,alpha=0.5) + 
  scale_colour_manual(values=c("C" = "darkgreen", "G" = "orange", "O"="brown")) +
  coord_equal() +
  theme_bw() + theme(legend.position = "none")+
  ylim(-0.005,0.005) + xlim(0,0.01)

pyro2 <- ggplot(data=SLpyrometa_mds, aes(x=NMDS1,y=NMDS2)) + geom_point(data=SLpyrometa_mds,aes(x=NMDS1,y=NMDS2,colour=compartment),size=1.5,alpha=0.5) + 
  scale_colour_manual(values=c("soil" = "darkmagenta", "litter" = "darkorange")) +
  coord_equal() +
  theme_bw() +
  labs(title = "") +
  ylim(-0.005,0.005) + xlim(0,0.01) + theme(legend.position = "none")

pyro3 <- ggplot(data=SLpyrometa_mds, aes(x=NMDS1,y=NMDS2)) + geom_point(data=SLpyrometa_mds,aes(x=NMDS1,y=NMDS2,colour=fire),size=1.5,alpha=0.5) + 
  scale_colour_manual(values=c("B" = "darkred", "UB" = "orange")) +
  coord_equal() +
  theme_bw() +
  labs(title = "Pyrofile NMDS by fire") +
  ylim(-0.005,0.005) + xlim(0,0.01)

pyro4 <- ggplot(data=SLpyrometa_mds, aes(x=NMDS1,y=NMDS2)) + geom_point(data=SLpyrometa_mds,aes(x=NMDS1,y=NMDS2,colour=ecoregion),size=1.5,alpha=0.5) + coord_equal() +
  theme_bw() + theme(legend.position = "none")+
  ylim(-0.005,0.005) + xlim(0,0.01)

pyroNMDS <- ggarrange(pyro1, pyro2, pyro3, pyro4, 
                      labels = c("A", "B", "C", "E"),
                      ncol = 2, nrow = 2) + bgcolor("white")
ggsave("pyroNMDS.png", path="figures/", width=8, height =6)


# habitat, compartment, and ecoregion fig ----
leg_fill <- get_legend(path1 + theme(legend.position = "bottom"))
leg2 <- get_legend(path2 + theme(legend.position = "bottom"))
leg_pattern <- get_legend(path4 + theme(legend.position = "bottom"))


allNMDS <- ggarrange(sap1, sap2, sap4,
                     path1, path2, path4,
                      ecto1, ecto2, ecto4,
                      pyro1, pyro2, pyro4,
                      labels = c("Saprotrophs","","", "Pathogens", "","", "Ectomycorrhizae", "","", "Pyrophiles", "",""), ncol = 3, nrow = 4, align = "hv") + bgcolor("white")
ggsave("allNMDS2.png", path="figures/", width=8, height =10)

final_plot <- ggarrange(leg_fill, leg2, leg_pattern, ncol=1, nrow=3) + bgcolor("white")

ggsave("allNMDS_leg2.png", path="figures/", width=7, height =4)
