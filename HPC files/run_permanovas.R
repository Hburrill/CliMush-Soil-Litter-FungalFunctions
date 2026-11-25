# Packages ----
library(vegan)
library(ggplot2)
library(dplyr)
library(base)
set.seed(123)

# import data ----
# this has functional group data (richness, diversity, RA)
metadata1 <- read.csv("CliMush-SoilLitter-metaWfunctions_13Aug25.csv", header=T, sep=",")
metadata1 <- metadata1[,-1]

metadata2 <- read.csv("final_data.csv", header=T, sep=",")
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


# import SL_FT
SL_FT <- read.csv("CliMush_SL_otu_functions1.csv", head=T, sep=",")
SL_FT <- SL_FT[,-1]


# pathogens ----
SLpath_filt <- SL_FT[which(SL_FT$ppath2=="1"),]
SLpath_asv <- SLpath_filt[,-c(1509:1548)]
rownames(SLpath_asv) <- SLpath_asv$OTU_ID
SLpath_asv$OTU_ID <- NULL
SLpath_asvT <- t(SLpath_asv)
SLpath_asvT <- SLpath_asvT[rownames(SLpath_asvT) != "illumina_litter_2023.10_D16_BG_09", ]
tot <- rowSums(SLpath_asvT)
SLpath_asvT <- SLpath_asvT[tot > 0, ]

SLpath_rclr <- transform(x=SLpath_asvT, "rclr")
SLpath_rclr2 <- mutate_all(SLpath_rclr, function(x) as.numeric(as.character(x)))
SLpath_rclr3 <- SLpath_rclr2[rowSums(SLpath_rclr2[, -1])>0, ]
SLpath_rclr4 <- SLpath_rclr3[-c(1422:1423),-1]
SLpath_ATCH <- vegdist(SLpath_rclr4, method="robust.aitchison")
#writeLines(rownames(SLpath_rclr4), "SLpath_rclr4_names2.txt")

pathnames <- read.table("SLpath_rclr4_names2.txt", header=F)
pathmeta_perm <- metadata[(metadata$OTU_ID %in% pathnames$V1),]

# check
# metadata$OTU_ID[!metadata$OTU_ID %in% pathnames$V1]
# pathnames$V1[!pathnames$V1 %in% metadata$OTU_ID]

#attr(SLpath_ATCH, "Labels") <- rownames(SLpath_rclr4)
#saveRDS(SLpath_ATCH, file = "SLpath_ATCH2.rds")

#write.csv(as.data.frame(SLpath_ATCH), "SLpath_ATCH.csv")

#str(SLpathmeta_filt)
#str(SLpath_ATCH)
#pathmod <- adonis2(SLpath_ATCH ~ ecoregion*habitat*fire + subplot + compartment + year + season + tmean_Annual + ppt_mean_Annual, data=pathmeta_perm, by="terms")
#pathmod

path_comp <- betadisper(SLpath_ATCH, pathmeta_perm$compartment)
anova(path_comp)
path_comp_tab <- data.frame(
  Compartment = names(path_comp$group.distances),
  distance_to_centroid = path_comp$group.distances
)
path_comp_tab
boxplot(path_comp, xlab="Compartment", ylab="Distance to centroid")
TukeyHSD(path_comp)

path_hab <- betadisper(SLpath_ATCH, pathmeta_perm$habitat)
anova(path_hab)
path_hab_tab <- data.frame(
  Compartment = names(path_hab$group.distances),
  distance_to_centroid = path_hab$group.distances
)
path_hab_tab
boxplot(path_hab, xlab="Habitat", ylab="Distance to centroid")
TukeyHSD(path_hab)

path_eco <- betadisper(SLpath_ATCH, pathmeta_perm$ecoregion)
anova(path_eco)
path_eco_tab <- data.frame(
  Compartment = names(path_eco$group.distances),
  distance_to_centroid = path_eco$group.distances
)
path_eco_tab
boxplot(path_eco, xlab="Ecoregion", ylab="Distance to centroid")
TukeyHSD(path_eco)

path_fire <- betadisper(SLpath_ATCH, pathmeta_perm$fire)
anova(path_fire)
path_fire_tab <- data.frame(
  Compartment = names(path_fire$group.distances),
  distance_to_centroid = path_fire$group.distances
)
path_fire_tab
boxplot(path_fire, xlab="Fire", ylab="Distance to centroid")
TukeyHSD(path_fire)

path_sub <- betadisper(SLpath_ATCH, pathmeta_perm$subplot)
anova(path_sub)
path_sub_tab <- data.frame(
  Compartment = names(path_sub$group.distances),
  distance_to_centroid = path_sub$group.distances
)
path_sub_tab
boxplot(path_sub, xlab="Subplot", ylab="Distance to centroid")
TukeyHSD(path_sub)

path_year <- betadisper(SLpath_ATCH, pathmeta_perm$year)
anova(path_year)
path_year_tab <- data.frame(
  Compartment = names(path_year$group.distances),
  distance_to_centroid = path_year$group.distances
)
path_year_tab
boxplot(path_year, xlab="Year", ylab="Distance to centroid")
TukeyHSD(path_year)

path_season <- betadisper(SLpath_ATCH, pathmeta_perm$season)
anova(path_season)
path_season_tab <- data.frame(
  Compartment = names(path_season$group.distances),
  distance_to_centroid = path_season$group.distances
)
path_season_tab
boxplot(path_season, xlab="Season", ylab="Distance to centroid")
TukeyHSD(path_season)

path_mat <- betadisper(SLpath_ATCH, pathmeta_perm$tmean_Annual)
anova(path_mat)
path_mat_tab <- data.frame(
  Compartment = names(path_mat$group.distances),
  distance_to_centroid = path_mat$group.distances
)
path_mat_tab
boxplot(path_mat, xlab="mat", ylab="Distance to centroid")
TukeyHSD(path_mat)

path_map <- betadisper(SLpath_ATCH, pathmeta_perm$ppt_mean_Annual)
anova(path_map)
path_map_tab <- data.frame(
  Compartment = names(path_map$group.distances),
  distance_to_centroid = path_map$group.distances
)
path_map_tab
boxplot(path_map, xlab="map", ylab="Distance to centroid")
TukeyHSD(path_map)


# saprotrophs ----
SL_FT$allsaps <- ifelse(SL_FT$primary_lifestyle=="dung_saprotroph" | SL_FT$primary_lifestyle=="litter_saprotroph" | SL_FT$primary_lifestyle=="soil_saprotroph" | SL_FT$primary_lifestyle=="wood_saprotroph", "1", "0")

SLsap_filt <- SL_FT[which(SL_FT$allsaps=="1"),]
SLsap_asv <- SLsap_filt[,-c(1509:1545)]
SLsap_asvT <- t(SLsap_asv)
colnames(SLsap_asvT) <- SLsap_asvT[1,]
SLsap_asvT <- SLsap_asvT[-1,]
SLsap_asvT <- SLsap_asvT[rownames(SLsap_asvT) != "illumina_litter_2023.10_D16_BG_09", ]
tot <- rowSums(SLsap_asvT)
SLsap_asvT <- SLsap_asvT[tot > 0, ]
SLsap_rclr <- transform(x=SLsap_asvT, "rclr")
SLsap_rclr <- SLsap_rclr[,-1]
SLsap_rclr2 <- mutate_all(SLsap_rclr, function(x) as.numeric(as.character(x)))
SLsap_rclr3 <- SLsap_rclr2[rowSums(SLsap_rclr2[, -1])>0, ]
SLsap_rclr4 <- SLsap_rclr3[-c(1422:1423),-1]
SLsap_ATCH <- vegdist(SLsap_rclr4, method="robust.aitchison")
#writeLines(rownames(SLsap_rclr4), "SLsap_rclr4_names2.txt")

sapnames <- read.table("SLsap_rclr4_names2.txt", header=F)
sapmeta_perm <- metadata[(metadata$OTU_ID %in% sapnames$V1),]

# check
# metadata$OTU_ID[!metadata$OTU_ID %in% sapnames$V1]
# sapnames$V1[!sapnames$V1 %in% metadata$OTU_ID]
#attr(SLsap_ATCH, "Labels") <- rownames(SLsap_rclr4)
#saveRDS(SLsap_ATCH, file = "SLsap_ATCH2.rds")

#sapSLmeta <- metadata[((metadata$allsapsdiv)>0),] # 1513
#sapSLmeta_filt <- sapSLmeta[(sapSLmeta$OTU_ID %in% rownames(SLsap_rclr4)),]

#sapmod <- adonis2(SLsap_ATCH ~ ecoregion*habitat*fire + subplot + compartment + year + season + tmean_Annual + ppt_mean_Annual, data=sapmeta_perm, by="terms")
#sapmod

sap_comp <- betadisper(betadisper(SLsap_ATCH, sapmeta_perm$compartment))
anova(sap_comp)
sap_comp_tab <- data.frame(
  Compartment = names(sap_comp$group.distances),
  distance_to_centroid = sap_comp$group.distances
)
boxplot(sap_comp, xlab="Compartment", ylab="Distance to centroid")
TukeyHSD(sap_comp)

sap_hab <- betadisper(SLsap_ATCH, sapmeta_perm$habitat)
anova(sap_hab)
sap_hab_tab <- data.frame(
  Compartment = names(sap_hab$group.distances),
  distance_to_centroid = sap_hab$group.distances
)
sap_hab_tab
boxplot(sap_hab, xlab="Habitat", ylab="Distance to centroid")
TukeyHSD(sap_hab)

sap_eco <- betadisper(SLsap_ATCH, sapmeta_perm$ecoregion)
anova(sap_eco)
sap_eco_tab <- data.frame(
  Compartment = names(sap_eco$group.distances),
  distance_to_centroid = sap_eco$group.distances
)
sap_eco_tab
boxplot(sap_eco, xlab="Ecoregion", ylab="Distance to centroid")
TukeyHSD(sap_eco)

sap_fire <- betadisper(SLsap_ATCH, sapmeta_perm$fire)
anova(sap_fire)
sap_fire_tab <- data.frame(
  Compartment = names(sap_fire$group.distances),
  distance_to_centroid = sap_fire$group.distances
)
sap_fire_tab
boxplot(sap_fire, xlab="Fire", ylab="Distance to centroid")
TukeyHSD(sap_fire)

sap_sub <- betadisper(SLsap_ATCH, sapmeta_perm$subplot)
anova(sap_sub)
sap_sub_tab <- data.frame(
  Compartment = names(sap_sub$group.distances),
  distance_to_centroid = sap_sub$group.distances
)
sap_sub_tab
boxplot(sap_sub, xlab="Subplot", ylab="Distance to centroid")
TukeyHSD(sap_sub)

sap_year <- betadisper(SLsap_ATCH, sapmeta_perm$year)
anova(sap_year)
sap_year_tab <- data.frame(
  Compartment = names(sap_year$group.distances),
  distance_to_centroid = sap_year$group.distances
)
sap_year_tab
boxplot(sap_year, xlab="Year", ylab="Distance to centroid")
TukeyHSD(sap_year)

sap_season <- betadisper(SLsap_ATCH, sapmeta_perm$season)
anova(sap_season)
sap_season_tab <- data.frame(
  Compartment = names(sap_season$group.distances),
  distance_to_centroid = sap_season$group.distances
)
sap_season_tab
boxplot(sap_season, xlab="Season", ylab="Distance to centroid")
TukeyHSD(sap_season)

sap_mat <- betadisper(SLsap_ATCH, sapmeta_perm$tmean_Annual)
anova(sap_mat)
sap_mat_tab <- data.frame(
  Compartment = names(sap_mat$group.distances),
  distance_to_centroid = sap_mat$group.distances
)
sap_mat_tab
boxplot(sap_mat, xlab="mat", ylab="Distance to centroid")
TukeyHSD(sap_mat)

sap_map <- betadisper(SLsap_ATCH, sapmeta_perm$ppt_mean_Annual)
anova(sap_map)
sap_map_tab <- data.frame(
  Compartment = names(sap_map$group.distances),
  distance_to_centroid = sap_map$group.distances
)
sap_map_tab
boxplot(sap_map, xlab="map", ylab="Distance to centroid")
TukeyHSD(sap_map)

# ectomycorrhizae ----
SLecto_filt <- SL_FT[which(SL_FT$ecto=="1"),]
SLecto_asv <- SLecto_filt[,-c(1509:1545)]
SLecto_asvT <- t(SLecto_asv)
colnames(SLecto_asvT) <- SLecto_asvT[1,]
SLecto_asvT <- SLecto_asvT[rownames(SLecto_asvT) != "illumina_litter_2023.10_D16_BG_09" & rownames(SLecto_asvT) != "illumina_soil_2022.10_D05_UC_02", ]
tot <- rowSums(SLecto_asvT)
SLecto_asvT <- SLecto_asvT[tot > 0, ]
SLecto_rclr <- transform(x=SLecto_asvT, "rclr")
SLecto_rclr2 <- SLecto_rclr[,-1]
SLecto_rclr3 <- mutate_all(SLecto_rclr2, function(x) as.numeric(as.character(x)))
SLecto_rclr4 <- SLecto_rclr2[rowSums(SLecto_rclr3)>0, ]
#SLecto_rclr5 <- mutate_all(SLecto_rclr4, function(x) as.numeric(as.character(x)))
SLecto_ATCH <- vegdist(SLecto_rclr4, method="robust.aitchison")
#writeLines(rownames(SLecto_rclr4), "SLecto_rclr4_names2.txt")

ectonames <- read.table("SLecto_rclr4_names2.txt", header=F)
ectometa_perm <- metadata[(metadata$OTU_ID %in% ectonames$V1),]

# check
#metadata$OTU_ID[!metadata$OTU_ID %in% ectonames$V1]
#ectonames$V1[!ectonames$V1 %in% metadata$OTU_ID]
#attr(SLecto_ATCH, "Labels") <- rownames(SLecto_rclr4)
#saveRDS(SLecto_ATCH, file = "SLecto_ATCH2.rds")

ectoSLmeta <- metadata[((metadata$ectorich)>0),] # 1513
ectoSLmeta_filt <- ectoSLmeta[(ectoSLmeta$OTU_ID %in% rownames(SLecto_asvT)),]

#ectomod <- adonis2(SLecto_ATCH ~ ecoregion*habitat*fire + subplot + compartment + year + season + tmean_Annual + ppt_mean_Annual, data=ectometa_perm, by="terms")
#ectomod

ecto_comp <- betadisper(SLecto_ATCH, ectometa_perm$compartment)
anova(ecto_comp)
ecto_comp_tab <- data.frame(
  Compartment = names(ecto_comp$group.distances),
  distance_to_centroid = ecto_comp$group.distances
)
boxplot(ecto_comp, xlab="Compartment", ylab="Distance to centroid")
TukeyHSD(ecto_comp)

ecto_hab <- betadisper(SLecto_ATCH, ectometa_perm$habitat)
anova(ecto_hab)
ecto_hab_tab <- data.frame(
  Compartment = names(ecto_hab$group.distances),
  distance_to_centroid = ecto_hab$group.distances
)
ecto_hab_tab
boxplot(ecto_hab, xlab="Habitat", ylab="Distance to centroid")
TukeyHSD(ecto_hab)

ecto_eco <- betadisper(SLecto_ATCH, ectometa_perm$ecoregion)
anova(ecto_eco)
ecto_eco_tab <- data.frame(
  Compartment = names(ecto_eco$group.distances),
  distance_to_centroid = ecto_eco$group.distances
)
ecto_eco_tab
boxplot(ecto_eco, xlab="Ecoregion", ylab="Distance to centroid")
TukeyHSD(ecto_eco)

ecto_fire <- betadisper(SLecto_ATCH, ectometa_perm$fire)
anova(ecto_fire)
ecto_fire_tab <- data.frame(
  Compartment = names(ecto_fire$group.distances),
  distance_to_centroid = ecto_fire$group.distances
)
ecto_fire_tab
boxplot(ecto_fire, xlab="Fire", ylab="Distance to centroid")
TukeyHSD(ecto_fire)

ecto_sub <- betadisper(SLecto_ATCH, ectometa_perm$subplot)
anova(ecto_sub)
ecto_sub_tab <- data.frame(
  Compartment = names(ecto_sub$group.distances),
  distance_to_centroid = ecto_sub$group.distances
)
ecto_sub_tab
boxplot(ecto_sub, xlab="Subplot", ylab="Distance to centroid")
TukeyHSD(ecto_sub)

ecto_year <- betadisper(SLecto_ATCH, ectometa_perm$year)
anova(ecto_year)
ecto_year_tab <- data.frame(
  Compartment = names(ecto_year$group.distances),
  distance_to_centroid = ecto_year$group.distances
)
ecto_year_tab
boxplot(ecto_year, xlab="Year", ylab="Distance to centroid")
TukeyHSD(ecto_year)

ecto_season <- betadisper(SLecto_ATCH, ectometa_perm$season)
anova(ecto_season)
ecto_season_tab <- data.frame(
  Compartment = names(ecto_season$group.distances),
  distance_to_centroid = ecto_season$group.distances
)
ecto_season_tab
boxplot(ecto_season, xlab="Season", ylab="Distance to centroid")
TukeyHSD(ecto_season)

ecto_mat <- betadisper(SLecto_ATCH, ectometa_perm$tmean_Annual)
anova(ecto_mat)
ecto_mat_tab <- data.frame(
  Compartment = names(ecto_mat$group.distances),
  distance_to_centroid = ecto_mat$group.distances
)
ecto_mat_tab
boxplot(ecto_mat, xlab="mat", ylab="Distance to centroid")
TukeyHSD(ecto_mat)

ecto_map <- betadisper(SLecto_ATCH, ectometa_perm$ppt_mean_Annual)
anova(ecto_map)
ecto_map_tab <- data.frame(
  Compartment = names(ecto_map$group.distances),
  distance_to_centroid = ecto_map$group.distances
)
ecto_map_tab
boxplot(ecto_map, xlab="map", ylab="Distance to centroid")
TukeyHSD(ecto_map)

# pyrofiles ----
SLpyro_filt <- SL_FT[which(SL_FT$pyro=="1"),]
SLpyro_asv <- SLpyro_filt[,-c(1509:1545)]
rownames(SLpyro_asv) <- SLpyro_asv$OTU_ID
SLpyro_asv$OTU_ID <- NULL
SLpyro_asvT <- t(SLpyro_asv)
SLpyro_asvT <- SLpyro_asvT[rownames(SLpyro_asvT) != "illumina_litter_2023.10_D16_BG_09", ]
tot <- rowSums(SLpyro_asvT)
SLpyro_asvT <- SLpyro_asvT[tot > 0, ]
SLpyro_rclr <- transform(x=SLpyro_asvT, "rclr")
SLpyro_rclr2 <- SLpyro_rclr[,-1]
SLpyro_rclr3 <- mutate_all(SLpyro_rclr2, function(x) as.numeric(as.character(x)))
SLpyro_rclr4 <- SLpyro_rclr2[rowSums(SLpyro_rclr3)>0, ]
SLpyro_ATCH <- vegdist(SLpyro_rclr2[,-1], method="robust.aitchison")
#writeLines(rownames(SLpyro_rclr4), "SLpyro_rclr4_names2.txt")

pyronames <- read.table("SLpyro_rclr4_names2.txt", header=F)
pyrometa_perm <- metadata[(metadata$OTU_ID %in% pyronames$V1),]

# check
#metadata$OTU_ID[!metadata$OTU_ID %in% pyronames$V1]
#pyronames$V1[!pyronames$V1 %in% metadata$OTU_ID]
#attr(SLpyro_ATCH, "Labels") <- rownames(SLpyro_rclr4)
#saveRDS(SLpyro_ATCH, file = "SLpyro_ATCH2.rds")

SLpyrometa <- metadata[((metadata$pyrorich)>0),]
SLpyrometa_filt <- SLpyrometa[(SLpyrometa$OTU_ID %in% rownames(SLpyro_asvT)),]

raremax <- min(rowSums(SLpyro_asvT))
pyro_r <- rrarefy(SLpyro_asvT, raremax)

#pyromod <- adonis2(SLpyro_ATCH ~ ecoregion*habitat*fire + subplot + compartment + year + season + tmean_Annual + ppt_mean_Annual, data=pyrometa_perm, by="terms")
#pyromod

pyro_comp <- betadisper(SLpyro_ATCH, pyrometa_perm$compartment)
anova(pyro_comp)
pyro_comp_tab <- data.frame(
  Compartment = names(pyro_comp$group.distances),
  distance_to_centroid = pyro_comp$group.distances
)
boxplot(pyro_comp, xlab="Compartment", ylab="Distance to centroid")
TukeyHSD(pyro_comp)

pyro_hab <- betadisper(SLpyro_ATCH, pyrometa_perm$habitat)
anova(pyro_hab)
pyro_hab_tab <- data.frame(
  Compartment = names(pyro_hab$group.distances),
  distance_to_centroid = pyro_hab$group.distances
)
pyro_hab_tab
boxplot(pyro_hab, xlab="Habitat", ylab="Distance to centroid")
TukeyHSD(pyro_hab)

pyro_eco <- betadisper(SLpyro_ATCH, pyrometa_perm$ecoregion)
anova(pyro_eco)
pyro_eco_tab <- data.frame(
  Compartment = names(pyro_eco$group.distances),
  distance_to_centroid = pyro_eco$group.distances
)
pyro_eco_tab
boxplot(pyro_eco, xlab="Ecoregion", ylab="Distance to centroid")
TukeyHSD(pyro_eco)

pyro_fire <- betadisper(SLpyro_ATCH, pyrometa_perm$fire)
anova(pyro_fire)
pyro_fire_tab <- data.frame(
  Compartment = names(pyro_fire$group.distances),
  distance_to_centroid = pyro_fire$group.distances
)
pyro_fire_tab
boxplot(pyro_fire, xlab="Fire", ylab="Distance to centroid")
TukeyHSD(pyro_fire)

pyro_sub <- betadisper(SLpyro_ATCH, pyrometa_perm$subplot)
anova(pyro_sub)
pyro_sub_tab <- data.frame(
  Compartment = names(pyro_sub$group.distances),
  distance_to_centroid = pyro_sub$group.distances
)
pyro_sub_tab
boxplot(pyro_sub, xlab="Subplot", ylab="Distance to centroid")
TukeyHSD(pyro_sub)

pyro_year <- betadisper(SLpyro_ATCH, pyrometa_perm$year)
anova(pyro_year)
pyro_year_tab <- data.frame(
  Compartment = names(pyro_year$group.distances),
  distance_to_centroid = pyro_year$group.distances
)
pyro_year_tab
boxplot(pyro_year, xlab="Year", ylab="Distance to centroid")
TukeyHSD(pyro_year)

pyro_season <- betadisper(SLpyro_ATCH, pyrometa_perm$season)
anova(pyro_season)
pyro_season_tab <- data.frame(
  Compartment = names(pyro_season$group.distances),
  distance_to_centroid = pyro_season$group.distances
)
pyro_season_tab
boxplot(pyro_season, xlab="Season", ylab="Distance to centroid")
TukeyHSD(pyro_season)

pyro_mat <- betadisper(SLpyro_ATCH, pyrometa_perm$tmean_Annual)
anova(pyro_mat)
pyro_mat_tab <- data.frame(
  Compartment = names(pyro_mat$group.distances),
  distance_to_centroid = pyro_mat$group.distances
)
pyro_mat_tab
boxplot(pyro_mat, xlab="mat", ylab="Distance to centroid")
TukeyHSD(pyro_mat)

pyro_map <- betadisper(SLpyro_ATCH, pyrometa_perm$ppt_mean_Annual)
anova(pyro_map)
pyro_map_tab <- data.frame(
  Compartment = names(pyro_map$group.distances),
  distance_to_centroid = pyro_map$group.distances
)
pyro_map_tab
boxplot(pyro_map, xlab="map", ylab="Distance to centroid")
TukeyHSD(pyro_map)