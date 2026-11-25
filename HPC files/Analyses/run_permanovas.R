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
pathmod <- adonis2(SLpath_ATCH ~ ecoregion*habitat*fire + subplot + compartment + year + season + tmean_Annual + ppt_mean_Annual, data=pathmeta_perm, by="terms")
pathmod

anova(betadisper(SLpath_ATCH, pathmeta_perm$compartment))
anova(betadisper(SLpath_ATCH, pathmeta_perm$habitat))
anova(betadisper(SLpath_ATCH, pathmeta_perm$ecoregion))
anova(betadisper(SLpath_ATCH, pathmeta_perm$fire))
anova(betadisper(SLpath_ATCH, pathmeta_perm$subplot))
anova(betadisper(SLpath_ATCH, pathmeta_perm$year))
anova(betadisper(SLpath_ATCH, pathmeta_perm$season))
anova(betadisper(SLpath_ATCH, pathmeta_perm$tmean_Annual))
anova(betadisper(SLpath_ATCH, pathmeta_perm$ppt_mean_Annual))


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

sapSLmeta <- metadata[((metadata$allsapsdiv)>0),] # 1513
sapSLmeta_filt <- sapSLmeta[(sapSLmeta$OTU_ID %in% rownames(SLsap_rclr4)),]

sapmod <- adonis2(SLsap_ATCH ~ ecoregion*habitat*fire + subplot + compartment + year + season + tmean_Annual + ppt_mean_Annual, data=sapmeta_perm, by="terms")
sapmod

anova(betadisper(SLsap_ATCH, sapmeta_perm$compartment))
anova(betadisper(SLsap_ATCH, sapmeta_perm$habitat))
anova(betadisper(SLsap_ATCH, sapmeta_perm$ecoregion))
anova(betadisper(SLsap_ATCH, sapmeta_perm$fire))
anova(betadisper(SLsap_ATCH, sapmeta_perm$subplot))
anova(betadisper(SLsap_ATCH, sapmeta_perm$year))
anova(betadisper(SLsap_ATCH, sapmeta_perm$season))
anova(betadisper(SLsap_ATCH, sapmeta_perm$tmean_Annual))
anova(betadisper(SLsap_ATCH, sapmeta_perm$ppt_mean_Annual))

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

ectomod <- adonis2(SLecto_ATCH ~ ecoregion*habitat*fire + subplot + compartment + year + season + tmean_Annual + ppt_mean_Annual, data=ectometa_perm, by="terms")
ectomod

anova(betadisper(SLecto_ATCH, ectometa_perm$compartment))
anova(betadisper(SLecto_ATCH, ectometa_perm$habitat))
anova(betadisper(SLecto_ATCH, ectometa_perm$ecoregion))
anova(betadisper(SLecto_ATCH, ectometa_perm$fire))
anova(betadisper(SLecto_ATCH, ectometa_perm$subplot))
anova(betadisper(SLecto_ATCH, ectometa_perm$year))
anova(betadisper(SLecto_ATCH, ectometa_perm$season))
anova(betadisper(SLecto_ATCH, ectometa_perm$tmean_Annual))
anova(betadisper(SLecto_ATCH, ectometa_perm$ppt_mean_Annual))


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

pyromod <- adonis2(SLpyro_ATCH ~ ecoregion*habitat*fire + subplot + compartment + year + season + tmean_Annual + ppt_mean_Annual, data=pyrometa_perm, by="terms")
pyromod

anova(betadisper(SLpyro_ATCH, pyrometa_perm$compartment))
anova(betadisper(SLpyro_ATCH, pyrometa_perm$habitat))
anova(betadisper(SLpyro_ATCH, pyrometa_perm$ecoregion))
anova(betadisper(SLpyro_ATCH, pyrometa_perm$fire))
anova(betadisper(SLpyro_ATCH, pyrometa_perm$subplot))
anova(betadisper(SLpyro_ATCH, pyrometa_perm$year))
anova(betadisper(SLpyro_ATCH, pyrometa_perm$season))
anova(betadisper(SLpyro_ATCH, pyrometa_perm$tmean_Annual))
anova(betadisper(SLpyro_ATCH, pyrometa_perm$ppt_mean_Annual))
