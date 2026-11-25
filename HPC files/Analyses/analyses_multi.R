# Packages ----
library(vegan)
library(ggplot2)
#library(microViz)
library(dplyr)
set.seed(123)

# import SL_FT
SL_FT <- read.csv("CliMush_SL_otu_functions1.csv", head=T, sep=",")
#head(SL_FT)
SL_FT <- SL_FT[,-1]

# pathogens ----
SLpath_filt <- SL_FT[which(SL_FT$ppath2=="1"),]
SLpath_asv <- SLpath_filt[,-c(1509:1548)]
rownames(SLpath_asv) <- SLpath_asv$OTU_ID
SLpath_asv$OTU_ID <- NULL
SLpath_asvT <- t(SLpath_asv)
tot <- rowSums(SLpath_asvT)
SLpath_asvT <- SLpath_asvT[tot > 0, ]
#SLpathmeta <- merged.meta2[((merged.meta2$pprich1)>0),]
#SLpathmeta_filt <- SLpathmeta[(SLpathmeta$OTU_ID %in% rownames(SLpath_asvT)),]

SLpath_rclr <- transform(x=SLpath_asvT, "rclr")
SLpath_rclr2 <- mutate_all(SLpath_rclr, function(x) as.numeric(as.character(x)))
SLpath_rclr3 <- SLpath_rclr2[rowSums(SLpath_rclr2[, -1])>0, ]
SLpath_rclr4 <- SLpath_rclr3[-c(1422:1423),-1]
SLpath_ATCH <- vegdist(SLpath_rclr4, method="robust.aitchison")

#pathmod <- adonis2(SLpath_ATCH ~ subplot.x + compartment.x + season.x:year.x + eco.region + fire.x*habitat.x, by="terms", data=SLpathmeta_filt)
#pathmod

#pathmod2 <- adonis2(SLpath_ATCH ~ compartment.x + iBut_MMT + Dew_mean + precip_7day + elevation + PercMoist + CN + pH, by="terms", data=SLpathmeta_filt)
#pathmod2

raremax <- min(rowSums(SLpath_asvT))
path_r <- rrarefy(SLpath_asvT, raremax)
path_nmds <- metaMDS(SLpath_rclr4, distance="bray")
stressplot(path_nmds)
plot(path_nmds)

# from https://chrischizinski.github.io/rstats/vegan-ggplot2/
path.site <- as.data.frame(scores(path_nmds)$sites)
path.site$site <- rownames(path.site)
head(path.site)
write.csv(path.site, "pathsite_scores.csv")

path.species <- as.data.frame(scores(path_nmds)$species)
path.species$OTU_ID <- rownames(path.species)
head(path.species)
write.csv(path.species, "pathOTU_scores.csv")

# saprotrophs ----
SL_FT$allsaps <- ifelse(SL_FT$primary_lifestyle=="dung_saprotroph" | SL_FT$primary_lifestyle=="litter_saprotroph" | SL_FT$primary_lifestyle=="soil_saprotroph" | SL_FT$primary_lifestyle=="wood_saprotroph" | SL_FT$Secondary_lifestyle=="dung_saprotroph" | SL_FT$Secondary_lifestyle=="litter_saprotroph" | SL_FT$Secondary_lifestyle=="soil_saprotroph" | SL_FT$Secondary_lifestyle=="wood_saprotroph", "1", "0")

SLsap_filt <- SL_FT[which(SL_FT$allsaps=="1"),]
SLsap_asv <- SLsap_filt[,-c(1509:1543)]
SLsap_asvT <- t(SLsap_asv)
colnames(SLsap_asvT) <- SLsap_asvT[1,]
SLsap_asvT <- SLsap_asvT[-1,]
tot <- rowSums(SLsap_asvT)
SLsap_asvT <- SLsap_asvT[tot > 0, ]
SLsap_rclr <- transform(x=SLsap_asvT, "rclr")
SLsap_rclr2 <- SLsap_rclr[,-1]
SLsap_rclr3 <- mutate_all(SLsap_rclr2, function(x) as.numeric(as.character(x)))
SLsap_rclr4 <- SLsap_rclr3[rowSums(SLsap_rclr3>0), ]
SLsap_ATCH <- vegdist(SLsap_rclr3, method="robust.aitchison")

#sapSLmeta <- merged.meta2[((merged.meta2$saprich)>0),] # 1513
#sapSLmeta_filt <- sapSLmeta[(sapSLmeta$OTU_ID %in% rownames(SLsap_asvT)),]

#sapmod <- adonis2(SLsap_ATCH ~ subplot.x + compartment.x + season.x:year.x + eco.region + fire.x*habitat.x, by="terms", data=sapSLmeta_filt)
#sapmod

#sapmod2 <- adonis2(SLsap_ATCH ~ compartment.x + iBut_MMT + Dew_mean + precip_7day + elevation + PercMoist + CN + pH, by="terms", data=sapSLmeta_filt)
#sapmod2

raremax <- min(rowSums(SLsap_asvT))
sap_r <- rrarefy(SLsap_asvT, raremax)
sap_nmds <- metaMDS(SLsap_rclr3, distance="bray")
sap_nmds
stressplot(sap_nmds)
plot(sap_nmds)

sap.site <- as.data.frame(scores(sap_nmds)$sites)
sap.site$site <- rownames(sap.site)
head(sap.site)
write.csv(sap.site, "sapsite_scores.csv")

sap.species <- as.data.frame(scores(sap_nmds)$species)
sap.species$OTU_ID <- rownames(sap.species)
head(sap.species)
write.csv(sap.species, "sapOTU_scores.csv")

# ectomycorrhizae ----
SLecto_filt <- SL_FT[which(SL_FT$ecto=="1"),]
SLecto_asv <- SLecto_filt[,-c(1509:1543)]
SLecto_asvT <- t(SLecto_asv)
colnames(SLecto_asvT) <- SLecto_asvT[1,]
tot <- rowSums(SLecto_asvT)
SLecto_asvT <- SLecto_asvT[tot > 0, ]
SLecto_rclr <- transform(x=SLecto_asvT, "rclr")
SLecto_rclr2 <- SLecto_rclr[,-1]
SLecto_rclr3 <- mutate_all(SLecto_rclr2, function(x) as.numeric(as.character(x)))
SLecto_rclr4 <- SLecto_rclr2[rowSums(SLecto_rclr3)>0, ]
#SLecto_rclr5 <- mutate_all(SLecto_rclr4, function(x) as.numeric(as.character(x)))
SLecto_ATCH <- vegdist(SLecto_rclr4, method="robust.aitchison")
#ectoSLmeta <- merged.meta2[((merged.meta2$ectorich)>0),] # 1513
#ectoSLmeta_filt <- ectoSLmeta[(ectoSLmeta$OTU_ID %in% rownames(SLecto_asvT)),]

#ectomod <- adonis2(SLecto_ATCH ~ subplot + compartment + season:year + eco.region + fire*habitat, by="terms", data=ectoSLmeta_filt)
#ectomod

#ectomod2 <- adonis2(SLecto_ATCH ~ compartment.x + iBut_MMT + Dew_mean + precip_7day + elevation + PercMoist + CN + pH, by="terms", data=ectoSLmeta_filt)
#ectomod2

#raremax <- min(rowSums(SLecto_asvT))
#ecto_r <- rrarefy(SLecto_asvT, raremax)
ecto_nmds <- metaMDS(SLecto_rclr4, distance="bray")
ecto_nmds
stressplot(ecto_nmds)
plot(ecto_nmds)

ecto.site <- as.data.frame(scores(ecto_nmds)$sites)
ecto.site$site <- rownames(ecto.site)
head(ecto.site)
write.csv(ecto.site, "ectosite_scores.csv")

ecto.species <- as.data.frame(scores(ecto_nmds)$species)
ecto.species$OTU_ID <- rownames(ecto.species)
head(ecto.species)
write.csv(ecto.species, "ectoOTU_scores.csv")

# pyrofiles ----
SLpyro_filt <- SL_FT[which(SL_FT$pyro=="1"),]
SLpyro_asv <- SLpyro_filt[,-c(1506:1543)]
rownames(SLpyro_asv) <- SLpyro_asv$OTU_ID
SLpyro_asv$OTU_ID <- NULL
SLpyro_asvT <- t(SLpyro_asv)
tot <- rowSums(SLpyro_asvT)
SLpyro_asvT <- SLpyro_asvT[tot > 0, ]
SLpyro_rclr <- transform(x=SLpyro_asvT, "rclr")
SLpyro_rclr2 <- mutate_all(SLpyro_rclr, function(x) as.numeric(as.character(x)))
SLpyro_rclr3 <- SLpyro_rclr2[rowSums(SLpyro_rclr2[, -1])>0, ]
SLpyro_ATCH <- vegdist(SLpyro_rclr2[,-1], method="robust.aitchison")
#SLpyrometa <- merged.meta2[((merged.meta2$pyrorich)>0),]
#SLpyrometa_filt <- SLpyrometa[(SLpyrometa$OTU_ID %in% rownames(SLpyro_asvT)),]

#raremax <- min(rowSums(SLpyro_asvT))
#pyro_r <- rrarefy(SLpyro_asvT, raremax)
pyro_nmds <- metaMDS(SLpyro_rclr2[,-1], distance="bray")
pyro_nmds
stressplot(pyro_nmds)
plot(pyro_nmds)

pyro.site <- as.data.frame(scores(pyro_nmds)$sites)
pyro.site$site <- rownames(pyro.site)
head(pyro.site)
write.csv(pyro.site, "pyrosite_scores.csv")

pyro.species <- as.data.frame(scores(pyro_nmds)$species)
pyro.species$OTU_ID <- rownames(pyro.species)
head(pyro.species)
write.csv(pyro.species, "pyroOTU_scores.csv")