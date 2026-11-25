# Packages ----
library(vegan)
library(ggplot2)
library(dplyr)
library(base)
library(indicspecies)
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

rownames(SL_FT) <- paste0("ASV", seq_len(nrow(SL_FT)))

SL_FT2 <- mutate_all(SL_FT, function(x) as.numeric(as.character(x)))
SL_FT3 <- SL_FT2[,-c(1,1509:1544)]
otu.table.filter <- SL_FT3[,colSums(SL_FT3)>10]

# allfung ----
SL_FT3T <- t(SL_FT3)
tot <- rowSums(SL_FT3T)
SL_FT3T <- SL_FT3T[tot > 0, ]
md_all <- metadata[(metadata$OTU_ID %in% rownames(SL_FT3T)),]

SL_rclr <- transform(x=SL_FT3T, "rclr")
SL_rclr2 <- mutate_all(SL_rclr, function(x) as.numeric(as.character(x)))
SL_rclr3 <- SL_rclr2[rowSums(SL_rclr2[, -1])>0, ]
SL_rclr4 <- SL_rclr3[,-1]
SL_ATCH <- vegdist(SL_rclr4, method="robust.aitchison")

adonis2(SL_ATCH ~ ecoregion + habitat + fire + subplot + compartment + year + season + tmean_Annual + ppt_mean_Annual, data=md_all, by="terms")


# indicator sp ----
# steps from https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html

## pathogens ----
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
writeLines(rownames(SLpath_rclr4), "SLpath_rclr4_names2.txt")

pathnames <- read.table("SLpath_rclr4_names2.txt", header=F)
pathmeta_perm <- metadata[(metadata$OTU_ID %in% pathnames$V1),]

indic_path <- cbind(pathmeta_perm, SLpath_asvT)

# columns 104:522 names == OTU_ID
# check lengths first
length(rownames(SLpath_asv))
length(104:5200)
# rename columns
colnames(indic_path)[104:5200] <- rownames(SLpath_asvT)

#indic_path2 <- indic_path[-1517,]
#indic_path3 <- indic_path2[indic_path2$compartment!="soil-litter",]
pathabund = indic_path[,104:ncol(indic_path)]
pathburn = indic_path$fire
pathinv = multipatt(as.matrix(sapply(pathabund, as.numeric)), pathburn, func = "r.g", control = how(nperm=999))

pathinv_df <- as.data.frame(pathinv$sign)
pathinv_df$OTU_ID <- rownames(SLpath_filt)
SLpath_info <- SLpath_filt[,c(1509:1518)]
SLpath_info$OTU_ID <- rownames(SLpath_info)
pathinv_df2 <- left_join(pathinv_df, SLpath_info, by="OTU_ID")
write.csv(pathinv_df2, "pathfire.csv")

patheco = indic_path$ecoregion
pathinv2 = multipatt(as.matrix(sapply(pathabund, as.numeric)), patheco, func = "r.g", control = how(nperm=999))
pathinv2_df <- as.data.frame(pathinv2$sign)
pathinv2_df$OTU_ID <- rownames(SLpath_filt)
SLpath2_info <- SLpath_filt[,c(1509:1518)]
SLpath2_info$OTU_ID <- rownames(SLpath2_info)
pathinv3_df <- left_join(pathinv2_df, SLpath2_info, by="OTU_ID")
write.csv(pathinv3_df, "patheco.csv")

pathhab = indic_path$habitat
pathinv3 = multipatt(as.matrix(sapply(pathabund, as.numeric)), pathhab, func = "r.g", control = how(nperm=999))
pathinv3_df <- as.data.frame(pathinv3$sign)
pathinv3_df$OTU_ID <- rownames(SLpath_filt)
SLpath3_info <- SLpath_filt[,c(1509:1518)]
SLpath3_info$OTU_ID <- rownames(SLpath3_info)
pathinv4_df <- left_join(pathinv3_df, SLpath3_info, by="OTU_ID")
write.csv(pathinv4_df, "pathhab.csv")


pathcomp = indic_path$compartment
pathinv4 = multipatt(as.matrix(sapply(pathabund, as.numeric)), pathcomp, func = "r.g", control = how(nperm=999))
pathinv4_df <- as.data.frame(pathinv4$sign)
pathinv4_df$OTU_ID <- rownames(SLpath_filt)
SLpath4_info <- SLpath_filt[,c(1509:1518)]
SLpath4_info$OTU_ID <- rownames(SLpath4_info)
pathinv5_df <- left_join(pathinv4_df, SLpath_info, by="OTU_ID")
write.csv(pathinv5_df, "pathcomp.csv")


pathmat = indic_path$tmean_Annual
pathinv5 = multipatt(as.matrix(sapply(pathabund, as.numeric)), pathmat, func = "r.g", control = how(nperm=999))
pathinv5_df <- as.data.frame(pathinv5$sign)
pathinv5_df$OTU_ID <- rownames(SLpath_filt)
SLpath5_info <- SLpath_filt[,c(1509:1518)]
SLpath5_info$OTU_ID <- rownames(SLpath5_info)
pathinv6_df <- left_join(pathinv5_df, SLpath_info, by="OTU_ID")
write.csv(pathinv6_df, "pathmat.csv")

pathmap = indic_path$ppt_mean_Annual
pathinv6 = multipatt(as.matrix(sapply(pathabund, as.numeric)), pathmap, func = "r.g", control = how(nperm=999))
pathinv6_df <- as.data.frame(pathinv6$sign)
pathinv6_df$OTU_ID <- rownames(SLpath_filt)
SLpath6_info <- SLpath_filt[,c(1509:1518)]
SLpath6_info$OTU_ID <- rownames(SLpath6_info)
pathinv7_df <- left_join(pathinv6_df, SLpath_info, by="OTU_ID")
write.csv(pathinv7_df, "pathmap.csv")

## saprotrophs ----
SL_FT$allsaps <- ifelse(SL_FT$primary_lifestyle=="dung_saprotroph" | SL_FT$primary_lifestyle=="litter_saprotroph" | SL_FT$primary_lifestyle=="soil_saprotroph" | SL_FT$primary_lifestyle=="wood_saprotroph", "1", "0")

SLsap_filt <- SL_FT[which(SL_FT$allsaps=="1"),]
SLsap_asv <- SLsap_filt[,-c(1509:1548)]
SLsap_asvT <- t(SLsap_asv)
colnames(SLsap_asvT) <- SLsap_asvT[1,]
SLsap_asvT <- SLsap_asvT[-1,]
tot <- rowSums(SLsap_asvT)
SLsap_asvT <- SLsap_asvT[tot > 0, ]
SLsap_rclr <- transform(x=SLsap_asvT, "rclr")
SLsap_rclr <- SLsap_rclr[,-1]
SLsap_rclr2 <- mutate_all(SLsap_rclr, function(x) as.numeric(as.character(x)))
SLsap_rclr3 <- SLsap_rclr2[rowSums(SLsap_rclr2[, -1])>0, ]
SLsap_rclr4 <- SLsap_rclr3[-c(1422:1423),-1]
SLsap_ATCH <- vegdist(SLsap_rclr4, method="robust.aitchison")

writeLines(rownames(SLsap_rclr4), "SLsap_rclr4_names2.txt")
sapnames <- read.table("SLsap_rclr4_names2.txt", header=F)
sapmeta_perm <- metadata[(metadata$OTU_ID %in% sapnames$V1),]
SLsap_asvT2 <- SLsap_asvT[(rownames(SLsap_asvT) %in% sapmeta_perm$OTU_ID),]

indic_sap <- cbind(sapmeta_perm, SLsap_asvT2)
sapabund = indic_sap[,104:ncol(indic_sap)]
sapburn = indic_sap$fire
sapinv = multipatt(as.matrix(sapply(sapabund, as.numeric)), sapburn, func = "r.g", control = how(nperm=999))
sapinv_df <- as.data.frame(sapinv$sign)
sapinv_df$OTU_ID <- rownames(SLsap_filt)
SLsap_info <- SLsap_filt[,c(1509:1518)]
SLsap_info$OTU_ID <- rownames(SLsap_info)
sapinv_df2 <- left_join(sapinv_df, SLsap_info, by="OTU_ID")
write.csv(sapinv_df2, "sapfire.csv")

sapeco = indic_sap$ecoregion
sapinv2 = multipatt(as.matrix(sapply(sapabund, as.numeric)), sapeco, func = "r.g", control = how(nperm=999))
sapinv2_df <- as.data.frame(sapinv2$sign)
sapinv2_df$OTU_ID <- rownames(SLsap_filt)
SLsap2_info <- SLsap_filt[,c(1509:1518)]
SLsap2_info$OTU_ID <- rownames(SLsap2_info)
sapinv3_df <- left_join(sapinv2_df, SLsap2_info, by="OTU_ID")
write.csv(sapinv3_df, "sapeco.csv")

saphab = indic_sap$habitat
sapinv3 = multipatt(as.matrix(sapply(sapabund, as.numeric)), saphab, func = "r.g", control = how(nperm=999))
sapinv3_df <- as.data.frame(sapinv3$sign)
sapinv3_df$OTU_ID <- rownames(SLsap_filt)
SLsap3_info <- SLsap_filt[,c(1509:1518)]
SLsap3_info$OTU_ID <- rownames(SLsap3_info)
sapinv4_df <- left_join(sapinv3_df, SLsap3_info, by="OTU_ID")
write.csv(sapinv4_df, "saphab.csv")


sapcomp = indic_sap$compartment
sapinv4 = multipatt(as.matrix(sapply(sapabund, as.numeric)), sapcomp, func = "r.g", control = how(nperm=9999))
sapinv4_df <- as.data.frame(sapinv4$sign)
sapinv4_df$OTU_ID <- rownames(SLsap_filt)
SLsap4_info <- SLsap_filt[,c(1509:1518)]
SLsap4_info$OTU_ID <- rownames(SLsap4_info)
sapinv5_df <- left_join(sapinv4_df, SLsap_info, by="OTU_ID")
write.csv(sapinv5_df, "sapcomp.csv")

sapmat = indic_sap$tmean_Annual
sapinv5 = multipatt(as.matrix(sapply(sapabund, as.numeric)), sapmat, func = "r.g", control = how(nperm=9999))
sapinv5_df <- as.data.frame(sapinv5$sign)
sapinv5_df$OTU_ID <- rownames(SLsap_filt)
SLsap5_info <- SLsap_filt[,c(1509:1518)]
SLsap5_info$OTU_ID <- rownames(SLsap5_info)
sapinv6_df <- left_join(sapinv5_df, SLsap_info, by="OTU_ID")
write.csv(sapinv6_df, "sapmat.csv")

sapmap = indic_sap$ppt_mean_Annual
sapinv6 = multipatt(as.matrix(sapply(sapabund, as.numeric)), sapmap, func = "r.g", control = how(nperm=9999))
sapinv6_df <- as.data.frame(sapinv6$sign)
sapinv6_df$OTU_ID <- rownames(SLsap_filt)
SLsap6_info <- SLsap_filt[,c(1509:1518)]
SLsap6_info$OTU_ID <- rownames(SLsap6_info)
sapinv7_df <- left_join(sapinv6_df, SLsap_info, by="OTU_ID")
write.csv(sapinv7_df, "sapmap.csv")


## ecto ----
SLecto_filt <- SL_FT[which(SL_FT$ecto=="1"),]
SLecto_asv <- SLecto_filt[,-c(1509:1548)]
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

writeLines(rownames(SLecto_rclr4), "SLecto_rclr4_names2.txt")
ectonames <- read.table("SLecto_rclr4_names2.txt", header=F)
ectometa_perm <- metadata[(metadata$OTU_ID %in% ectonames$V1),]
SLecto_asvT2 <- SLecto_asvT[(rownames(SLecto_asvT) %in% ectometa_perm$OTU_ID),]

indic_ecto <- cbind(ectometa_perm, SLecto_asvT2)
ectoabund = indic_ecto[,104:ncol(indic_ecto)]
ectoburn = indic_ecto$fire
ectoinv = multipatt(as.matrix(sapply(ectoabund, as.numeric)), ectoburn, func = "r.g", control = how(nperm=999))
ectoinv_df <- as.data.frame(ectoinv$sign)
ectoinv_df$OTU_ID <- rownames(SLecto_filt)
SLecto_info <- SLecto_filt[,c(1509:1518)]
SLecto_info$OTU_ID <- rownames(SLecto_info)
ectoinv_df2 <- left_join(ectoinv_df, SLecto_info, by="OTU_ID")
write.csv(ectoinv_df2, "ectofire.csv")


ectoeco = indic_ecto$ecoregion
ectoinv2 = multipatt(as.matrix(sapply(ectoabund, as.numeric)), ectoeco, func = "r.g", control = how(nperm=999))
ectoinv2_df <- as.data.frame(ectoinv2$sign)
ectoinv2_df$OTU_ID <- rownames(SLecto_filt)
SLecto2_info <- SLecto_filt[,c(1509:1518)]
SLecto2_info$OTU_ID <- rownames(SLecto2_info)
ectoinv3_df <- left_join(ectoinv2_df, SLecto2_info, by="OTU_ID")
write.csv(ectoinv3_df, "ectoeco.csv")

ectohab = indic_ecto$habitat
ectoinv3 = multipatt(as.matrix(sapply(ectoabund, as.numeric)), ectohab, func = "r.g", control = how(nperm=999))
ectoinv3_df <- as.data.frame(ectoinv3$sign)
ectoinv3_df$OTU_ID <- rownames(SLecto_filt)
SLecto3_info <- SLecto_filt[,c(1509:1518)]
SLecto3_info$OTU_ID <- rownames(SLecto3_info)
ectoinv4_df <- left_join(ectoinv3_df, SLecto3_info, by="OTU_ID")
write.csv(ectoinv4_df, "ectohab.csv")

ectocomp = indic_ecto$compartment
ectoinv4 = multipatt(as.matrix(sapply(ectoabund, as.numeric)), ectocomp, func = "r.g", control = how(nperm=999))
ectoinv4_df <- as.data.frame(ectoinv4$sign)
ectoinv4_df$OTU_ID <- rownames(SLecto_filt)
SLecto4_info <- SLecto_filt[,c(1509:1518)]
SLecto4_info$OTU_ID <- rownames(SLecto4_info)
ectoinv5_df <- left_join(ectoinv4_df, SLecto_info, by="OTU_ID")
write.csv(ectoinv5_df, "ectocomp.csv")

ectomat = indic_ecto$tmean_Annual
ectoinv5 = multipatt(as.matrix(sapply(ectoabund, as.numeric)), ectomat, func = "r.g", control = how(nperm=999))
ectoinv5_df <- as.data.frame(ectoinv5$sign)
ectoinv5_df$OTU_ID <- rownames(SLecto_filt)
SLecto5_info <- SLecto_filt[,c(1509:1518)]
SLecto5_info$OTU_ID <- rownames(SLecto5_info)
ectoinv6_df <- left_join(ectoinv5_df, SLecto_info, by="OTU_ID")
write.csv(ectoinv6_df, "ectomat.csv")

ectomap = indic_ecto$ppt_mean_Annual
ectoinv6 = multipatt(as.matrix(sapply(ectoabund, as.numeric)), ectomap, func = "r.g", control = how(nperm=999))
ectoinv6_df <- as.data.frame(ectoinv6$sign)
ectoinv6_df$OTU_ID <- rownames(SLecto_filt)
SLecto6_info <- SLecto_filt[,c(1509:1518)]
SLecto6_info$OTU_ID <- rownames(SLecto6_info)
ectoinv7_df <- left_join(ectoinv6_df, SLecto_info, by="OTU_ID")
write.csv(ectoinv7_df, "ectomap.csv")

## pyro ----
SLpyro_filt <- SL_FT[which(SL_FT$pyro=="1"),]
SLpyro_asv <- SLpyro_filt[,-c(1509:1548)]
rownames(SLpyro_asv) <- SLpyro_asv$OTU_ID
SLpyro_asv$OTU_ID <- NULL
SLpyro_asvT <- t(SLpyro_asv)
tot <- rowSums(SLpyro_asvT)
SLpyro_asvT <- SLpyro_asvT[tot > 0, ]
SLpyro_rclr <- transform(x=SLpyro_asvT, "rclr")
SLpyro_rclr2 <- SLpyro_rclr[,-1]
SLpyro_rclr3 <- mutate_all(SLpyro_rclr2, function(x) as.numeric(as.character(x)))
SLpyro_rclr4 <- SLpyro_rclr2[rowSums(SLpyro_rclr3)>0, ]
SLpyro_ATCH <- vegdist(SLpyro_rclr2[,-1], method="robust.aitchison")

writeLines(rownames(SLpyro_rclr4), "SLpyro_rclr4_names2.txt")
pyronames <- read.table("SLpyro_rclr4_names2.txt", header=F)
pyrometa_perm <- metadata[(metadata$OTU_ID %in% pyronames$V1),]
SLpyro_asvT2 <- SLpyro_asvT[(rownames(SLpyro_asvT) %in% pyrometa_perm$OTU_ID),]

indic_pyro <- cbind(pyrometa_perm, SLpyro_asvT2)
pyroabund = indic_pyro[,104:ncol(indic_pyro)]
pyroburn = indic_pyro$fire
pyroinv = multipatt(as.matrix(sapply(pyroabund, as.numeric)), pyroburn, func = "r.g", control = how(nperm=999))
pyroinv_df <- as.data.frame(pyroinv$sign)
pyroinv_df$OTU_ID <- rownames(SLpyro_filt)
SLpyro_info <- SLpyro_filt[,c(1509:1518)]
SLpyro_info$OTU_ID <- rownames(SLpyro_info)
pyroinv_df2 <- left_join(pyroinv_df, SLpyro_info, by="OTU_ID")
write.csv(pyroinv_df2, "pyrofire.csv")

pyroeco = indic_pyro$ecoregion
pyroinv2 = multipatt(as.matrix(sapply(pyroabund, as.numeric)), pyroeco, func = "r.g", control = how(nperm=999))
pyroinv2_df <- as.data.frame(pyroinv2$sign)
pyroinv2_df$OTU_ID <- rownames(SLpyro_filt)
SLpyro2_info <- SLpyro_filt[,c(1509:1518)]
SLpyro2_info$OTU_ID <- rownames(SLpyro2_info)
pyroinv3_df <- left_join(pyroinv2_df, SLpyro2_info, by="OTU_ID")
write.csv(pyroinv3_df, "pyroeco.csv")

pyrohab = indic_pyro$habitat
pyroinv3 = multipatt(as.matrix(sapply(pyroabund, as.numeric)), pyrohab, func = "r.g", control = how(nperm=999))
pyroinv3_df <- as.data.frame(pyroinv3$sign)
pyroinv3_df$OTU_ID <- rownames(SLpyro_filt)
SLpyro3_info <- SLpyro_filt[,c(1509:1518)]
SLpyro3_info$OTU_ID <- rownames(SLpyro3_info)
pyroinv4_df <- left_join(pyroinv3_df, SLpyro3_info, by="OTU_ID")
write.csv(pyroinv4_df, "pyrohab.csv")

pyrocomp = indic_pyro$compartment
pyroinv4 = multipatt(as.matrix(sapply(pyroabund, as.numeric)), pyrocomp, func = "r.g", control = how(nperm=999))
pyroinv4_df <- as.data.frame(pyroinv4$sign)
pyroinv4_df$OTU_ID <- rownames(SLpyro_filt)
SLpyro4_info <- SLpyro_filt[,c(1509:1518)]
SLpyro4_info$OTU_ID <- rownames(SLpyro4_info)
pyroinv5_df <- left_join(pyroinv4_df, SLpyro_info, by="OTU_ID")
write.csv(pyroinv5_df, "pyrocomp.csv")


pyromat = indic_pyro$tmean_Annual
pyroinv5 = multipatt(as.matrix(sapply(pyroabund, as.numeric)), pyromat, func = "r.g", control = how(nperm=999))
pyroinv5_df <- as.data.frame(pyroinv5$sign)
pyroinv5_df$OTU_ID <- rownames(SLpyro_filt)
SLpyro5_info <- SLpyro_filt[,c(1509:1518)]
SLpyro5_info$OTU_ID <- rownames(SLpyro5_info)
pyroinv6_df <- left_join(pyroinv5_df, SLpyro_info, by="OTU_ID")
write.csv(pyroinv6_df, "pyromat.csv")


pyromap = indic_pyro$ppt_mean_Annual
pyroinv6 = multipatt(as.matrix(sapply(pyroabund, as.numeric)), pyromap, func = "r.g", control = how(nperm=999))
pyroinv6_df <- as.data.frame(pyroinv6$sign)
pyroinv6_df$OTU_ID <- rownames(SLpyro_filt)
SLpyro6_info <- SLpyro_filt[,c(1509:1518)]
SLpyro6_info$OTU_ID <- rownames(SLpyro6_info)
pyroinv7_df <- left_join(pyroinv6_df, SLpyro_info, by="OTU_ID")
write.csv(pyroinv7_df, "pyromap.csv")
