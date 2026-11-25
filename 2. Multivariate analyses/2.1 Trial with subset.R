# Packages ----
library(vegan)
library(ggplot2)
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
set.seed(123)

# import data ----
setwd("output/")
# this has functional group data (richness, diversity, RA)
metadata1 <- read.csv("CliMush-SoilLitter-metaWfunctions_14Aug25.csv", header=T, sep=",")
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

# import SL_FT
#SL_FT <- read.csv("CliMush_SL_otu_functions1.csv", head=T, sep=",")
#View(SL_FT)
#SL_FT <- SL_FT[,-1]
## subset for tests
#SL_FT <- SL_FT[-c(20000:821300),]
#write.csv(SL_FT, "SL_FT_subset.csv")
SL_FT <- read.csv("SL_FT_subset.csv", head=T, sep=",")
head(SL_FT)
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


# pathogens ----
SLpath_filt <- SL_FT[which(SL_FT$ppath2=="1"),]
SLpath_asv <- SLpath_filt[,-c(1509:1548)]
rownames(SLpath_asv) <- SLpath_asv$OTU_ID
SLpath_asv$OTU_ID <- NULL
SLpath_asvT <- t(SLpath_asv)
#SLpath_asvT <- SLpath_asvT[rownames(SLpath_asvT) != "illumina_litter_2023.10_D16_BG_09", ]
tot <- rowSums(SLpath_asvT)
SLpath_asvT <- SLpath_asvT[tot > 0, ]
SLpathmeta <- metadata1[((metadata1$pprich1)>0),]
SLpathmeta_filt <- SLpathmeta[(SLpathmeta$OTU_ID %in% rownames(SLpath_asvT)),]

SLpath_rclr <- transform(x=SLpath_asvT, "rclr")
SLpath_rclr2 <- mutate_all(SLpath_rclr, function(x) as.numeric(as.character(x)))
SLpath_rclr3 <- SLpath_rclr2[rowSums(SLpath_rclr2[, -1])>0, ]
SLpath_rclr4 <- SLpath_rclr3[-c(1422:1423),-1]
SLpath_ATCH <- vegdist(SLpath_rclr4, method="robust.aitchison")

#write.csv(as.data.frame(SLpath_ATCH), "SLpath_ATCH.csv")
pathnames <- rownames(SLpath_rclr4)
pathmeta_perm <- metadata[(metadata$OTU_ID %in% pathnames),]
head(pathmeta_perm)


anova(betadisper(SLpath_ATCH, pathmeta_perm$compartment))
anova(betadisper(SLpath_ATCH, pathmeta_perm$habitat))
anova(betadisper(SLpath_ATCH, pathmeta_perm$ecoregion))
anova(betadisper(SLpath_ATCH, pathmeta_perm$fire))

# check
metadata$OTU_ID[!metadata$OTU_ID %in% pathnames]
pathnames[!pathnames %in% metadata$OTU_ID]

colSums(is.na(pathmeta_perm))
# Check for any NaN or Inf
sapply(pathmeta_perm, function(x) sum(is.nan(x) | is.infinite(x)))
# Check for empty strings in character/factor columns
sapply(pathmeta_perm, function(x) sum(trimws(as.character(x)) == ""))
pathmeta_perm_noNA <- pathmeta_perm[!is.na(pathmeta_perm$CN), ]
NAlist <- pathnames[!pathnames %in% pathmeta_perm_noNA$OTU_ID]
SLpath_asv_noNA <- SLpath_asv[!rownames(SLpath_asv) %in% NAlist,]
SLpath_asvT_noNA <- t(SLpath_asv_noNA)
tot <- rowSums(SLpath_asvT_noNA)
SLpath_asvT_noNA <- SLpath_asvT_noNA[tot > 0, ]
SLpath_rclr_noNA <- transform(x=SLpath_asvT_noNA, "rclr")
SLpath_rclr2_noNA <- mutate_all(SLpath_rclr_noNA, function(x) as.numeric(as.character(x)))
SLpath_rclr3_noNA <- SLpath_rclr2_noNA[rowSums(SLpath_rclr2_noNA[, -1])>0, ]
SLpath_rclr4_noNA <- SLpath_rclr3_noNA[-c(1422:1423),-1]
SLpath_ATCH_noNA <- vegdist(SLpath_rclr4_noNA, method="robust.aitchison")

adonis2(SLpath_ATCH_noNA ~ ecoregion + habitat + fire + subplot + compartment + year + season + tmean_Annual + ppt_mean_Annual + elevation + pH + CN + moisture, data=pathmeta_perm_noNA, by="terms")

#raremax <- min(rowSums(SLpath_asvT))
#path_r <- rrarefy(SLpath_asvT, raremax)
path_nmds <- metaMDS(SLpath_rclr4, distance="bray", maxit=80, trymax=80)
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
SL_FT$allsaps <- ifelse(SL_FT$primary_lifestyle=="dung_saprotroph" | SL_FT$primary_lifestyle=="litter_saprotroph" | SL_FT$primary_lifestyle=="soil_saprotroph" | SL_FT$primary_lifestyle=="wood_saprotroph", "1", "0")

SLsap_filt <- SL_FT[which(SL_FT$allsaps=="1"),]
SLsap_asv <- SLsap_filt[,-c(1509:1545)]
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

sapSLmeta <- metadata[((metadata$allsapsdiv)>0),] # 1513
sapSLmeta_filt <- sapSLmeta[(sapSLmeta$OTU_ID %in% rownames(SLsap_rclr4)),]

sapmod <- adonis2(SLsap_ATCH ~ compartment + year + season + fire + habitat + subplot + eco.region, by="terms", data=sapSLmeta_filt)
sapmod

anova(betadisper(SLsap_ATCH, sapSLmeta_filt$compartment))
anova(betadisper(SLsap_ATCH, sapSLmeta_filt$habitat))
anova(betadisper(SLsap_ATCH, sapSLmeta_filt$eco.region))
anova(betadisper(SLsap_ATCH, sapSLmeta_filt$fire))

#raremax <- min(rowSums(SLsap_asvT))
#sap_r <- rrarefy(SLsap_asvT, raremax)
sap_nmds <- metaMDS(SLsap_rclr4, distance="bray", maxit=80, trymax=80)
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
SLecto_asv <- SLecto_filt[,-c(1509:1545)]
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
ectoSLmeta <- metadata[((metadata$ectorich)>0),] # 1513
ectoSLmeta_filt <- ectoSLmeta[(ectoSLmeta$OTU_ID %in% rownames(SLecto_asvT)),]

ectomod <- adonis2(SLecto_ATCH ~ compartment + year + season + fire + habitat + subplot + eco.region, by="terms", data=ectoSLmeta_filt)
ectomod

anova(betadisper(SLecto_ATCH, ectoSLmeta_filt$compartment))
anova(betadisper(SLecto_ATCH, ectoSLmeta_filt$habitat))
anova(betadisper(SLecto_ATCH, ectoSLmeta_filt$eco.region))
anova(betadisper(SLecto_ATCH, ectoSLmeta_filt$fire))


#raremax <- min(rowSums(SLecto_asvT))
#ecto_r <- rrarefy(SLecto_asvT, raremax)
ecto_nmds <- metaMDS(SLecto_rclr4, distance="bray", maxit=80, trymax=80)
ecto_nmds
stressplot(ecto_nmds)
plot(ecto_nmds)

ecto.site <- as.data.frame(scores(ecto_nmds)$sites)
ecto.site$site <- rownames(ecto.site)
head(ecto.site)
write.csv(ecto.site, "ectosite_scores.csv")

# pyrofiles ----
SLpyro_filt <- SL_FT[which(SL_FT$pyro=="1"),]
SLpyro_asv <- SLpyro_filt[,-c(1509:1545)]
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
SLpyrometa <- metadata[((metadata$pyrorich)>0),]
SLpyrometa_filt <- SLpyrometa[(SLpyrometa$OTU_ID %in% rownames(SLpyro_asvT)),]

#raremax <- min(rowSums(SLpyro_asvT))
#pyro_r <- rrarefy(SLpyro_asvT, raremax)

pyromod <- adonis2(SLpyro_ATCH ~ compartment + year + season + fire + habitat + subplot + eco.region, by="terms", data=SLpyrometa_filt)
pyromod

anova(betadisper(SLpyro_ATCH, SLpyrometa_filt$compartment))
anova(betadisper(SLpyro_ATCH, SLpyrometa_filt$habitat))
anova(betadisper(SLpyro_ATCH, SLpyrometa_filt$eco.region))
anova(betadisper(SLpyro_ATCH, SLpyrometa_filt$fire))

pyro_nmds <- metaMDS(SLpyro_rclr4, distance="bray", maxit=80, trymax=80)
pyro_nmds
stressplot(pyro_nmds)
plot(pyro_nmds)

pyro.site <- as.data.frame(scores(pyro_nmds)$sites)
pyro.site$site <- rownames(pyro.site)
head(pyro.site)
write.csv(pyro.site, "pyrosite_scores.csv")


# old ----
# co-occurrence network
library(igraph)
library(Hmisc) 
library(Matrix)
# Read taxonomy file associated with OTU table into new object: otu_taxonomy.csv
tax<-read.csv("illumina_isl-2407.otus.taxonomy3.txt", header=T, sep="\t")
head(tax)

matchedtax <- tax[(tax$ASV_no %in% rownames(otu.table.filter)),]

otu.cor<-rcorr(as.matrix(otu.table.filter), type="spearman")
otu.pval <- forceSymmetric(otu.cor$P) # Self-correlation as NA
# Select only the taxa for the filtered OTUs by using rownames of otu.pval
sel.tax <-matchedtax[rownames(otu.pval),,drop=FALSE]
# Sanity check 
all.equal(rownames(sel.tax), rownames(otu.pval))

p.yes<-otu.pval<0.001
# Select the r values for p.yes
r.val=otu.cor$r # select all the correlation values 
p.yes.r<-r.val*p.yes # only select correlation values based on p-value criterion
# Select OTUs by level of correlation
p.yes.r<-abs(p.yes.r)>0.75 # output is logical vector
p.yes.rr<-p.yes.r*r.val # use logical vector for subscripting

# Create an adjacency matrix
adjm<-as.matrix(p.yes.rr)
# Add taxonomic information associated with adjacency matrix
colnames(adjm)<-as.vector(sel.tax$Family)
rownames(adjm)<-as.vector(sel.tax$Family)
# Create an adjacency matrix in igraph format 
net.grph=graph_from_adjacency_matrix(adjm,mode="undirected",weighted=TRUE,
                                     diag=FALSE)
# Calculate edge weight == level of correlation
edgew<-E(net.grph)$weight
# Identify isolated nodes
bad.vs<-V(net.grph)[degree(net.grph) == 0] 
# Remove isolated nodes
net.grph2 <-delete_vertices(net.grph, bad.vs)
# Plot the graph object
plot(net.grph,
     vertex.size=4,
     vertex.frame.color="black",
     edge.curved=F,
     edge.width=1.5,
     layout=layout.fruchterman.reingold,
     edge.color=ifelse(edgew<0,"red","blue"),
     vertex.label=NA,
     vertex.label.color="black",
     vertex.label.family="Times New Roman",
     vertex.label.font=2)

