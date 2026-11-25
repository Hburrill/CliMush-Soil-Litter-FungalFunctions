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
setwd("C:/Users/haley/University of Oregon Dropbox/Haley Burrill/CLIMUSH NSF-Neon/Analysis/SoilLitterFunctionPaper_HB/Analyses/output")
# this has functional group data (richness, diversity, RA)
metadata <- read.csv("CliMush-SoilLitter-metaWfunctions_6Aug25.csv", header=T, sep=",")
head(metadata)
metadata <- metadata[,-1]

# import SL_FT
#SL_FT <- read.csv("CliMush_SL_otu_functions1.csv", head=T, sep=",")
#head(SL_FT)
#SL_FT <- SL_FT[,-1]
## subset for tests
#SL_FT <- SL_FT[-c(20000:821300),]
#write.csv(SL_FT, "SL_FT_subset.csv")
SL_FT <- read.csv("SL_FT_subset.csv", head=T, sep=",")
head(SL_FT)
SL_FT <- SL_FT[,-1]

SL_FT2 <- mutate_all(SL_FT, function(x) as.numeric(as.character(x)))
SL_FT3 <- SL_FT2[,-c(1,1509:1543)]
otu.table.filter <- SL_FT3[,colSums(SL_FT3)>10]



# co-occurrence network
library(igraph)
library(Hmisc) 
library(Matrix)
# Read taxonomy file associated with OTU table into new object: otu_taxonomy.csv
tax<-read.csv(file.choose(),header=T, row.names = 1)

otu.cor<-rcorr(as.matrix(otu.table.filter), type="spearman")
otu.pval <- forceSymmetric(otu.cor$P) # Self-correlation as NA
# Select only the taxa for the filtered OTUs by using rownames of otu.pval
sel.tax <-tax[rownames(otu.pval),,drop=FALSE]
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
net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,
                         diag=FALSE)
# Calculate edge weight == level of correlation
edgew<-E(net.grph)$weight
# Identify isolated nodes
bad.vs<-V(net.grph)[degree(net.grph) == 0] 
# Remove isolated nodes
net.grph <-delete.vertices(net.grph, bad.vs)
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


# pathogens ----
SLpath_filt <- SL_FT[which(SL_FT$ppath2=="1",]
SLpath_asv <- SLpath_filt[,-c(1509:1548)]
rownames(SLpath_asv) <- SLpath_asv$OTU_ID
SLpath_asv$OTU_ID <- NULL
SLpath_asvT <- t(SLpath_asv)
tot <- rowSums(SLpath_asvT)
SLpath_asvT <- SLpath_asvT[tot > 0, ]
SLpathmeta <- merged.meta2[((merged.meta2$pprich1)>0),]
SLpathmeta_filt <- SLpathmeta[(SLpathmeta$OTU_ID %in% rownames(SLpath_asvT)),]

SLpath_rclr <- transform(x=SLpath_asvT, "rclr")
SLpath_rclr2 <- mutate_all(SLpath_rclr, function(x) as.numeric(as.character(x)))
SLpath_rclr3 <- SLpath_rclr2[rowSums(SLpath_rclr2[, -1])>0, ]
SLpath_rclr4 <- SLpath_rclr3[-c(1422:1423),-1]
SLpath_ATCH <- vegdist(SLpath_rclr4, method="robust.aitchison")
#write.csv(as.data.frame(SLpath_ATCH), "SLpath_ATCH.csv")

#pathmod <- adonis2(SLpath_ATCH ~ subplot.x + compartment.x + season.x:year.x + eco.region + fire.x*habitat.x, by="terms", data=SLpathmeta_filt)
#pathmod

#pathmod2 <- adonis2(SLpath_ATCH ~ compartment.x + iBut_MMT + Dew_mean + precip_7day + elevation + PercMoist + CN + pH, by="terms", data=SLpathmeta_filt)
#pathmod2

#raremax <- min(rowSums(SLpath_asvT))
#path_r <- rrarefy(SLpath_asvT, raremax)
path_nmds <- metaMDS(SLpath_rclr4, distance="bray", maxit=200, trymax=100)
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
SLsap_asv <- SLsap_filt[,-c(1509:1543)]
SLsap_asvT <- t(SLsap_asv)
colnames(SLsap_asvT) <- SLsap_asvT[1,]
SLsap_asvT <- SLsap_asvT[-1,]
tot <- rowSums(SLsap_asvT)
SLsap_asvT <- SLsap_asvT[tot > 0, ]
SLsap_rclr <- transform(x=SLsap_asvT, "rclr")
SLsap_rclr2 <- SLsap_rclr[,-1]
SLsap_rclr3 <- mutate_all(SLsap_rclr2, function(x) as.numeric(as.character(x)))
#SLsap_rclr4 <- SLsap_rclr3[rowSums(SLsap_rclr3>0), ]
SLsap_ATCH <- vegdist(SLsap_rclr3, method="robust.aitchison")

sapSLmeta <- merged.meta2[((merged.meta2$saprich)>0),] # 1513
sapSLmeta_filt <- sapSLmeta[(sapSLmeta$OTU_ID %in% rownames(SLsap_asvT)),]

#sapmod <- adonis2(SLsap_ATCH ~ subplot.x + compartment.x + season.x:year.x + eco.region + fire.x*habitat.x, by="terms", data=sapSLmeta_filt)
#sapmod

#sapmod2 <- adonis2(SLsap_ATCH ~ compartment.x + iBut_MMT + Dew_mean + precip_7day + elevation + PercMoist + CN + pH, by="terms", data=sapSLmeta_filt)
#sapmod2

#raremax <- min(rowSums(SLsap_asvT))
#sap_r <- rrarefy(SLsap_asvT, raremax)
sap_nmds <- metaMDS(SLsap_asvT, distance="bray")
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
ectoSLmeta <- merged.meta2[((merged.meta2$ectorich)>0),] # 1513
ectoSLmeta_filt <- ectoSLmeta[(ectoSLmeta$OTU_ID %in% rownames(SLecto_asvT)),]

#ectomod <- adonis2(SLecto_ATCH ~ subplot + compartment + season:year + eco.region + fire*habitat, by="terms", data=ectoSLmeta_filt)
#ectomod

#ectomod2 <- adonis2(SLecto_ATCH ~ compartment.x + iBut_MMT + Dew_mean + precip_7day + elevation + PercMoist + CN + pH, by="terms", data=ectoSLmeta_filt)
#ectomod2

#raremax <- min(rowSums(SLecto_asvT))
#ecto_r <- rrarefy(SLecto_asvT, raremax)
ecto_nmds <- metaMDS(SLecto_asvT, distance="bray")
ecto_nmds
stressplot(ecto_nmds)
plot(ecto_nmds)

ecto.scores <- as.data.frame(scores(ecto_nmds))
ecto.scores$OTU_ID <- rownames(ecto.scores)
head(ecto.scores)
write.csv(ecto.scores, "ectoscores.csv")

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
SLpyrometa <- merged.meta2[((merged.meta2$pyrorich)>0),]
SLpyrometa_filt <- SLpyrometa[(SLpyrometa$OTU_ID %in% rownames(SLpyro_asvT)),]

#raremax <- min(rowSums(SLpyro_asvT))
#pyro_r <- rrarefy(SLpyro_asvT, raremax)
pyro_nmds <- metaMDS(SLpyro_asvT, distance="bray")
pyro_nmds
stressplot(pyro_nmds)
plot(pyro_nmds)

pyro.scores <- as.data.frame(scores(pyro_nmds))
pyro.scores$OTU_ID <- rownames(pyro.scores)
pyro.scores$grp <- grp
write.csv(pyro.scores, "pyroscores.csv")


#pyromod <- adonis2(SLpyro_ATCH ~ subplot + compartment + season:year + eco.region + fire*habitat, by="terms", data=SLpyrometa_filt)
#pyromod

#pyromod2 <- adonis2(SLpyro_ATCH ~ compartment.x + iBut_MMT + Dew_mean + precip_7day + elevation + PercMoist + CN + pH, by="terms", data=SLpyrometa_filt)
#pyromod2

# indicator sp ----
# steps from https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html
## pathogens ----
#indic_path <- cbind(merged.meta2, SLpath_asvT[-(1518:1520),])
#indic_path2 <- indic_path[-1517,]
#indic_path3 <- indic_path2[indic_path2$compartment!="soil-litter",]
#pathabund = indic_path3[,31:ncol(indic_path3)]
#pathburn = indic_path3$fire
#pathinv = multipatt(as.matrix(sapply(pathabund, as.numeric)), pathburn, func = "r.g", control = how(nperm=9999))
#summary(pathinv)
#sink(file = "pathinv.txt")
#summary(pathinv)
#sink(file = NULL)

#pathveg = indic_path3$habitat
#pathinv_veg = multipatt(as.matrix(sapply(pathabund, as.numeric)), pathveg, func = "r.g", control = how(nperm=999))
#summary(pathinv_veg)
#sink(file = "pathinv_veg.txt")
#summary(pathinv_veg)
#sink(file = NULL)

#sub_indic_path = data.frame(Burn = indic_path3$fire, ASV23227 = indic_path3$ASV23227, ASV12919 = indic_path3$ASV12919, ASV11629 = indic_path3$ASV11629, ASV7060 = indic_path3$ASV7060)
#sub_indic_path_m = melt(sub_indic_path, id = "Burn")

#ggplot(sub_indic_path_m, aes(x = variable, y = value, fill = Burn)) + geom_boxplot(colour = "black", position = position_dodge(0.5)) +  geom_vline(xintercept = c(1.5,2.5,3.5), colour = "grey85", size = 1.2) + theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10, face = "bold"), legend.position = "right", axis.text.x = element_text(face = "bold",colour = "black", size = 12), axis.text.y = element_text(face = "bold", size = 12, colour = "black"), axis.title.y = element_text(face = "bold", size = 14, colour = "black"),     panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),       legend.key=element_blank()) + labs(x= "", y = "Relative Abundance (%)", fill = "Fire") + scale_fill_manual(values = c("steelblue2", "steelblue4"))

## saprotrophs ----

## ecto ----


## pyro ----

