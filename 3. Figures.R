# Figures for SL ----

library(vegan)
library(phyloseq)
library(ggplot2)
library(microViz)
library(AICcPermanova)
library(dplyr)
library(ggordiplots)
library(base)
library(glmulti)
library(lme4)
library(Matrix)
library(CCA)
library(reshape2)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,  conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop, .fun = function(xx, col) {
    c(N = length2(xx[[col]], na.rm=na.rm), mean = mean   (xx[[col]], na.rm=na.rm), sd   = sd     (xx[[col]], na.rm=na.rm))}, measurevar)
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}


## pathogens ----
### univariate ----
# import SLpathmeta_noSL
meanpathdiv <- summarySE(SLpathmeta_noSL, measurevar="pathdiv", groupvars="habitat")
ggplot(meanpathdiv, aes(x=habitat, y=pathdiv, color=habitat)) + geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=pathdiv-se,ymax=pathdiv+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="Fungal pathogen diversity", x="Veg type", y = "Shannon-Index Diversity")+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=09), legend.position = "bottom") + scale_color_manual(name="habitat", values=c("darkgreen","brown","orange")) # lowest in c, g & o same

meanpathdivfire <- summarySE(SLpathmeta_noSL, measurevar="pathdiv", groupvars="fire")
ggplot(meanpathdivfire, aes(x=fire, y=pathdiv, color=fire)) + geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=pathdiv-se,ymax=pathdiv+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="Fungal pathogen diversity", x="Burn", y = "Shannon-Index Diversity")+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=09), legend.position = "bottom") + scale_color_manual(name="fire", values=c("darkred","gold")) # pathogen diversity higher in burned, in general

meanpathcomp <- summarySE(SLpathmeta_noSL, measurevar="pathdiv", groupvars="compartment")
ggplot(meanpathcomp, aes(x=compartment, y=pathdiv, color=compartment)) +  geom_point(size =4) + geom_errorbar(aes(ymin=pathdiv-se,ymax=pathdiv+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="Fungal pathogen diversity", x="Compartment", y = "Shannon-Index Diversity")+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=09), legend.position = "bottom") + facet_grid(scales = "free_x", space = "free_x") + scale_color_manual(name="compartment", values=c("turquoise","darkviolet"))

### indicator sp ----
# import sub_indic_path_m

ggplot(sub_indic_path_m, aes(x = variable, y = value, fill = Burn)) + geom_boxplot(colour = "black", position = position_dodge(0.5)) +  geom_vline(xintercept = c(1.5,2.5,3.5), colour = "grey85", size = 1.2) + theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10, face = "bold"), legend.position = "right", axis.text.x = element_text(face = "bold",colour = "black", size = 12), axis.text.y = element_text(face = "bold", size = 12, colour = "black"), axis.title.y = element_text(face = "bold", size = 14, colour = "black"),     panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),       legend.key=element_blank()) + labs(x= "", y = "Relative Abundance (%)", fill = "Fire") + 
  scale_fill_manual(values = c("steelblue2", "steelblue4"))

### composition ----
# import SLpathmeta
SLpathmeta <- read.csv("Soil-Litter-pathogens-metadata.csv")

# remove soil-litter
SLpathmeta2 <- SLpathmeta[-(SLpathmeta$compartment=="soil-litter"),]
SLpathmeta2 <- SLpathmeta[-c(743:753),]


ggplot(SLpathmeta2, aes(x = pathpc1, y = pathpc2, col = habitat)) + geom_point(alpha = 0.3) + stat_ellipse(geom="polygon", alpha = 0, show.legend = T, level = 0.95) + theme_minimal() + labs(title="Fungal pathogen composition",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "top") + facet_wrap(~compartment, labeller = labeller(season = c("litter" = "Litter", "soil" = "Soil"))) + theme_classic() + scale_color_manual(name="habitat", values=c("darkgreen","darkorange","brown"))

ggplot(SLpathmeta2, aes(x = pathpc1, y = pathpc2, col = eco.region)) + geom_point(alpha = 0.3) + stat_ellipse(geom="polygon", alpha = 0, show.legend = T, level = 0.95) + theme_minimal() + labs(title="Fungal pathogen composition",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "top")

nofirepath <- SLpathmeta2[(SLpathmeta2$eco.region!="SROC" & SLpathmeta2$eco.region!="NE"),]
# nofirepath <- SLpathmeta2[-c(116:132,  287:316,  535:544,  570:579,  644:653,  683:692,  874:891, 1028:1043, 1222:1227, 1253:1262, 1326:1335, 1366:1372),]

# fire pc2 & 4
ggplot(nofirepath, aes(x = pathpc1, y = pathpc2, col = fire)) + geom_point(alpha = 0.4) + stat_ellipse(geom="polygon", aes(fill = fire), alpha = 0.05, show.legend = T, level = 0.95) + theme_minimal() + labs(title="Fungal pathogen composition",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "top") + facet_wrap(~eco.region) + theme_classic() + scale_color_manual(name="fire", values=c("maroon","darkturquoise"))

ggplot(nofirepath, aes(x = pathpc1, y = pathpc2, col = fire)) + geom_point(alpha = 0.4) + stat_ellipse(geom="polygon", aes(fill = fire), alpha = 0.05, show.legend = T, level = 0.95) + theme_minimal() + labs(title="Fungal pathogen composition",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "top") + facet_wrap(~habitat) + theme_classic() + scale_color_manual(name="fire", values=c("maroon","darkturquoise"))

## saprotrophs ----
### univariate ----
# import sapSLmeta_noSL

meansapdiv <- summarySE(sapSLmeta_noSL, measurevar="sapdiv", groupvars=c("compartment","habitat"))
ggplot(meansapdiv, aes(x=compartment, y=sapdiv, color=habitat)) + facet_wrap(~habitat) +  geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=sapdiv-se,ymax=sapdiv+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="Fungal saprotroph diversity", x="Eco region", y = "Shannon-Index Diversity")+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=09), legend.position = "bottom") # greatest difference in oaks, litter saprotroph diversity higher than soil

meansaprichcomp <- summarySE(sapSLmeta_noSL, measurevar="saprich", groupvars="compartment")
ggplot(meansaprichcomp, aes(x=compartment, y=saprich, color=compartment)) + geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=saprich-se,ymax=saprich+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="Fungal saprotroph richness", x="Compartment", y = "Richness")+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=09), legend.position = "bottom") + scale_color_manual(name="compartment", values=c("turquoise","darkviolet"))

meansaprichhab <- summarySE(sapSLmeta_noSL, measurevar="saprich", groupvars="habitat")
ggplot(meansaprichhab, aes(x=habitat, y=saprich, color=habitat)) + geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=saprich-se,ymax=saprich+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="Fungal saprotroph richness", x="Vegetation type", y = "Richness")+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=09), legend.position = "bottom") + scale_color_manual(name="habitat", values=c("darkgreen","brown","orange"))

### indicator sp ----
# import 


### composition ----
# import sapSLmeta2
sapSLmeta2 <- read.csv("Soil-Litter-saprotrophs-metadata.csv")
# remove soil-litter
sapSLmeta3 <- sapSLmeta2[-c(720:730),]
write.csv(sapSLmeta3, "Soil-Litter-saprotrophs-metadata.csv")

sap1 <- ggplot(sapSLmeta3, aes(x = sappc1, y = sappc2, col = habitat)) + geom_point(alpha = 0.4) + stat_ellipse(geom="polygon", aes(fill = habitat), alpha = 0.05, show.legend = T, level = 0.95) + theme_classic() + theme(plot.title = element_text(size=13), legend.position = "top") + facet_wrap(~compartment, labeller = labeller(season = c("litter" = "Litter", "soil" = "Soil"))) + labs(x="PC1", y="PC2")
sap1

sap2 <- ggplot(sapSLmeta3, aes(x = sappc1, y = sappc2, col = habitat)) + geom_point(alpha = 0.4) + stat_ellipse(geom="polygon", aes(fill = habitat), alpha = 0.05, show.legend = T, level = 0.95) + theme_classic() + theme(plot.title = element_text(size=13), legend.position = "top") + facet_wrap(~season) + labs(x="PC1", y="PC2")
sap2

ggplot(sapSLmeta3, aes(x = sappc1, y = sappc2, col = eco.region)) + geom_point(alpha = 0.3) + stat_ellipse(geom="polygon", alpha = 0, show.legend = T, level = 0.95) + theme_minimal() + labs(title="Saprotroph fungal composition",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "top") + facet_wrap(~season) + theme_classic()

nofiresap <- sapSLmeta3[-c(sapSLmeta3$eco.region=="SROC"| sapSLmeta3$eco.region=="NE"),]
which(sapSLmeta3$eco.region=="SROC"| sapSLmeta3$eco.region=="NE")
nofiresap <- sapSLmeta3[-c(116:133, 291:320, 513:522, 548:557, 622:631, 661:670, 852:869, 1031:1058, 1267:1276, 1302:1311, 1377:1386, 1417:1426),]

ggplot(nofiresap, aes(x = sappc1, y = sappc2, col = fire)) + geom_point(alpha = 0.4) + stat_ellipse(geom="polygon", aes(fill = fire), alpha = 0.05, show.legend = T, level = 0.95) + theme_minimal() + labs(title="Saprotroph composition",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "top") + facet_wrap(~habitat) + theme_classic() + scale_color_manual(name="fire", values=c("maroon","darkturquoise"))


## ectomycorrhizae ----
### univariate ----
# import ectoSLmeta_noSL

meanectodiv <- summarySE(ectoSLmeta_noSL, measurevar="ectodiv", groupvars="habitat")
ggplot(meanectodiv, aes(x=habitat, y=ectodiv, color=habitat)) + geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=ectodiv-se,ymax=ectodiv+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="", x="Veg type", y = "Shannon-Index Diversity")+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=09), legend.position = "none") + scale_color_manual(name="habitat", values=c("darkgreen","brown","orange"))

meanectofire <- summarySE(ectoSLmeta_noSL, measurevar="ectodiv", groupvars="fire")
ggplot(meanectofire, aes(x=fire, y=ectodiv, color=fire)) +  geom_point(size =4) + geom_errorbar(aes(ymin=ectodiv-se,ymax=ectodiv+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="Fungal ectomycorrhizal diversity", x="Fire", y = "Shannon-Index Diversity")+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=09), legend.position = "bottom") + facet_grid(scales = "free_x", space = "free_x") + scale_color_manual(name="fire", values=c("darkred","gold"))

meanectofirehab <- summarySE(ectoSLmeta_noSL, measurevar="ectodiv", groupvars=c("fire","habitat"))
ggplot(meanectofirehab, aes(x=habitat, y=ectodiv, color=fire)) +  geom_point(size =4) + geom_errorbar(aes(ymin=ectodiv-se,ymax=ectodiv+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="Fungal ectomycorrhizal diversity", x="Fire", y = "Shannon-Index Diversity")+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=09), legend.position = "bottom") + facet_grid(scales = "free_x", space = "free_x") + scale_color_manual(name="fire", values=c("darkred","gold"))

meanectorichcomp <- summarySE(ectoSLmeta_noSL, measurevar="ectorich", groupvars="compartment")
ggplot(meanectorichcomp, aes(x=compartment, y=ectorich, color=compartment)) + geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=ectorich-se,ymax=ectorich+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="Fungal ectomycorrhizae richness", x="Compartment", y = "Richness")+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=09), legend.position = "bottom") + scale_color_manual(name="compartment", values=c("turquoise","darkviolet"))

### indicator sp ----
# import


### composition ----
#import SLectometa

# remove soil-litter
SLectometa2 <- ectoSLmeta2[ectoSLmeta2$compartment!="soil-litter",]
#write.csv(sapSLmeta3, "Soil-Litter-saprotrophs-metadata.csv")

ecto1 <- ggplot(SLectometa2, aes(x = ectopc1, y = ectopc2, col = habitat)) + geom_point(alpha = 0.4) + stat_ellipse(geom="polygon", aes(fill = habitat), alpha = 0.05, show.legend = T, level = 0.95) + theme_classic() + theme(plot.title = element_text(size=13), legend.position = "top") + facet_wrap(~compartment, labeller = labeller(season = c("litter" = "Litter", "soil" = "Soil"))) + labs(x="PC1", y="PC2")

## pyrofiles ----
### univariate ----
# import SLpyrometa_noSL

meanpyrodiv <- summarySE(SLpyrometa_noSL, measurevar="pyrodiv", groupvars="habitat")
ggplot(meanpyrodiv, aes(x=habitat, y=pyrodiv, color=habitat)) +  geom_point(size =4) + geom_errorbar(aes(ymin=pyrodiv-se,ymax=pyrodiv+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="Fungal pyrophile diversity", x="Veg type", y = "Shannon-Index Diversity")+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=09), legend.position = "none") + facet_grid(scales = "free_x", space = "free_x") + scale_color_manual(name="habitat", values=c("darkgreen","brown","orange")) # in general, pyrophile diversity higher in soil, greatest difference in conifer then grassland then oaks

meanpyrofire <- summarySE(SLpyrometa_noSL, measurevar="pyrodiv", groupvars="fire")
ggplot(meanpyrofire, aes(x=fire, y=pyrodiv, color=fire)) +  geom_point(size =4) + geom_errorbar(aes(ymin=pyrodiv-se,ymax=pyrodiv+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="Fungal pyrophile diversity", x="Fire", y = "Shannon-Index Diversity")+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=09), legend.position = "bottom") + facet_grid(scales = "free_x", space = "free_x") + scale_color_manual(name="fire", values=c("darkred","gold"))

meanpyrocomp <- summarySE(SLpyrometa_noSL, measurevar="pyrodiv", groupvars="compartment")
ggplot(meanpyrocomp, aes(x=compartment, y=pyrodiv, color=compartment)) +  geom_point(size =4) + geom_errorbar(aes(ymin=pyrodiv-se,ymax=pyrodiv+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="Fungal pyrophile diversity", x="Compartment", y = "Shannon-Index Diversity")+ theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=09), legend.position = "bottom") + facet_grid(scales = "free_x", space = "free_x") + scale_color_manual(name="compartment", values=c("turquoise","darkviolet"))


### indicator sp ----
# import

### composition ----
# import SLpyrometa

# remove soil-litter
SLpyrometa2 <- SLpyrometa[SLpyrometa$compartment!="soil-litter",]
#write.csv(sapSLmeta3, "Soil-Litter-saprotrophs-metadata.csv")

pyro1 <- ggplot(SLpyrometa2, aes(x = pyropc1, y = pyropc2, col = habitat)) + geom_point(alpha = 0.4) + stat_ellipse(geom="polygon", aes(fill = habitat), alpha = 0.05, show.legend = T, level = 0.95) + theme_classic() + theme(plot.title = element_text(size=13), legend.position = "top") + facet_wrap(~compartment, labeller = labeller(season = c("litter" = "Litter", "soil" = "Soil"))) + labs(x="PC1", y="PC2")
pyro1

pyro2 <- ggplot(SLpyrometa2, aes(x = pyropc1, y = pyropc2, col = habitat)) + geom_point(alpha = 0.4) + stat_ellipse(geom="polygon", aes(fill = habitat), alpha = 0.05, show.legend = T, level = 0.95) + theme_classic() + theme(plot.title = element_text(size=13), legend.position = "top") + facet_wrap(~season) + labs(x="PC1", y="PC2")
pyro2

ggplot(SLpyrometa2, aes(x = pyropc1, y = pyropc2, col = eco.region)) + geom_point(alpha = 0.3) + stat_ellipse(geom="polygon", alpha = 0, show.legend = T, level = 0.95) + theme_minimal() + labs(title="Saprotroph fungal composition",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "top") + facet_wrap(~season) + theme_classic()

nofirepyro <- SLpyrometa2[(SLpyrometa2$eco.region!="SROC" & SLpyrometa2$eco.region!="NE"),]

ggplot(nofirepyro, aes(x = pyropc1, y = pyropc2, col = fire)) + geom_point(alpha = 0.4) + stat_ellipse(geom="polygon", aes(fill = fire), alpha = 0.05, show.legend = T, level = 0.95) + theme_minimal() + labs(title="Saprotroph composition",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "top") + facet_wrap(~habitat) + theme_classic() + scale_color_manual(name="fire", values=c("maroon","darkturquoise"))


ggplot(nofirepyro, aes(x = pyropc1, y = pyropc2, col = fire)) + geom_point(alpha = 0.4) + stat_ellipse(geom="polygon", aes(fill = fire), alpha = 0.05, show.legend = T, level = 0.95) + theme_minimal() + labs(title="Saprotroph composition",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "top") + facet_wrap(~eco.region) + theme_classic() + scale_color_manual(name="fire", values=c("maroon","darkturquoise"))

ggplot(nofirepyro, aes(x = pyropc1, y = pyropc2, col = fire)) + geom_point(alpha = 0.4) + stat_ellipse(geom="polygon", aes(fill = fire), alpha = 0.05, show.legend = T, level = 0.95) + theme_minimal() + labs(title="Saprotroph composition",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "top") + theme_classic() + scale_color_manual(name="fire", values=c("maroon","darkturquoise"))