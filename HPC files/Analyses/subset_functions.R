# Subset functional groups ----
library(vegan)
library(dplyr)

pyrolist <- read.csv("Pyrophilous_list-SL19Mar25.csv", head=T)

# read in SL_FT & metadata
SL_FT <- read.csv("CliMush_SL_otuFT.csv", head=T, sep=",")
head(SL_FT)
SL_FT <- SL_FT[,-1]

metadata <- read.csv("CliMush-SoilLitter-metadata-25Feb25.csv", head=T, sep=",")

SL_asv <- SL_FT[,-c(1509:1544)]

# Select numeric sample columns
d <- colnames(SL_asv[, sapply(SL_asv, is.numeric)])
n <- length(d)

# Prepare output data frame
columns <- c("Sample", "fungrel_ab", "fungdiv", "fungrich")
outputfung <- data.frame(matrix(nrow = n, ncol = length(columns)))
colnames(outputfung) <- columns

# Loop through each sample
for (i in 1:n) {
  sample_name <- d[i]
  total_fungi <- sum(SL_asv[, sample_name], na.rm = TRUE)
  total_reads <- sum(SL_asv, na.rm = TRUE)
  rel_ab <- total_fungi / total_reads
  div_val <- diversity(SL_asv[, sample_name], index = "shannon")
  rich_val <- specnumber(SL_asv[, sample_name])
  
  # Store results
  outputfung$Sample[i]    <- sample_name
  outputfung$fungrel_ab[i] <- rel_ab
  outputfung$fungdiv[i]    <- div_val
  outputfung$fungrich[i]   <- rich_val
}


metadata1 <- left_join(metadata, outputfung, by=c("OTU_ID" = "Sample"))


## pathogens ----
SL_FT$ppath1 <- ifelse(SL_FT$primary_lifestyle=="plant_pathogen", "1", "0")

SL_FT$ppath2 <- ifelse(SL_FT$primary_lifestyle=="plant_pathogen" | SL_FT$Secondary_lifestyle=="plant_pathogen", "1", "0")

## calculate relative abundance of pathogens/sample
d <- colnames(SL_FT[,sapply(SL_FT, is.numeric)])
n<-length(d)
columns = c("pprel_ab1","ppdiv1","pprich1","pprel_ab2","ppdiv2","pprich2") 
outputpath = data.frame(matrix(nrow = n, ncol = length(columns)))
colnames(outputpath) = columns

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$ppath1=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputpath$pprel_ab1[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$ppath1=="1"), d[i]])
  df <- as.data.frame(s)
  outputpath$ppdiv1[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$ppath1=="1"), d[i]])
  df <- as.data.frame(s)
  outputpath$pprich1[i] <- df[[1]]}

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$ppath2=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputpath$pprel_ab2[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$ppath2=="1"), d[i]])
  df <- as.data.frame(s)
  outputpath$ppdiv2[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$ppath2=="1"), d[i]])
  df <- as.data.frame(s)
  outputpath$pprich2[i] <- df[[1]]}

outputpath$Sample <- d # now you have a dataframe with sample ID as one column and relative abundance of rhizobia reads as the other column

metadata2 <- left_join(metadata1, outputpath, by=c("OTU_ID" = "Sample"))

## sap tot ----
SL_FT$allsaps1 <- ifelse(SL_FT$primary_lifestyle=="dung_saprotroph" | SL_FT$primary_lifestyle=="litter_saprotroph" | SL_FT$primary_lifestyle=="soil_saprotroph" | SL_FT$primary_lifestyle=="wood_saprotroph", "1", "0")

SL_FT$allsaps2 <- ifelse(SL_FT$primary_lifestyle=="dung_saprotroph" | SL_FT$primary_lifestyle=="litter_saprotroph" | SL_FT$primary_lifestyle=="soil_saprotroph" | SL_FT$primary_lifestyle=="wood_saprotroph" | SL_FT$Secondary_lifestyle=="dung_saprotroph" | SL_FT$Secondary_lifestyle=="litter_saprotroph" | SL_FT$Secondary_lifestyle=="soil_saprotroph" | SL_FT$Secondary_lifestyle=="wood_saprotroph", "1", "0")


## calculate relative abundance of pathogens/sample
d <- colnames(SL_FT[,sapply(SL_FT, is.numeric)])
p<-length(d)
columns = c("allsaps_ab1","allsapsdiv1","allsapsrich1","allsaps_ab2","allsapsdiv2","allsapsrich2") 
outputallsaps = data.frame(matrix(nrow = p, ncol = length(columns)))
colnames(outputallsaps) = columns

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$allsaps1=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputallsaps$allsaps_ab1[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$allsaps1=="1"), d[i]])
  df <- as.data.frame(s)
  outputallsaps$allsapsdiv1[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$allsaps1=="1"), d[i]])
  df <- as.data.frame(s)
  outputallsaps$allsapsrich1[i] <- df[[1]]}

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$allsaps2=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputallsaps$allsaps_ab2[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$allsaps2=="1"), d[i]])
  df <- as.data.frame(s)
  outputallsaps$allsapsdiv2[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$allsaps2=="1"), d[i]])
  df <- as.data.frame(s)
  outputallsaps$allsapsrich2[i] <- df[[1]]}

outputallsaps$Sample <- d # now you have a dataframe with sample ID as one column and relative abundance of rhizobia reads as the other column

metadata3 <- left_join(metadata2, outputallsaps, by=c("OTU_ID" = "Sample"))

## dung saps ----
SL_FT$dungsap1 <- ifelse(SL_FT$primary_lifestyle=="dung_saprotroph", "1", "0")
SL_FT$dungsap2 <- ifelse(SL_FT$primary_lifestyle=="dung_saprotroph"| SL_FT$Secondary_lifestyle=="dung_saprotroph", "1", "0")

## calculate relative abundance of pathogens/sample
d <- colnames(SL_FT[,sapply(SL_FT, is.numeric)])
p<-length(d)
columns = c("dungsaprel_ab1","dungsapdiv1","dungsaprich1","dungsaprel_ab2","dungsapdiv2","dungsaprich2") 
outputdungsap = data.frame(matrix(nrow = p, ncol = length(columns)))
colnames(outputdungsap) = columns

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$dungsap1=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputdungsap$dungsaprel_ab1[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$dungsap1=="1"), d[i]])
  df <- as.data.frame(s)
  outputdungsap$dungsapdiv1[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$dungsap1=="1"), d[i]])
  df <- as.data.frame(s)
  outputdungsap$dungsaprich1[i] <- df[[1]]}

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$dungsap2=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputdungsap$dungsaprel_ab2[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$dungsap2=="1"), d[i]])
  df <- as.data.frame(s)
  outputdungsap$dungsapdiv2[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$dungsap2=="1"), d[i]])
  df <- as.data.frame(s)
  outputdungsap$dungsaprich2[i] <- df[[1]]}

outputdungsap$Sample <- d # now you have a dataframe with sample ID as one column and relative abundance of rhizobia reads as the other column

metadata4 <- left_join(metadata3, outputdungsap, by=c("OTU_ID" = "Sample"))

## litter saps ----
SL_FT$litsap1 <- ifelse(SL_FT$primary_lifestyle=="litter_saprotroph", "1", "0")
SL_FT$litsap2 <- ifelse(SL_FT$primary_lifestyle=="litter_saprotroph"| SL_FT$Secondary_lifestyle=="litter_saprotroph", "1", "0")

## calculate relative abundance of pathogens/sample
d <- colnames(SL_FT[,sapply(SL_FT, is.numeric)])
p<-length(d)
columns = c("litsaprel_ab1","litsapdiv1","litsaprich1","litsaprel_ab2","litsapdiv2","litsaprich2") 
outputlitsap = data.frame(matrix(nrow = p, ncol = length(columns)))
colnames(outputlitsap) = columns

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$litsap1=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputlitsap$litsaprel_ab1[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$litsap1=="1"), d[i]])
  df <- as.data.frame(s)
  outputlitsap$litsapdiv1[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$litsap1=="1"), d[i]])
  df <- as.data.frame(s)
  outputlitsap$litsaprich1[i] <- df[[1]]}

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$litsap2=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputlitsap$litsaprel_ab2[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$litsap2=="1"), d[i]])
  df <- as.data.frame(s)
  outputlitsap$litsapdiv2[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$litsap2=="1"), d[i]])
  df <- as.data.frame(s)
  outputlitsap$litsaprich2[i] <- df[[1]]}

outputlitsap$Sample <- d # now you have a dataframe with sample ID as one column and relative abundance of rhizobia reads as the other column

metadata5 <- left_join(metadata4, outputlitsap, by=c("OTU_ID" = "Sample"))

## soil saps ----
SL_FT$soilsap1 <- ifelse(SL_FT$primary_lifestyle=="soil_saprotroph", "1", "0")
SL_FT$soilsap2 <- ifelse(SL_FT$primary_lifestyle=="soil_saprotroph" | SL_FT$Secondary_lifestyle=="soil_saprotroph", "1", "0")

## calculate relative abundance of pathogens/sample
d <- colnames(SL_FT[,sapply(SL_FT, is.numeric)])
p<-length(d)
columns = c("soilsaprel_ab1","soilsapdiv1","soilsaprich1","soilsaprel_ab2","soilsapdiv2","soilsaprich2") 
outputsoilsap = data.frame(matrix(nrow = p, ncol = length(columns)))
colnames(outputsoilsap) = columns

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$soilsap1=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputsoilsap$soilsaprel_ab1[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$soilsap1=="1"), d[i]])
  df <- as.data.frame(s)
  outputsoilsap$soilsapdiv1[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$soilsap1=="1"), d[i]])
  df <- as.data.frame(s)
  outputsoilsap$soilsaprich1[i] <- df[[1]]}

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$soilsap2=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputsoilsap$soilsaprel_ab2[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$soilsap2=="1"), d[i]])
  df <- as.data.frame(s)
  outputsoilsap$soilsapdiv2[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$soilsap2=="1"), d[i]])
  df <- as.data.frame(s)
  outputsoilsap$soilsaprich2[i] <- df[[1]]}

outputsoilsap$Sample <- d # now you have a dataframe with sample ID as one column and relative abundance of rhizobia reads as the other column

metadata6 <- left_join(metadata5, outputsoilsap, by=c("OTU_ID" = "Sample"))

## wood saps ----
SL_FT$woodsap1 <- ifelse(SL_FT$primary_lifestyle=="wood_saprotroph", "1", "0")
SL_FT$woodsap2 <- ifelse(SL_FT$primary_lifestyle=="wood_saprotroph" | SL_FT$Secondary_lifestyle=="wood_saprotroph", "1", "0")

## calculate relative abundance of pathogens/sample
d <- colnames(SL_FT[,sapply(SL_FT, is.numeric)])
p<-length(d)
columns = c("woodsaprel_ab1","woodsapdiv1","woodsaprich1","woodsaprel_ab2","woodsapdiv2","woodsaprich2") 
outputwoodsap = data.frame(matrix(nrow = p, ncol = length(columns)))
colnames(outputwoodsap) = columns

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$woodsap1=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputwoodsap$woodsaprel_ab1[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$woodsap1=="1"), d[i]])
  df <- as.data.frame(s)
  outputwoodsap$woodsapdiv1[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$woodsap1=="1"), d[i]])
  df <- as.data.frame(s)
  outputwoodsap$woodsaprich1[i] <- df[[1]]}

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$woodsap2=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputwoodsap$woodsaprel_ab2[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$woodsap2=="1"), d[i]])
  df <- as.data.frame(s)
  outputwoodsap$woodsapdiv2[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$woodsap2=="1"), d[i]])
  df <- as.data.frame(s)
  outputwoodsap$woodsaprich2[i] <- df[[1]]}

outputwoodsap$Sample <- d # now you have a dataframe with sample ID as one column and relative abundance of rhizobia reads as the other column

metadata7 <- left_join(metadata6, outputwoodsap, by=c("OTU_ID" = "Sample"))

## ectomycorrhizae ----
SL_FT$ecto <- ifelse(SL_FT$primary_lifestyle=="ectomycorrhizal", "1", "0")

d <- colnames(SL_FT[,sapply(SL_FT, is.numeric)])
n<-length(d)
columns = c("ectorel_ab","ectodiv","ectorich") 
outputecto = data.frame(matrix(nrow = n, ncol = length(columns)))
colnames(outputecto) = columns

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$ecto=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputecto$ectorel_ab[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$ecto=="1"), d[i]])
  df <- as.data.frame(s)
  outputecto$ectodiv[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$ecto=="1"), d[i]])
  df <- as.data.frame(s)
  outputecto$ectorich[i] <- df[[1]]}

outputecto$Sample <- d # now you have a dataframe with sample ID as one column and relative abundance of rhizobia reads as the other column

metadata8 <- left_join(metadata7, outputecto, by=c("OTU_ID" = "Sample"))

## pyrofiles ----
pyrosp <- pyrolist$Pyrophilous_Species

# match with species
SL_full <-read.csv("CliMush_SLfung_gensp.csv", header=T)
head(SL_full)
SL_full$pyro_present <- SL_full$gensp %in% pyrosp
SL_full %>% count(pyro_present)

SL_FT$pyro <- ifelse(SL_full$pyro_present=="TRUE", "1", "0")

## calculate relative abundance of pathogens/sample
d <- colnames(SL_FT[,sapply(SL_FT, is.numeric)])
n<-length(d)
columns = c("pyrorel_ab","pyrodiv","pyrorich") 
outputpyro = data.frame(matrix(nrow = n, ncol = length(columns)))
colnames(outputpyro) = columns

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$pyro=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputpyro$pyrorel_ab[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$pyro=="1"), d[i]])
  df <- as.data.frame(s)
  outputpyro$pyrodiv[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$pyro=="1"), d[i]])
  df <- as.data.frame(s)
  outputpyro$pyrorich[i] <- df[[1]]}

outputpyro$Sample <- d # now you have a dataframe with sample ID as one column and relative abundance of rhizobia reads as the other column

metadata9 <- left_join(metadata8, outputpyro, by=c("OTU_ID" = "Sample"))

# export data ----
write.csv(metadata9, "CliMush-SoilLitter-metaWfunctions_14Aug25.csv")
write.csv(SL_FT, "CliMush_SL_otu_functions1.csv")
