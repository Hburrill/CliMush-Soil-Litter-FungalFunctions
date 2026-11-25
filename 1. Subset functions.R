# Subset functional groups ----
library(vegan)
library(dplyr)
library(data.table) # fread/fwrite: much faster and lighter on memory than read.csv/write.csv

# Run from the Analyses folder (open Analyses.Rproj, or the HPC script cd's there).
# Paths are relative to it, so there's no setwd() and the script can be re-run in the same session.
pyrolist <- read.csv("raw/Pyrophilous_list-SL19Mar25.csv", head=T)

# read in SL_FT & metadata
SL_FT <- fread("output/CliMush_SL_otuFT.csv", data.table = FALSE, check.names = TRUE)
head(SL_FT[, c(1:3, (ncol(SL_FT)-21):ncol(SL_FT))])
SL_FT <- SL_FT[,-1]

# Sample columns, picked by name so this keeps working if samples are added/removed
sample_cols <- grep("^illumina_", names(SL_FT), value = TRUE)
length(sample_cols) # 1508 as of Nov 2026

# how many total ASVs?
# 821317

# how many pathogens(primary)?
sum(SL_FT$primary_lifestyle == "plant_pathogen", na.rm = TRUE) # 150398

# saps?
sum(SL_FT$primary_lifestyle=="dung_saprotroph" | SL_FT$primary_lifestyle=="litter_saprotroph" | SL_FT$primary_lifestyle=="soil_saprotroph" | SL_FT$primary_lifestyle=="wood_saprotroph", na.rm = TRUE) # 266322

# ecto?
sum(SL_FT$primary_lifestyle == "ectomycorrhizal", na.rm = TRUE) # 28521


## subset for tests
# Commented out: mutate_all() makes a full numeric copy of the table (~10 GB), which
# doesn't fit in a laptop's memory, and otu.table.filter isn't used later in this script.
# SL_FT2 <- SL_FT[which(SL_FT$confidence>=0.9), sample_cols]
# otu.table.filter <- SL_FT2[,colSums(SL_FT2)>10]

metadata <- read.csv("output/CliMush-SoilLitter-metadata-25Feb25.csv", head=T, sep=",")

## allfung ----
d <- sample_cols
n <- length(d)

# Total reads across all samples: computed once, one column at a time (sum() on the whole
# data frame inside the loop copied the full table on every pass)
total_reads <- sum(vapply(SL_FT[d], function(x) as.numeric(sum(x, na.rm = TRUE)), numeric(1)))

# Prepare output data frame
columns <- c("Sample", "fungrel_ab", "fungdiv", "fungrich")
outputfung <- data.frame(matrix(nrow = n, ncol = length(columns)))
colnames(outputfung) <- columns

# Loop through each sample
for (i in 1:n) {
  sample_name <- d[i]
  total_fungi <- sum(SL_FT[, sample_name], na.rm = TRUE)
  rel_ab <- total_fungi / total_reads
  div_val <- vegan::diversity(SL_FT[, sample_name])
  rich_val <- vegan::specnumber(SL_FT[, sample_name])

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
d <- sample_cols
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
d <- sample_cols
n<-length(d)
columns = c("allsaps_ab1","allsapsdiv1","allsapsrich1","allsaps_ab2","allsapsdiv2","allsapsrich2") 
outputallsaps = data.frame(matrix(nrow = n, ncol = length(columns)))
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
d <- sample_cols
n<-length(d)
columns = c("dungsaprel_ab1","dungsapdiv1","dungsaprich1","dungsaprel_ab2","dungsapdiv2","dungsaprich2") 
outputdungsap = data.frame(matrix(nrow = n, ncol = length(columns)))
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
d <- sample_cols
n<-length(d)
columns = c("litsaprel_ab1","litsapdiv1","litsaprich1","litsaprel_ab2","litsapdiv2","litsaprich2") 
outputlitsap = data.frame(matrix(nrow = n, ncol = length(columns)))
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
d <- sample_cols
n<-length(d)
columns = c("soilsaprel_ab1","soilsapdiv1","soilsaprich1","soilsaprel_ab2","soilsapdiv2","soilsaprich2") 
outputsoilsap = data.frame(matrix(nrow = n, ncol = length(columns)))
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
d <- sample_cols
n<-length(d)
columns = c("woodsaprel_ab1","woodsapdiv1","woodsaprich1","woodsaprel_ab2","woodsapdiv2","woodsaprich2") 
outputwoodsap = data.frame(matrix(nrow = n, ncol = length(columns)))
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

d <- sample_cols
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
# Only the genus/species columns are needed, so read just those (the full file is ~2.5 GB)
SL_full <- fread("output/CliMush_SLfung_gensp.csv", select = c("genus", "gensp"), data.table = FALSE)
head(SL_full)
SL_full$pyro_present <- SL_full$gensp %in% pyrosp
SL_full %>% count(pyro_present) # 5889

# This file has no shared OTU ID with SL_FT, so rows are matched by position.
# Stop if the two files are out of sync rather than silently mislabeling taxa.
stopifnot(nrow(SL_full) == nrow(SL_FT),
          identical(as.character(SL_full$genus), as.character(SL_FT$genus)))

SL_FT$pyro <- ifelse(SL_full$pyro_present=="TRUE", "1", "0")
# pyro?


## calculate relative abundance of pathogens/sample
d <- sample_cols
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

## pyro_subset ----
# Listed pyrophiles (pyro == 1) in the genera Pyronema, Neurospora and Geopyxis only,
# to test whether overall pyrophile trends hold for these core groups.
pyro_genera <- c("Pyronema", "Neurospora", "Geopyxis")
SL_FT$pyro_subset <- ifelse(SL_FT$pyro=="1" & SL_FT$genus %in% pyro_genera, "1", "0")
sum(SL_FT$pyro_subset=="1") # 24 OTUs (P. domesticum, N. crassa, G. majalis/alpina/carbonaria)

d <- sample_cols
n<-length(d)
columns = c("pyrosubrel_ab","pyrosubdiv","pyrosubrich")
outputpyrosub = data.frame(matrix(nrow = n, ncol = length(columns)))
colnames(outputpyrosub) = columns

for(i in 1:n){
  s <- sum(SL_FT[which(SL_FT$pyro_subset=="1"), d[i]])/sum(SL_FT[,d[i]])
  df <- as.data.frame(s)
  outputpyrosub$pyrosubrel_ab[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(SL_FT[which(SL_FT$pyro_subset=="1"), d[i]])
  df <- as.data.frame(s)
  outputpyrosub$pyrosubdiv[i] <- df[[1]]}
for(i in 1:n){
  s <- specnumber(SL_FT[which(SL_FT$pyro_subset=="1"), d[i]])
  df <- as.data.frame(s)
  outputpyrosub$pyrosubrich[i] <- df[[1]]}

outputpyrosub$Sample <- d

metadata10 <- left_join(metadata9, outputpyrosub, by=c("OTU_ID" = "Sample"))

# export data ----
write.csv(metadata10, "output/CliMush-SoilLitter-metaWfunctions_24Nov26.csv")
# fwrite is far faster for this ~2.6 GB table; row.names=TRUE keeps the same layout as write.csv
fwrite(SL_FT, "output/CliMush_SL_otu_functions_24Nov26.csv", row.names = TRUE)
