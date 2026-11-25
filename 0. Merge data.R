# Merging datasets ----
library(dplyr)
library(base)
library(Matrix)
library(reshape2)
library(tidyr)

## OTU table ----
setwd("raw/")
tocsv <- read.delim("OTU_Table_Soil-Litter_taxonomy.txt", header=T)
head(tocsv)

FT <- read.csv("S1_2020_466_MOESM4_ESM.csv")
FT_search <- FT[,-c(1:7)]

# filter out non fungi
tocsv$taxonomy
tocsv_fung <- tocsv[tocsv$taxonomy!="NonFungal_Unassigned",]
#write.csv(tocsv_fung, "OTU_Table_Soil-Litter_taxonomy.csv")
tocsv_fung$taxonomy

# filter out singletons
tocsv_fungno1 <- tocsv_fung %>%
  filter(rowSums(tocsv_fung[,-1509]) >= 5)
#write.csv(tocsv_fungno1, "OTU_Table_Soil-Litter_taxonomyfilt5.csv")

# Split by comma, then separate prefix and values
tocsv_fung2 <- tocsv_fung %>% separate(taxonomy, c("domain", "phylum", "class", "order", "family", "genus","species"), ",[a-z]:")

SLfung_genus <- tocsv_fung2[,-c(1509:1513,1515)]
SLfung_genus_clean <- SLfung_genus %>%
  extract(genus, into = c("genus", "confidence"), regex = "(.+)\\((.+)\\)")

SLfung_sp <- tocsv_fung2[,-c(1509:1513)]
SLfung_sp_clean <- SLfung_sp %>%
  extract(species, into = c("species", "confidence"), regex = "(.+)\\((.+)\\)") %>%
  extract(genus, into = c("genus", "confidence"), regex = "(.+)\\((.+)\\)")
SLfung_gensp <- SLfung_sp_clean %>%
  mutate(gensp = paste(genus, species, sep=" "))
write.csv(SLfung_gensp, "output/CliMush_SLfung_gensp.csv")

SLfung_genus_clean$GENUS <- SLfung_genus_clean$genus
# fungal traits
matchedFT <- left_join(SLfung_genus_clean, FT_search, by="GENUS")
rownames(matchedFT) <- paste0("OTU", seq_len(nrow(matchedFT)))

# export
SL_asvs <- matchedFT[,-c(1509:1530)]
write.csv(SL_asvs, "output/CliMush_SL_otus.csv")

# export
write.csv(matchedFT, "output/CliMush_SL_otuFT.csv")

## subplot data ----
metadata <- read.csv("CliMush-SoilLitter-metadata-25Feb25.csv", head=T, sep=",")
head(metadata)
metadata$year <- as.integer(metadata$year)
metadata$year <- metadata$year - 2000
metadata <- metadata %>% mutate(month = sub(".*-(.*)$", "\\1", date))
metadata$month <- as.integer(metadata$month)

## climate data ----
# first run Data/WeatherData/CLIMUSH_climate_data_03242025.R
# export climate data, import here
subplot_temp <- read.csv("subplot_tempdata_6Aug25_HB.csv", head=T, sep=",")
head(subplot_temp)
subplot_temp <- subplot_temp[,-1]
subplot_temp$subplot <- as.character(subplot_temp$subplot)

# need to summarize T by day for each month
subplot_temp_month <- subplot_temp %>% 
  group_by(site, month, year, treatment, subplot) %>%
  summarize(iBut_MMT = mean(subplot_T))


metadata_sptemp <- left_join(metadata, subplot_temp_month, by=c("site", "year", "month","subplot"), relationship = "many-to-many") # what's missing?


# export
write.csv(metadata_sptemp, "output/metadata_wSPtemp.csv")

climate.dat <- read.csv("climate_dat_6Aug25_HB.csv", head=T, sep=",")
head(climate.dat)
climate.dat <- climate.dat[,-1]

# summarize by month
climate.dat_month <- climate.dat %>% 
  group_by(site, month, year) %>%
  summarize(across(c(RH_mean, precip_7day, Dew_mean), mean, na.rm = TRUE))
  


# merge climate data with metadata
metadata_tot <- left_join(metadata_sptemp, climate.dat_month, by=c("site", "year", "month"))
# since climate data is on the site level (all subplots at each site have same values), subset by just that subplot (5)
metadata_tot5 <- metadata_tot[metadata_tot$subplot=="5",]
metadata_tot5_otus <- metadata_tot5[!is.na(metadata_tot5$OTU_ID),]
metadata_tot5_otus <- metadata_tot5_otus[!is.na(metadata_tot5_otus$MAT),]

# export
write.csv(metadata_tot, "output/CliMush_climatedat.csv")

# add predictors ----
# for C:N & pH, subset into 2022 spring and merge by subplot
pred_dat <- read.csv("pred_dat_16June25.csv", header=T)
head(pred_dat)
pred_dat$subplot <- as.character(pred_dat$subplot)

otu <- dplyr::select(metadata, site, habitat, fire, subplot, compartment, year, season, OTU_ID)
pred_dat_wotu <- left_join(pred_dat, otu, by=c("site", "habitat", "fire", "subplot"))

soil_pred_dat_wotu <- pred_dat_wotu[pred_dat_wotu$compartment=="soil",]
# soil
soilpred <- data.frame(
  OTU_ID <- (soil_pred_dat_wotu$OTU_ID),
  year <- (soil_pred_dat_wotu$year),
  season <- (soil_pred_dat_wotu$season),
  site <- (soil_pred_dat_wotu$eco.region),
  subplot <- (soil_pred_dat_wotu$subplot),
  habitat <- (soil_pred_dat_wotu$habitat),
  fire <- (soil_pred_dat_wotu$fire),
  PercMoist <- (soil_pred_dat_wotu$spring_soil_moist),
  CN <- (soil_pred_dat_wotu$soil_CN),
  pH <- (soil_pred_dat_wotu$soil_pH),
  elevation <- (soil_pred_dat_wotu$Elevation.m),
  tree_sp <- (soil_pred_dat_wotu$tree_sp),
  tree_fam <- (soil_pred_dat_wotu$tree_fam),
  nonveg <- (soil_pred_dat_wotu$nonveg_perc),
  plnt1_sp <- (soil_pred_dat_wotu$plntsp1),
  plnt1_fam <- (soil_pred_dat_wotu$plntsp1_fam),
  plnt1_perc <- (soil_pred_dat_wotu$plntsp1_perc),
  plnt2_sp <- (soil_pred_dat_wotu$plntsp2),
  plnt2_fam <- (soil_pred_dat_wotu$plntsp2_fam),
  plnt2_perc <- (soil_pred_dat_wotu$plntsp2_perc),
  plnt3_sp <- (soil_pred_dat_wotu$plntsp3),
  plnt3_fam <- (soil_pred_dat_wotu$plntsp3_fam),
  plnt3_perc <- (soil_pred_dat_wotu$plntsp3_perc))
#View(soilpred)
colnames(soilpred) <- c("OTU_ID", "year", "season", "site", "subplot", "habitat", "fire", "PercMoist", "CN", "pH", "elevation", "tree_sp", "tree_fam", "nonveg", "plnt1_sp", "plnt1_fam", "plntsp1_perc", "plnt2_sp", "plnt2_fam", "plntsp2_perc", "plnt3_sp", "plnt3_fam", "plntsp3_perc")
write.csv(soilpred, "output/Soil_dat.csv")

litter_pred_dat_wotu <- pred_dat_wotu[pred_dat_wotu$compartment=="litter",]
# litter
litterpred <- data.frame(
  OTU_ID <- (litter_pred_dat_wotu$OTU_ID),
  year <- (litter_pred_dat_wotu$year),
  season <- (litter_pred_dat_wotu$season),
  site <- (litter_pred_dat_wotu$eco.region),
  subplot <- (litter_pred_dat_wotu$subplot),
  habitat <- (litter_pred_dat_wotu$habitat),
  fire <- (litter_pred_dat_wotu$fire),
  PercMoist <- (litter_pred_dat_wotu$spring_litter_moist),
  CN <- (litter_pred_dat_wotu$litter_CN),
  pH <- (litter_pred_dat_wotu$litter_pH),
  elevation <- (litter_pred_dat_wotu$Elevation.m),
  tree_sp <- (litter_pred_dat_wotu$tree_sp),
  tree_fam <- (litter_pred_dat_wotu$tree_fam),
  nonveg <- (litter_pred_dat_wotu$nonveg_perc),
  plnt1_sp <- (litter_pred_dat_wotu$plntsp1),
  plnt1_fam <- (litter_pred_dat_wotu$plntsp1_fam),
  plnt1_perc <- (litter_pred_dat_wotu$plntsp1_perc),
  plnt2_sp <- (litter_pred_dat_wotu$plntsp2),
  plnt2_fam <- (litter_pred_dat_wotu$plntsp2_fam),
  plnt2_perc <- (litter_pred_dat_wotu$plntsp2_perc),
  plnt3_sp <- (litter_pred_dat_wotu$plntsp3),
  plnt3_fam <- (litter_pred_dat_wotu$plntsp3_fam),
  plnt3_perc <- (litter_pred_dat_wotu$plntsp3_perc))
#View(litterpred)
colnames(litterpred) <- c("OTU_ID","year", "season", "site", "subplot", "habitat", "fire", "PercMoist", "CN", "pH", "elevation", "tree_sp", "tree_fam", "nonveg", "plnt1_sp", "plnt1_fam", "plntsp1_perc", "plnt2_sp", "plnt2_fam", "plntsp2_perc", "plnt3_sp", "plnt3_fam", "plntsp3_perc")

write.csv(litterpred, "output/Litter_dat.csv")