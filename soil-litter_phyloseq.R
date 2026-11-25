library(phyloseq)
library(dplyr)
library(tidyr)

setwd("C:/Users/haley/University of Oregon Dropbox/Haley Burrill/CLIMUSH NSF-Neon/Analysis/SoilLitterFunctionPaper_HB/Analyses/output")
# Extract OTU table into a df
ps <- read.csv("CliMush_SL_otuFT.csv", header=T)
otu_tab <- ps[,-c(1510:1531)]
rownames(otu_tab) <- otu_tab$X.1
otu_tab$X.1 <- NULL
otu_df <- as.data.frame(otu_table(otu_tab, taxa_are_rows = T))

# Extract taxonomy
taxa <- read.csv("CliMush_SLfung_gensp.csv", header=T)
taxa_gensp <- taxa[,c(1510,1512)]
tax_df <- as.data.frame(tax_table(taxa_gensp))
rownames(tax_df) <- rownames(otu_tab)

# Merge OTU counts + taxonomy
df <- otu_df
df$Genus <- tax_df$ta1
df$Species <- tax_df$ta2
df$primary_function <- ps$primary_lifestyle

pyrolist <- read.csv("C:/Users/haley/University of Oregon Dropbox/Haley Burrill/CLIMUSH NSF-Neon/Analysis/SoilLitterFunctionPaper_HB/Analyses/raw/Pyrophilous_list-SL19Mar25.csv", head=T)
pyrosp <- pyrolist$Pyrophilous_Species

# match with species
taxa$pyro_present <- taxa$gensp %in% pyrosp
taxa %>% count(pyro_present) # 5889

df$pyro <- taxa$pyro_present

## overlap ----
otu_pa <- otu_table(otu_df)
if (!taxa_are_rows(ps_pa)) otu_pa <- t(otu_pa)

meta <- sample_data(ps_pa) %>% as.data.frame()

soil_samples   <- rownames(meta)[meta$Compartment == "Soil"]
litter_samples <- rownames(meta)[meta$Compartment == "Litter"]

soil_otus <- rownames(otu_pa)[rowSums(otu_pa[, soil_samples]) > 0]
litter_otus <- rownames(otu_pa)[rowSums(otu_pa[, litter_samples]) > 0]


## list ----
# Sum abundance across all samples
df2 <- df %>%
  mutate(total_abundance = rowSums(df[,-c(1509:1512)]))

# Filter only the groups you care about
target_groups <- c("plant_pathogen", "litter_saprotroph", "wood_saprotroph", "soil_saprotroph", "ectomycorrhizal")

df_filtered <- df2 %>% 
  filter(primary_function %in% target_groups)

write.csv(df_filtered, "C:/Users/haley/University of Oregon Dropbox/Haley Burrill/CLIMUSH NSF-Neon/Analysis/SoilLitterFunctionPaper_HB/Analyses/output/functionalgroups_only.csv")

df_long <- df_filtered |>
  pivot_longer(
    cols = -all_of(tax_cols),
    names_to = "sample_id",
    values_to = "abundance"
  )


litter_df <- df_filtered[, grepl("^illumina_litter_", colnames(df_filtered))]

soil_df <- df_filtered[, grepl("^illumina_soil_", colnames(df_filtered))]


# Get top 10 in each primary function
top10 <- df_filtered %>%
  group_by(primary_function, pyro) %>%
  slice_max(order_by = total_abundance, n = 10)

write.csv(top10[,-c(1:1508)], "top10taxa_perfunc.csv")

top10lit <- litter_df %>%
  group_by(primary_function, pyro) %>%
  slice_max(order_by = total_abundance, n = 10)

write.csv(top10lit[,-c(1:1508)], "top10littertaxa_perfunc.csv")

top10soil <- soil_df %>%
  group_by(primary_function, pyro) %>%
  slice_max(order_by = total_abundance, n = 10)

write.csv(top10soil[,-c(1:1508)], "top10soiltaxa_perfunc.csv")
