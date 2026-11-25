# Packages ----
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(base)
library(glmulti)
library(lme4)
library(Matrix)
library(CCA)
library(reshape2)
library(indicspecies)
library(usmap)
library(ggpubr)
set.seed(123)

# import data ----
setwd("C:/Users/haley/University of Oregon Dropbox/Haley Burrill/CLIMUSH NSF-Neon/Analysis/SoilLitterFunctionPaper_HB/Analyses/output")
# this has functional group data (richness, diversity, RA)
metadata1 <- read.csv("CliMush-SoilLitter-metaWfunctions_13Aug25.csv", header=T, sep=",")
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
#head(SL_FT)
#SL_FT <- SL_FT[,-1]
## subset for tests
#SL_FT <- SL_FT[-c(20000:821300),]
#write.csv(SL_FT, "SL_FT_subset.csv")
#SL_FT <- read.csv("SL_FT_subset.csv", head=T, sep=",")
#head(SL_FT)
#SL_FT <- SL_FT[,-1]

#SL_FT2 <- mutate_all(SL_FT, function(x) as.numeric(as.character(x)))
#SL_FT3 <- SL_FT2[,-c(1,1509:1543)]
#otu.table.filter <- SL_FT3[,colSums(SL_FT3)>10]


metadata$prop_ppath <- metadata$pprich2/metadata$fungrich
metadata$prop_sap <- metadata$allsapsrich2/metadata$fungrich
metadata$prop_ecto <- metadata$ectorich/metadata$fungrich
metadata$prop_pyro <- metadata$pyrorich/metadata$fungrich
metadata$prop_other <- 1 - metadata$prop_ppath - metadata$prop_sap - metadata$prop_ecto - metadata$prop_pyro

# reshape to long format
propdat <- metadata %>%
  pivot_longer(
    cols = c(prop_ppath, prop_sap, prop_ecto, prop_pyro, prop_other),
    names_to = "funct",
    values_to = "prop"
  ) %>%
  filter(prop >= 0)

propdat$funct <- factor(propdat$funct, levels=c("prop_ppath", "prop_sap", "prop_ecto", "prop_pyro", "prop_other"))
#write.csv(propdat, "prop_data.csv")

propdat <- read.csv("prop_data3.csv", head=T)
head(propdat)
propdat$funct <- factor(propdat$funct, levels=c("prop_sap", "prop_ppath", "prop_ecto", "prop_pyro", "prop_other"))

# between compartments
compartments <- ggplot(propdat, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "bottom", plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values=c("darkgoldenrod4","red","darkgreen","coral","blue")) + facet_wrap(~compartment, scales = "free_y") + labs(title="")
ggsave("compartment_pies.png", path="figures/", width=6, height =4)


# between ecoregions
propdat_DSW <- propdat[propdat$state=="AZ",]
propdat_GRTL <- propdat[propdat$state=="MN",]
propdat_NE <- propdat[propdat$state=="MA",]
propdat_PNW <- propdat[propdat$ecoregion=="Pacific Northwest",]
propdat_PRAI <- propdat[propdat$state=="KS",]
propdat_SE <- propdat[propdat$state=="FL",]
propdat_SROC <- propdat[propdat$state=="CO",]
propdat_TAIG <- propdat[propdat$state=="AK",]

# pie chart ----
# all sites by habitat & fire
DSW_B <- propdat_DSW[propdat_DSW$fire=="B",]
DSW_UB <- propdat_DSW[propdat_DSW$fire=="UB",]

DSWB <- ggplot(DSW_B, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values=c("darkgoldenrod4","red","darkgreen","coral","blue")) + facet_wrap(~habitat, scales = "free_y") + labs(title="DSW - burned") 
DSWUB <- ggplot(DSW_UB, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values=c("darkgoldenrod4","red","darkgreen","coral","blue")) + facet_wrap(~habitat, scales = "free_y") + labs(title="unburned") 

NEUB <- ggplot(propdat_NE, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values=c("darkgoldenrod4","red","darkgreen","coral","blue")) + facet_wrap(~habitat, scales = "free_y")  + labs(title="NE - unburned only") 

PNW_B <- propdat_PNW[propdat_PNW$fire=="B",]
PNW_UB <- propdat_PNW[propdat_PNW$fire=="UB",]

PNWB <- ggplot(PNW_B, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values=c("darkgoldenrod4","red","darkgreen","coral","blue")) + facet_wrap(~habitat, scales = "free_y")  + labs(title="PNW - burned") 
PNWUB <- ggplot(PNW_UB, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values=c("darkgoldenrod4","red","darkgreen","coral","blue")) + facet_wrap(~habitat, scales = "free_y")  + labs(title="unburned") 

PRAI_B <- propdat_PRAI[propdat_PRAI$fire=="B",]
PRAI_UB <- propdat_PRAI[propdat_PRAI$fire=="UB",]

PRAIB <- ggplot(PRAI_B, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values=c("darkgoldenrod4","red","darkgreen","coral","blue")) + facet_wrap(~habitat, scales = "free_y")  + labs(title="PRAI - burned") 
PRAIUB <- ggplot(PRAI_UB, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values=c("darkgoldenrod4","red","darkgreen","coral","blue")) + facet_wrap(~habitat, scales = "free_y")  + labs(title="unburned") 

SE_B <- propdat_SE[propdat_SE$fire=="B",]
SE_UB <- propdat_SE[propdat_SE$fire=="UB",]

SEB <- ggplot(SE_B, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values=c("darkgoldenrod4","red","darkgreen","coral","blue")) + facet_wrap(~habitat, scales = "free_y")  + labs(title="SE - burned") 
SEUB <- ggplot(SE_UB, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values=c("darkgoldenrod4","red","darkgreen","coral","blue")) + facet_wrap(~habitat, scales = "free_y")  + labs(title="unburned") 

SROCUB <- ggplot(propdat_SROC, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values=c("darkgoldenrod4","red","darkgreen","coral","blue")) + facet_wrap(~habitat, scales = "free_y") + labs(title="SROC - unburned only") 

TAIG_B <- propdat_TAIG[propdat_TAIG$fire=="B",]
TAIG_UB <- propdat_TAIG[propdat_TAIG$fire=="UB",]

TAIGB <- ggplot(TAIG_B, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values=c("darkgoldenrod4","red","darkgreen","coral","blue")) + facet_wrap(~habitat, scales = "free_y")  + labs(title="TAIG - burned") 
TAIGUB <- ggplot(TAIG_UB, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values=c("darkgoldenrod4","red","darkgreen","coral","blue")) + facet_wrap(~habitat, scales = "free_y")  + labs(title="unburned") 

p1 <- ggplot() + theme_void()

allfirecomb <- ggarrange(DSWB, DSWUB, p1, NEUB, PNWB, PNWUB, PRAIB, PRAIUB, SEB, SEUB, p1, SROCUB, TAIGB, TAIGUB, ncol = 2, nrow = 8, common.legend = TRUE, legend="top") + bgcolor("white")
ggsave("allpiescomb.png", path="figures/", width=6, height =12)




ggplot(propdat_DSW, aes(x="", y=prop, fill=funct)) +
  geom_bar(stat="identity", mapping = aes(x = habitat)) +
  coord_polar() +
  theme_minimal() +
  theme(legend.position="top") +
  xlab("") + ylab("")

DSW <- ggplot(propdat_DSW, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "mm")) + scale_fill_manual(values=c("darkred","orange","darkgreen","magenta","darkgrey")) + facet_wrap(~habitat, scales = "free_y")

ggsave("DSW_pie.png", path="figures/", width = 6, height = 4, dpi = 300)

GRTL <- ggplot(propdat_GRTL, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "mm")) + scale_fill_manual(values=c("darkred","orange","darkgreen","magenta","darkgrey")) + facet_wrap(~habitat, scales = "free_y")

ggsave("GRTL_pie.png", path="figures/", bg = "white", width = 6, height = 4, dpi = 300)
  
NE <- ggplot(propdat_NE, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "mm")) + scale_fill_manual(values=c("darkred","orange","darkgreen","magenta","darkgrey")) + facet_wrap(~habitat, scales = "free_y")

ggsave("NE_pie.png", path="figures/", bg = "white", width = 6, height = 4, dpi = 300)
  
PNW <- ggplot(propdat_PNW, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "mm")) + scale_fill_manual(values=c("darkred","orange","darkgreen","magenta","darkgrey")) + facet_wrap(~habitat, scales = "free_y")

ggsave("PNW_pie.png", path="figures/", bg = "white", width = 6, height = 4, dpi = 300)
  
PRAI <- ggplot(propdat_PRAI, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "mm")) + scale_fill_manual(values=c("darkred","orange","darkgreen","magenta","darkgrey")) + facet_wrap(~habitat, scales = "free_y")

ggsave("PRAI_pie.png", path="figures/", bg = "white", width = 6, height = 4, dpi = 300)
  
SE <- ggplot(propdat_SE, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "mm")) + scale_fill_manual(values=c("darkred","orange","darkgreen","magenta","darkgrey")) + facet_wrap(~habitat, scales = "free_y")

ggsave("SE_pie.png", path="figures/", bg = "white", width = 6, height = 4, dpi = 300)
  
SROC <- ggplot(propdat_SROC, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "mm")) + scale_fill_manual(values=c("darkred","orange","darkgreen","magenta","darkgrey")) + facet_wrap(~habitat, scales = "free_y")

ggsave("SROC_pie.png", path="figures/", bg = "white", width = 6, height = 4, dpi = 300)
  
TAIG <- ggplot(propdat_TAIG, aes(x = 1, y = prop, fill = funct)) + geom_bar(stat = "identity", width = 1) + coord_polar(theta = "y") + theme_void() + theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "mm")) + scale_fill_manual(values=c("darkred","orange","darkgreen","magenta","darkgrey")) + facet_wrap(~habitat, scales = "free_y")

ggsave("TAIG_pie.png", path="figures/", bg = "white", width = 6, height = 4, dpi = 300)

eco_pie <- ggarrange(DSW, GRTL, NE, PNW, PRAI, SE, SROC, TAIG, labels = c("DSW", "GRTL", "NE", "PNW", "PRAI", "SE", "SROC", "TAIG"), ncol = 4, nrow = 2, common.legend = TRUE, legend="bottom") + bgcolor("white")
ggsave("eco_pie.png", path="figures/", width=8, height =6)

ggplot(propdat, aes(x="", y=prop, fill=funct)) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="top") +
  facet_wrap(vars(habitat))  +
  facet_wrap(vars(eco.region))

ggplot(propdat, aes(fill = funct, x = habitat, y = prop)) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme(legend.position="top") +
  facet_wrap(vars(eco.region)) +
  ylim(0,100)

# map ----
# Get US states map data
us_states <- map_data("state")

# States you want to highlight
highlight_states <- c("oregon", "alaska", "arizona", 
                      "kansas", "minnesota", "colorado", 
                      "massachusetts", "florida")

# Add a column for highlighting
us_states <- us_states %>%
  mutate(highlight = ifelse(region %in% highlight_states, "highlight", "normal"))

ggplot(us_states, aes(long, lat, group = group, fill = highlight)) +
  geom_polygon(color = "black") +
  coord_quickmap() +                          # use coord_quickmap for geographic data
  scale_fill_manual(values = c("highlight" = "steelblue", "normal" = "gray90")) +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Functional group richness by eco-region")


states_to_highlight <- data.frame(
  state = c("OR","AK","AZ","KS","MN","CO","MA","FL"),  # use state abbreviations
  highlight = "yes"
)

blankmap <- plot_usmap(regions = "states", data = states_to_highlight, values = "highlight") +
  scale_fill_manual(values = c("yes" = "steelblue"), na.value = "gray90") +
  theme(legend.position = "none")
ggsave("blankmap.png", path="figures/", width=8, height =6)

# combine
# Example pie data for your states
# Suppose each state has proportions of A, B, and C
pie_data <- data.frame(
  state = c("OR","AK","AZ","KS","MN","CO","MA","FL"),
  A = c(0.2, 0.6, 0.3, 0.4, 0.5, 0.3, 0.7, 0.2),
  B = c(0.5, 0.2, 0.4, 0.4, 0.2, 0.5, 0.2, 0.4),
  C = c(0.3, 0.2, 0.3, 0.2, 0.3, 0.2, 0.1, 0.4)
)

if (!require("pacman")) install.packages("pacman")
p_load(ggplot2, ggtree, dplyr, tidyr, sp, maps, pipeR, grid, XML, gtable)

getLabelPoint <- function(county) {Polygon(county[c('long', 'lat')])@labpt}


# Get state centroids from usmap
centroids <- by(us_states, us_states$subregion, getLabelPoint)     # Returns list
centroids <- do.call("rbind.data.frame", centroids)  # Convert to Data Frame
names(centroids) <- c('long', 'lat')                 # Appropriate Header

centroids <- usmap::us_map(regions = "states") %>%
  group_by(full) %>%
  summarise(x = mean(range(x)), y = mean(range(y)))

centroids$geom

# Join centroids with pie_data
pie_locs <- pie_data %>%
  left_join(centroids, by = c("state" = "full"))

# Base map with highlighted states
map <- plot_usmap(regions = "states",
                  data = data.frame(state = pie_data$state, highlight="yes"),
                  values = "highlight") +
  scale_fill_manual(values = c("yes" = "gray90"), na.value = "gray95") +
  theme(legend.position = "none")

# Add pies
map +
  geom_scatterpie(data = pie_locs,
                  aes(x = x, y = y, r = 1),   # r controls pie radius
                  cols = c("A","B","C")) +
  coord_fixed() +
  scale_fill_manual(values = c("A"="steelblue","B"="orange","C"="forestgreen")) +
  labs(title = "US States with Pie Charts")


# from github ex ----
library(ggplot2)
library(dplyr)
library(maps)
library(grid)
library(gtable)
library(purrr)
if (!require("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
library(ggtree)

# Get US map polygons
us_map <- map_data("state")

# Define states of interest
states_of_interest <- c("oregon", "arizona", "colorado",
 "kansas", "minnesota", "massachusetts",  "florida", "alaska")

# Filter map
us_map_sub <- us_map %>% filter(region %in% states_of_interest)

# Compute centroids for those states
get_centroid <- function(df) {
  data.frame(long = mean(range(df$long)), 
             lat = mean(range(df$lat)))
}

centroids <- us_map_sub %>%
  group_by(region) %>%
  group_modify(~ get_centroid(.x)) %>%
  ungroup()

centroids <- bind_rows(centroids,
                       data.frame(region = "alaska",
                                  long = -150,
                                  lat  = 64))

# Merge fungal proportions with centroid coordinates
propdat2 <- propdat %>%
  mutate(state = tolower(state.name[match(state, state.abb)]))

propdat_with_coords <- propdat2 %>%
  mutate(state = tolower(state)) %>%
  left_join(centroids, by = c("state" = "region"))

propdat_with_coords <- propdat_with_coords %>%
  mutate(long_plot = ifelse(state=="alaska", long + 80, long),
         lat_plot  = ifelse(state=="alaska", lat - 30, lat))

# Split data by state
pie_data_list <- split(propdat_with_coords, propdat_with_coords$state)

pies <- map(pie_data_list, function(df) {
  ggplot(df, aes(x = 1, y = prop, fill = funct)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = unit(c(0,0,0,0), "mm")) +
    scale_fill_manual(values=c("darkred","orange","darkgreen","magenta","darkgrey"))
})

# Extract legend from one pie
p_leg <- ggplot(pie_data_list[[1]], aes(x=1, y=prop, fill=funct)) +
  geom_bar(stat="identity", width=1) +
  coord_polar(theta="y") +
  theme_void() + theme(legend.position = "right")

leg_grob <- gtable_filter(ggplot_gtable(ggplot_build(p_leg)), "guide-box")

##retry----
# Compute a scaling factor based on the range of the main map
x_range <- diff(range(us_map$long))
y_range <- diff(range(us_map$lat))
base_size <- 0.02  # fraction of map width/height

pie_overlays <- lapply(names(pies), function(st) {
  dat <- pie_data_list[[st]][1, ]
  size <- ifelse(st == "alaska", 6, 2)
  annotation_custom(
    grob = ggplotGrob(pies[[st]] +
    theme(plot.title = element_text(hjust = 0.0),
    legend.text = element_text(size = 12),
    legend.spacing.y = unit(0.0, 'cm'),
    legend.box.spacing = unit(-2, "pt"))+ facet_wrap(~habitat, strip.position = "bottom") + scale_fill_manual(values=c("darkred","orange","darkgreen","magenta","darkgrey"))),
    xmin = dat$long - size, xmax = dat$long + size,
    ymin = dat$lat - size, ymax = dat$lat + size
  )
})

# Extract legend from one pie
p_leg <- ggplot(pie_data_list[[1]], aes(x=1, y=prop, fill=funct)) +
  geom_bar(stat="identity", width=1) +
  coord_polar(theta="y") +
  theme_void() +
  scale_fill_manual(values=c("darkred","orange","darkgreen","magenta","darkgrey")) +
  theme(legend.position=c(0.75,0.5), legend.title=element_blank(), strip.text = element_text(size = 8))

extractLegend <- function(gg) {
  grobs <- ggplot_gtable(ggplot_build(gg))
  foo <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
  grobs$grobs[[foo]]
}

# Extract wanted legend
wantedLegend <- extractLegend(p_leg)

p <- ggplot(us_map, aes(long, lat, group=group)) +
  geom_polygon(fill="white", color="black") +
  coord_fixed(1.3) +
  theme_minimal() +
  pie_overlays + # add list of grobs
  wantedLegend
p
ggsave("us_fungi_pies_wh.jpg", p, path="figures/",
       width = 12, height = 8, units = "in", dpi = 300)

p2 <- p + coord_fixed(1.3, xlim = c(-160, -65), ylim = c(20, 70))
p2

ggsave("us_fungi_pies_wAK.jpg", p2, path="figures/",
       width = 12, height = 8, units = "in", dpi = 300)


library(cowplot)

map_plot <- ggplot(us_map, aes(long, lat, group=group)) +
  geom_polygon(fill="white", color="black") +
  coord_fixed(1.3, xlim = c(-180, -65), ylim = c(20, 70)) +
  theme_minimal()

p3 <- ggdraw(map_plot)

# Overlay pies using relative coordinates (0-1)
for (st in names(pies)) {
  dat <- pie_data_list[[st]][1, ]
  
  # rescale long/lat to [0,1] relative to xlim/ylim
  x <- (dat$long_plot - (-180)) / ( -65 - (-180))
  y <- (dat$lat_plot  - 20) / (70 - 20)
  
  size <- ifelse(st=="alaska", 0.05, 0.08)  # fraction of plot
  p3 <- p3 + draw_plot(pies[[st]], x=x, y=y, width=size, height=size)
}
p3
