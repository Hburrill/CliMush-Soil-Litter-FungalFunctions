setwd("C:/Users/haley/University of Oregon Dropbox/Haley Burrill/CLIMUSH NSF-Neon/Analysis/SoilLitterFunctionPaper_HB/Analyses/output/figures")
guilds   <- c("Saprotrophs","Ectomycorrhizae","Pathogens","Pyrophiles")
habitats <- c("Conifer","Grassland","Oak")
fires    <- c("Burned","Unburned")
comps    <- c("litter","soil")

# Schematic predictions (0–1 scale), same logic as before
df <- expand.grid(
  guild = guilds,
  habitat = habitats,
  fire = fires,
  compartment = comps) %>%
  as_tibble() %>%
  mutate(  # build a predicted diversity value (pred) for every combination of guild × habitat × fire × compartment.
    hab_mult = case_when(
      guild == "Pathogens"        & habitat == "Oak"       ~ 1.00,
      guild == "Pathogens"        & habitat == "Grassland" ~ 0.65,
      guild == "Pathogens"        & habitat == "Conifer"   ~ 0.45,
      
      guild == "Saprotrophs"      & habitat == "Oak"       ~ 1.00,
      guild == "Saprotrophs"      & habitat == "Conifer"   ~ 0.75,
      guild == "Saprotrophs"      & habitat == "Grassland" ~ 0.55,
      
      guild == "Ectomycorrhizae"     & habitat == "Conifer"   ~ 1.00,
      guild == "Ectomycorrhizae"     & habitat == "Oak"       ~ 0.65,
      guild == "Ectomycorrhizae"     & habitat == "Grassland" ~ 0.15,
      
      guild == "Pyrophiles"            & habitat == "Conifer"   ~ 1.00,
      guild == "Pyrophiles"            & habitat == "Oak"       ~ 0.75,
      guild == "Pyrophiles"            & habitat == "Grassland" ~ 0.55,
      TRUE ~ 0.7
    ),
    fire_adj = case_when(
      guild %in% c("Pathogens","Saprotrophs","Ectomycorrhizae") & fire == "Unburned" ~ +0.20,
      guild %in% c("Pathogens","Saprotrophs","Ectomycorrhizae") & fire == "Burned"   ~ -0.05,
      guild == "Pyrophiles"                                                        & fire == "Burned"   ~ +0.25,
      guild == "Pyrophiles"                                                        & fire == "Unburned" ~ -0.05,
      TRUE ~ 0
    ),
    comp_adj = case_when(
      guild %in% c("Pathogens","Saprotrophs") & compartment == "litter" ~ +0.20,
      guild %in% c("Pathogens","Saprotrophs") & compartment == "soil"   ~ -0.05,
      
      guild == "Ectomycorrhizae"                         & compartment == "soil"   ~ +0.25,
      guild == "Ectomycorrhizae"                         & compartment == "litter" ~ -0.10,
      
      guild == "Pyrophiles" & fire == "Burned"             & compartment == "soil"   ~ +0.20,
      guild == "Pyrophiles" & fire == "Burned"             & compartment == "litter" ~ -0.05,
      guild == "Pyrophiles" & fire == "Unburned"           & compartment == "litter" ~ 0.10,
      guild == "Pyrophiles" & fire == "Unburned"           & compartment == "soil"   ~ -0.25,
      TRUE ~ 0
    ),
    pred = pmin(1, pmax(0, 0.35 + 0.45*hab_mult + fire_adj + comp_adj)),
    guild = factor(guild, levels = guilds),
    habitat = factor(habitat, levels = habitats),
    fire = factor(fire, levels = fires),
    compartment = factor(compartment, levels = comps)
  )

# apply absolute overrides ----
# Absolute overrides with wildcards: leave any key as NA to mean "all"
overrides_abs <- tibble::tribble(
  ~guild,         ~habitat,   ~fire,       ~compartment, ~abs,
  # EXAMPLES:
  "Pyrophiles",    NA,         "Burned",    NA,           0.90  # all Pyrophiles when Burned
  # ,"Pathogens",     "Oak",      NA,          NA,           0.80  # Pathogens in Oak (any fire/compartment)
  # , "Saprotrophs", "Grassland","Unburned", "soil",     0.0 # test with extreme value
)

# Helper to apply absolute overrides; NA = wildcard; later rows overwrite earlier ones
apply_abs_overrides <- function(dat, ov) {
  if (is.null(ov) || nrow(ov) == 0) return(dat)
  keys <- c("guild","habitat","fire","compartment")
  for (i in seq_len(nrow(ov))) {
    row <- ov[i, ]
    idx <- rep(TRUE, nrow(dat))
    for (k in keys) {
      v <- row[[k]][[1]]
      if (!is.na(v)) idx <- idx & as.character(dat[[k]]) == as.character(v)
    }
    if (!is.na(row$abs[[1]])) dat$pred[idx] <- row$abs[[1]]
  }
  dat$pred <- pmin(1, pmax(0, dat$pred))  # keep in [0,1]
  dat
}

# Apply overrides (last row wins)
df <- apply_abs_overrides(df, overrides_abs)


lab_fun <-  labeller(guild = label_value, habitat = label_value)

comp_cols <- c(litter = "#E69F00", soil = "#56B4E9")

df$guild <- factor(df$guild, levels=c("Saprotrophs","Pathogens","Ectomycorrhizae","Pyrophiles"))
p <- ggplot(df, aes(x = fire, y = pred, color = compartment, group = compartment)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.8) +
  # geom_line(position = position_dodge(width = 0.38), size = 1.2) +
  # geom_point(position = position_dodge(width = 0.38), size = 2.8) +
  scale_color_manual(values = comp_cols, name = "Compartment") +
  scale_y_continuous("Hypothesized Relative Richness", limits = c(0,1), breaks = c(0, .5, 1)) +
  scale_x_discrete("Fire status") +
  facet_grid(rows = vars(guild), cols = vars(habitat), switch = "y",
             labeller = lab_fun) +
  # theme_minimal(base_size = 10) #+
theme_minimal() +
theme(
  strip.placement = "outside",
  # strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", lineheight = 1.0),
  strip.text.x = element_text(face = "bold"),
  legend.position = "bottom",
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_line(color = "grey90"),
  panel.spacing.x = unit(0.4, "lines"),
  # panel.spacing.y = unit(0.25, "lines"),
  plot.title = element_text(face = "bold")
)

p
ggsave("hypothesesHB.jpg", width=6, height=6)
