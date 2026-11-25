
library(brms)
library(posterior)   # for as_draws_df
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(ggridges)
library(patchwork)
library(corrplot)

theme_set(theme_minimal())

options(mc.cores = max(1, parallel::detectCores() - 1))

# 1) paths, labels, etc ----
results_path  <- "results_plots/exclusive1"
models_path   <- file.path(results_path, "models")

# pretty labels
response_labels <- c(
  pyrorich="Pyrophile Richness", pyrodiv="Pyrophile Diversity",
  fungrich="Total Fungal Richness", fungdiv="Total Fungal Diversity",
  ppdiv1="Pathogen Diversity", pprich1="Pathogen Richness",
  allsapsdiv1="Saprotroph Diversity", allsapsrich1="Saprotroph Richness",
  ectorich="Ectomycorrhizal Richness", ectodiv="Ectomycorrhizal Diversity"
)

# pretty labels
param_labels <- c(
  moisture_c="Moisture", pH_c="pH", CN_c="C:N Ratio", elevation_c="Elevation",
  tmean_Annual_c="Mean Annual Temp", ppt_mean_Annual_c="Mean Annual Precip",
  compartmentsoil="Compartment: Soil",
  habitatConiferous="Habitat: Coniferous", habitatGrassland="Habitat: Grassland", habitatOak="Habitat: Oak",
  fireBurned="Fire: Burned", fireUnburned="Fire: Unburned",
  seasonSpring="Season: Spring", seasonSummer="Season: Summer",
  seasonFall="Season: Fall", seasonWinter="Season: Winter", year="Year"
)

# parameter groups (for selective plotting)
env_params       <- c("moisture_c","pH_c","CN_c","elevation_c","tmean_Annual_c","ppt_mean_Annual_c")
treatment_params <- c("compartmentsoil","habitatConiferous","habitatGrassland","habitatOak","fireBurned","fireUnburned")
temporal_params  <- c("seasonSpring","seasonSummer","seasonFall","seasonWinter","year")

# 2) read & prepare data ----

dat <- read.csv("processed data/final_data.csv", header = TRUE)

dat <- dat %>%
  mutate(
    moisture_c        = scale(moisture)[,1],
    pH_c              = scale(pH)[,1],
    CN_c              = scale(CN)[,1],
    elevation_c       = scale(elevation)[,1],
    tmean_Annual_c    = scale(tmean_Annual)[,1],
    ppt_mean_Annual_c = scale(ppt_mean_Annual)[,1]
  )

head(dat, 2)

# 3) load fitted models ----

model_files <- list.files(models_path, pattern = "_model\\.rds$", full.names = TRUE)

models <- model_files |>
  setNames(nm = str_remove(basename(model_files), "_model\\.rds$")) |>
  map(readRDS)

responses_all <- names(models)
responses_all

## check----
models$pyrorich
models$fungrich


# 4) correlation plot  ----

diversity_cols <- names(response_labels)
env_cols <- c("moisture","CN","pH","elevation","tmean_Annual","ppt_mean_Annual")

cor_vars  <- intersect(colnames(dat), c(diversity_cols, env_cols))
cor_data  <- dat |> dplyr::select(all_of(cor_vars)) |> dplyr::select(where(is.numeric))
cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

# pretty rename
to_rename <- intersect(rownames(cor_matrix), diversity_cols)
colnames(cor_matrix)[colnames(cor_matrix) %in% to_rename] <- response_labels[to_rename]
rownames(cor_matrix)[rownames(cor_matrix) %in% to_rename] <- response_labels[to_rename]

pdf(file.path(results_path, "variable_correlations.pdf"), width = 12, height = 10)
corrplot(cor_matrix, method="color", type="upper", order="hclust",
         tl.cex=0.7, tl.srt=45, tl.col="black",
         col=colorRampPalette(c("#8B0000","white","#00008B"))(100),
         addCoef.col="black", number.cex=0.5, diag=FALSE,
         title="Correlations: Response and Predictor Variables",
         mar=c(0,0,2,0))
dev.off()


# 5) PLOTTING FUNCTIONS ----

## fixed-effect summaries- single model ----

summarize_fixed <- function(fit, response_name, keep_params = NULL, drop_params = NULL) {
  draws <- as_draws_df(fit)
  fixed_cols <- grep("^b_(?!Intercept)", names(draws), value = TRUE, perl = TRUE)
  
  if (!is.null(keep_params)) {
    fixed_cols <- fixed_cols[vapply(fixed_cols, function(x) {
      p <- sub("^b_", "", x); p %in% keep_params
    }, logical(1))]
  }
  if (!is.null(drop_params)) {
    fixed_cols <- fixed_cols[!vapply(fixed_cols, function(x) {
      p <- sub("^b_", "", x); p %in% drop_params
    }, logical(1))]
  }
  if (length(fixed_cols) == 0) return(tibble())
  
  draws |>
    dplyr::select(all_of(fixed_cols)) |>
    pivot_longer(everything(), names_to = "param_raw", values_to = "value") |>
    mutate(param = sub("^b_", "", param_raw)) |>
    group_by(param) |>
    summarise(
      mean  = mean(value),
      median= median(value),
      q2.5  = quantile(value, 0.025),
      q5    = quantile(value, 0.05),
      q95   = quantile(value, 0.95),
      q97.5 = quantile(value, 0.975),
      prob_positive = mean(value > 0),
      prob_negative = mean(value < 0),
      .groups = "drop"
    ) |>
    mutate(
      response        = response_name,
      prob_direction  = pmax(prob_positive, prob_negative),
      significant_95  = (q2.5 > 0 | q97.5 < 0),
      significant_90  = (q5 > 0  | q95  < 0),
      param_label     = ifelse(param %in% names(param_labels),
                               param_labels[param],
                               str_to_title(gsub("_"," ",param))),
      response_label  = ifelse(response %in% names(response_labels),
                               response_labels[response], response)
    )
}

## summarize multiple models ----

summarize_models <- function(models, responses = names(models),
                             keep_params = NULL, drop_params = NULL) {
  map_dfr(responses, function(rsp) {
    if (!rsp %in% names(models)) return(tibble())
    summarize_fixed(models[[rsp]], rsp, keep_params, drop_params)
  }) %>%# put variables in correct order
    mutate(response_label = factor(response_label, levels = rev(desired_order)))
}


## forest panels by param group ----

plot_fixed_group <- function(fixed_df, title = NULL, ncol = 3) {
  # if (nrow(fixed_df) == 0) return(ggplot() + theme_void())
  
  effect_colors <- c(
    "Positive (95% CI)"="#00008B", "Positive (90% CI)"="#4169E1",
    "Negative (95% CI)"="#8B0000", "Negative (90% CI)"="#CD5C5C",
    "Not significant"="#808080"
  )
  
  df <- fixed_df |>
    mutate(
      sig_category = case_when(
        significant_95 & mean > 0 ~ "Positive (95% CI)",
        significant_95 & mean < 0 ~ "Negative (95% CI)",
        significant_90 & mean > 0 ~ "Positive (90% CI)",
        significant_90 & mean < 0 ~ "Negative (90% CI)",
        TRUE ~ "Not significant"
      ),
      sig_category = factor(sig_category, levels = names(effect_colors))
    )
  
  ggplot(df, aes(x = mean, y = response_label)) +
    geom_linerange(aes(xmin = q2.5, xmax = q97.5), size = 0.5, alpha = 0.5) +
    geom_linerange(aes(xmin = q5, xmax = q95, color = sig_category), size = 1.1, alpha = 0.9) +
    geom_point(aes(color = sig_category), size = 2.2) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    facet_wrap(~ param_label, scales = "free_x", ncol = ncol) +
    scale_color_manual(values = effect_colors, name = "Effect") +
    labs(x = "Effect size", y = NULL, title = title) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "bottom"
    )
}


## heatmap summary of fixed effects ----

plot_fixed_heatmap <- function(fixed_df) {
  if (nrow(fixed_df) == 0) return(ggplot() + theme_void())
  
  hm <- fixed_df |>
    mutate(
      effect_strength = case_when(
        significant_95 ~ mean,
        significant_90 ~ 0.7 * mean,
        TRUE           ~ 0.3 * mean
      ),
      param_label = factor(param_label)
    )
  
  lim <- max(abs(hm$mean))
  ggplot(hm, aes(x = param_label, y = response_label)) +
    geom_tile(aes(fill = effect_strength), color = "white", size = 0.4) +
    scale_fill_gradient2(low="#8B0000", mid="white", high="#00008B",
                         midpoint=0, limits=c(-lim, lim), name="Effect\nSize") +
    labs(title="Effect Summary Heatmap", x=NULL, y=NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle=45, hjust=1, size=9),
      axis.text.y = element_text(size=9),
      panel.grid  = element_blank(),
      legend.position = "right"
    )
}



# 6) build prediction grid ----

build_pred_grid <- function(data) {
  comp_lv <- if(is.factor(data$compartment)) levels(data$compartment) else sort(unique(na.omit(data$compartment)))
  hab_lv  <- if(is.factor(data$habitat))     levels(data$habitat)     else sort(unique(na.omit(data$habitat)))
  fire_lv <- if(is.factor(data$fire))        levels(data$fire)        else sort(unique(na.omit(data$fire)))
  season1 <- if(is.factor(data$season))      levels(data$season)[1]   else sort(unique(na.omit(data$season)))[1]
  year_md <- median(data$year, na.rm = TRUE)
  
  grid <- expand.grid(
    compartment        = comp_lv,
    habitat            = hab_lv,
    fire               = fire_lv,
    moisture_c         = 0, pH_c = 0, CN_c = 0, elevation_c = 0,
    tmean_Annual_c     = 0, ppt_mean_Annual_c = 0,
    season             = season1,
    year               = year_md,
    stringsAsFactors   = FALSE
  )
  
  # Align factor levels if needed
  as_model_factor <- function(x, ref) if(is.factor(ref)) factor(x, levels = levels(ref)) else x
  grid$compartment <- as_model_factor(grid$compartment, data$compartment)
  grid$habitat     <- as_model_factor(grid$habitat,     data$habitat)
  grid$fire        <- as_model_factor(grid$fire,        data$fire)
  grid$season      <- as_model_factor(grid$season,      data$season)
  grid
}

# 7) posterior draws & summaries for a single model ----

posterior_draws_grid <- function(fit, grid) {
  # Matrix: iterations x nrow(grid)
  posterior_epred(fit, newdata = grid, re_formula = NA)
}

posterior_summary_grid <- function(draws_mat, grid, response_name) {
  tibble(
    mean   = apply(draws_mat, 2, mean),
    median = apply(draws_mat, 2, median),
    q2.5   = apply(draws_mat, 2, quantile, 0.025),
    q10    = apply(draws_mat, 2, quantile, 0.10),
    q90    = apply(draws_mat, 2, quantile, 0.90),
    q97.5  = apply(draws_mat, 2, quantile, 0.975)
  ) |>
    bind_cols(grid) |>
    mutate(response = response_name,
           response_label = ifelse(response_name %in% names(response_labels),
                                   response_labels[response_name], response_name))
}

# 8) contrast probabilities from posterior draws ----

contrast_probs_from_draws <- function(draws_mat, grid, response_name) {
  # Ensure expected level names exist or adapt here
  if (!all(c("litter","soil") %in% as.character(grid$compartment))) {
    stop("contrast_probs_from_draws(): expected 'litter' and 'soil' in compartment levels.")
  }
  if (!all(c("Burned","Unburned") %in% as.character(grid$fire))) {
    stop("contrast_probs_from_draws(): expected 'Burned' and 'Unburned' in fire levels.")
  }
  
  habitats <- unique(grid$habitat)
  res_list <- vector("list", length(habitats))
  
  for (i in seq_along(habitats)) {
    h <- habitats[i]
    idx_lb <- which(grid$compartment=="litter" & grid$habitat==h & grid$fire=="Burned")
    idx_lu <- which(grid$compartment=="litter" & grid$habitat==h & grid$fire=="Unburned")
    idx_sb <- which(grid$compartment=="soil"   & grid$habitat==h & grid$fire=="Burned")
    idx_su <- which(grid$compartment=="soil"   & grid$habitat==h & grid$fire=="Unburned")
    
    if (length(idx_lb)==1 && length(idx_lu)==1 && length(idx_sb)==1 && length(idx_su)==1) {
      # Vectorized differences across posterior draws (rows)
      fire_diff_litter <- draws_mat[, idx_lb] - draws_mat[, idx_lu]
      fire_diff_soil   <- draws_mat[, idx_sb] - draws_mat[, idx_su]
      fire_diff_avg    <- (fire_diff_litter + fire_diff_soil)/2
      
      comp_diff_burned   <- draws_mat[, idx_sb] - draws_mat[, idx_lb]
      comp_diff_unburned <- draws_mat[, idx_su] - draws_mat[, idx_lu]
      comp_diff_avg      <- (comp_diff_burned + comp_diff_unburned)/2
      
      interaction <- comp_diff_burned - comp_diff_unburned
      
      # Two-tailed "probability of a difference" as max(P>0, P<0)
      p_fire <- mean(fire_diff_avg > 0);        prob_fire_diff <- max(p_fire, 1 - p_fire)
      p_comp <- mean(comp_diff_avg > 0);        prob_comp_diff <- max(p_comp, 1 - p_comp)
      p_int  <- mean(interaction   > 0);        prob_interaction_diff <- max(p_int, 1 - p_int)
      
      res_list[[i]] <- tibble(
        response       = response_name,
        response_label = ifelse(response_name %in% names(response_labels),
                                response_labels[response_name], response_name),
        habitat        = as.character(h),
        fire_effect        = mean(fire_diff_avg),
        comp_effect        = mean(comp_diff_avg),
        interaction_effect = mean(interaction),
        fire_direction = ifelse(mean(fire_diff_avg) > 0, "positive", "negative"),
        comp_direction = ifelse(mean(comp_diff_avg) > 0, "positive", "negative"),
        int_direction  = ifelse(mean(interaction)   > 0, "positive", "negative"),
        prob_fire_diff        = prob_fire_diff,
        prob_comp_diff        = prob_comp_diff,
        prob_interaction_diff = prob_interaction_diff
      )
    } else {
      res_list[[i]] <- tibble()
    }
  }
  bind_rows(res_list)
}

# 9) build predictions + probabilities for a subset of responses ----

compute_predictions_and_probs <- function(models, data, responses) {
  grid <- build_pred_grid(data)
  
  # For each response, compute posterior draws, summaries, and contrast probs
  out <- map(responses, function(rsp) {
    fit <- models[[rsp]]
    if (is.null(fit)) return(NULL)
    
    draws_mat <- posterior_draws_grid(fit, grid)
    pred_sum  <- posterior_summary_grid(draws_mat, grid, rsp)
    probs     <- contrast_probs_from_draws(draws_mat, grid, rsp)
    
    list(response = rsp, grid = grid, draws = draws_mat, pred = pred_sum, probs = probs)
  })
  
  # Combine
  preds <- bind_rows(map(out, "pred"))
  probs <- bind_rows(map(out, "probs"))
  
  list(grid = grid, predictions = preds, probabilities = probs, per_response = out)
}

# 10) plot function: Fire × Compartment × Habitat with probability text ----

plot_fire_comp_interaction <- function(pred_df, probs_df, out_file = NULL, width=10, height=8) {
  # Text annotations per response × habitat
  prob_text <- probs_df |>
    mutate(
      fire_text = paste0("Fire: ", sprintf("%.0f%%", prob_fire_diff*100)),
      comp_text = paste0("Comp: ", sprintf("%.0f%%", prob_comp_diff*100)),
      int_text  = paste0("Int: ",  sprintf("%.0f%%", prob_interaction_diff*100)),
      fire_color = ifelse(prob_fire_diff >= 0.90 & fire_direction=="positive", "blue",
                          ifelse(prob_fire_diff >= 0.90 & fire_direction=="negative", "red", "black")),
      comp_color = ifelse(prob_comp_diff >= 0.90 & comp_direction=="positive", "blue",
                          ifelse(prob_comp_diff >= 0.90 & comp_direction=="negative", "red", "black")),
      int_color  = ifelse(prob_interaction_diff >= 0.90 & int_direction=="positive", "blue",
                          ifelse(prob_interaction_diff >= 0.90 & int_direction=="negative", "red", "black"))
    )
  
  p <- ggplot(pred_df, aes(x = fire, y = mean)) +
    geom_errorbar(aes(ymin = q10, ymax = q90, color = compartment),
                  width = 0.15, alpha = 0.6, position = position_dodge(0.3), size = 0.8) +
    geom_point(aes(color = compartment), size = 2.8, position = position_dodge(0.3)) +
    geom_text(data = prob_text, aes(x = 1.2, y = Inf, label = fire_text, color = fire_color),
              hjust = 0, vjust = 1.5, size = 3, inherit.aes = FALSE) +
    geom_text(data = prob_text, aes(x = 1.2, y = Inf, label = comp_text, color = comp_color),
              hjust = 0, vjust = 3.0, size = 3, inherit.aes = FALSE) +
    geom_text(data = prob_text, aes(x = 1.2, y = Inf, label = int_text,  color = int_color),
              hjust = 0, vjust = 4.5, size = 3, inherit.aes = FALSE) +
    facet_grid(response_label ~ habitat, scales = "free_y", switch = "y") +
    scale_color_manual(values = c(litter = "#E69F00", soil = "#56B4E9",
                                  blue = "blue", red = "red", black = "black"),
                       breaks = c("litter","soil"), name = "Compartment") +
    labs(
      # title = "Compartment × Habitat × Fire Interactions",
      # subtitle = "Probabilities: blue = positive ≥90%, red = negative ≥90%",
      x = "Fire Status", y = NULL
    ) +
    theme_minimal() +
    theme(
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8, face = "bold"),
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.placement = "outside",
      legend.position = "bottom"
    )
  
    ggsave(file.path(results_path, out_file), p, width = width, height = height)
  p
}


# 11) heatmap ----

plot_probability_heatmap <- function(probs_df, out_file = NULL, width=8, height=6) {
  df <- probs_df |>
    pivot_longer(cols = c(prob_fire_diff, prob_comp_diff, prob_interaction_diff),
                 names_to = "prob_type", values_to = "probability") |>
    pivot_longer(cols = c(fire_direction, comp_direction, int_direction),
                 names_to = "dir_type", values_to = "direction") |>
    filter((prob_type == "prob_fire_diff"        & dir_type == "fire_direction") |
             (prob_type == "prob_comp_diff"        & dir_type == "comp_direction") |
             (prob_type == "prob_interaction_diff" & dir_type == "int_direction")) |>
    mutate(
      column_label = recode(prob_type,
                            prob_fire_diff="Fire",prob_comp_diff="Compartment",prob_interaction_diff="Fire × Comp"),
      sig_level = case_when(
        probability > 0.95 & direction=="positive" ~ "Pos >95%",
        probability >=0.90 & direction=="positive" ~ "Pos 90-95%",
        probability > 0.95 & direction=="negative" ~ "Neg >95%",
        probability >=0.90 & direction=="negative" ~ "Neg 90-95%",
        TRUE ~ "NS"
      ),
      sig_level = factor(sig_level, levels = c("Pos >95%","Pos 90-95%","NS","Neg 90-95%","Neg >95%"))
    )
  
  p <- ggplot(df, aes(x = column_label, y = response_label, fill = sig_level)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.0f%%", probability*100)), size = 3) +
    facet_grid(~ habitat, scales = "free_x", space = "free_x") +
    scale_fill_manual(
      values = c("Pos >95%"="#08519c","Pos 90-95%"="#6baed6","NS"="#d9d9d9","Neg 90-95%"="#fb6a4a","Neg >95%"="#a50f15"),
      name = "Effect Direction\n& Significance", drop = FALSE
    ) +
    labs(x=NULL, y=NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle=45, hjust=1, size=9),
      axis.text.y = element_text(size=9),
      strip.text  = element_text(face="bold", size=10),
      legend.position = "bottom",
      panel.grid = element_blank(),
      plot.subtitle = element_text(size=10)
    )
  
  if (!is.null(out_file)) {
    dir.create(file.path(results_path), showWarnings = FALSE)
    ggsave(file.path(results_path, out_file), p, width = width, height = height)
  }
  p
}


# 12) * Run workflow for chosen subset ----

# responses_to_plot <- c("fungdiv","fungrich")  # 
responses_to_plot <- c("fungrich", "pprich1","allsapsrich1", "ectorich","pyrorich")

# get predictions and probabilities 
pred_prob <- compute_predictions_and_probs(models, dat, responses_to_plot)

head(pred_prob$predictions)
head(pred_prob$probabilities)

# haley: Change order to saprotrophs, pathogens, ectomycorrhizae, pyrophiles, all fungi
desired_order <- c(  "Saprotroph Richness", "Pathogen Richness","Ectomycorrhizal Richness","Pyrophile Richness","Total Fungal Richness")
  
pred_prob$predictions$response_label  <- factor(
  pred_prob$predictions$response_label,   levels = desired_order
)
pred_prob$probabilities$response_label <- factor(
  pred_prob$probabilities$response_label,  levels = desired_order
)
# Fire × Compartment × Habitat plot with posterior-based probabilities
p_fc <- plot_fire_comp_interaction(
  pred_df = pred_prob$predictions,
  probs_df = pred_prob$probabilities,
  out_file = "group_pred.pdf",
  width = 10, height = 8
)

# heatmap
p_heat_prob <- plot_probability_heatmap(
  pred_prob$probabilities,
  out_file = "heatmap.pdf",
  width = 8, height = 6
)

# * PLOT fixed-effect summaries and plots ----
fixed_env       <- summarize_models(models, responses = responses_to_plot, keep_params = env_params)
fixed_treatment <- summarize_models(models, responses = responses_to_plot, keep_params = treatment_params)
fixed_temporal  <- summarize_models(models, responses = responses_to_plot, keep_params = temporal_params)

p_env       <- plot_fixed_group(fixed_env) # ,       "Environmental Variables"
p_treatment <- plot_fixed_group(fixed_treatment,ncol=2) # , "Habitat / Fire / Compartment Variables
p_temporal  <- plot_fixed_group(fixed_temporal) # ,  "Temporal Variables"

combined <- p_env / p_treatment / p_temporal + plot_layout(heights = c(2,1.5,1))
ggsave(file.path(results_path, "param_intervals.pdf"), combined, width=12, height=12)


# plot subsets - easier to view
p_env 
ggsave(file.path(results_path, "env.pdf"), width=7, height=5)

p_temporal
ggsave(file.path(results_path, "temporal.pdf"), width=7, height=5)

p_treatment
ggsave(file.path(results_path, "fire_comp_habitat.pdf"), width=10, height=8)


# 16) run workflow for ALL vairables - div and rich ----

# Choose which responses to include in predictions/contrasts
ALL_responses_to_plot <- c("fungrich","fungdiv","ppdiv1", "pprich1","allsapsdiv1","allsapsrich1","ectorich","ectodiv","ectorich","pyrorich","pyrodiv")

# Compute predictions and probabilities (fully Bayesian)
pred_prob_ALL <- compute_predictions_and_probs(models, dat, ALL_responses_to_plot)

# Fire × Compartment × Habitat plot with posterior-based probabilities
p_fc <- plot_fire_comp_interaction(
  pred_df = pred_prob_ALL$predictions,
  probs_df = pred_prob_ALL$probabilities,
  out_file = "group_pred_ALL.pdf",
  width = 10, height = 15
)

# Heatmap of directional probabilities
p_heat_prob <- plot_probability_heatmap(
  pred_prob_ALL$probabilities,
  out_file = "heatmap_ALL.pdf",
  width = 8, height = 6
)

# Example fixed-effect summaries and plots
fixed_env       <- summarize_models(models, responses = ALL_responses_to_plot, keep_params = env_params)
fixed_treatment <- summarize_models(models, responses = ALL_responses_to_plot, keep_params = treatment_params)
fixed_temporal  <- summarize_models(models, responses = ALL_responses_to_plot, keep_params = temporal_params)

p_env       <- plot_fixed_group(na.omit(fixed_env),       "Environmental Variables")
p_treatment <- plot_fixed_group(na.omit(fixed_treatment), "Habitat / Fire / Compartment Variables", ncol=2)
p_temporal  <- plot_fixed_group(na.omit(fixed_temporal),  "Temporal Variables")

combined <- p_env / p_treatment / p_temporal + plot_layout(heights = c(2,1.5,1))
ggsave(file.path(results_path, "param_intervals_ALL2.pdf"), combined, width=14, height=16)
ggsave(file.path(results_path, "param_intervals_ALL2.jpg"), combined, width=10, height=14)

combined2 <- p_env / p_treatment + plot_layout(heights = c(2,1.5))
ggsave(file.path(results_path, "param_intervals2.jpg"), combined2, width=10, height=8)

# plot subsets - easier to view
p_env 
ggsave(file.path(results_path, "env_ALL.pdf"), width=8, height=7)

p_temporal
ggsave(file.path(results_path, "temporal_ALL.pdf"), width=7, height=5)

p_treatment
ggsave(file.path(results_path, "fire_comp_habitat_ALL.pdf"), width=10, height=8)

# compare raw inclusive vs exclusive ----

p1 <- ggplot(dat, aes(x = pprich1, y = pprich2)) +
  geom_abline(intercept=0, slope=1) + geom_point() + labs(title="pathogens")

p2 <- ggplot(dat, aes(x = allsapsrich1, y = allsapsrich2)) +
  geom_abline(intercept=0, slope=1) + geom_point() + labs(title="all sapro")

p1 | p2
ggsave(file.path(results_path, "1 vs 2 raw.pdf"), width=8, height=4)
