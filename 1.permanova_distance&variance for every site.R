library(vegan)
library(tidyverse)
library(ggpubr)
library(viridis)

# data ----

sites <- read.csv("./my_data/selected_sites.csv",
                  sep = ",")

cover_dat <- read.csv("./nutnet_data/full-cover-2025-12-09.csv",
                      sep = ";")

metadat <- read.csv("./nutnet_data/comb-by-plot-clim-soil-diversity-2025-12-09.csv",
                    sep = ",")

cover_df <- sites %>%
  left_join(cover_dat,by="site_code") %>%
  filter(trt %in% c("Control","NPK"))

meta_df <- sites %>%
  left_join(metadat,by="site_code") %>%
  distinct(site_code,block,plot, .keep_all = TRUE)


# analysis for all sites ----


# List for collecting final results
clasification_results <- list()


sites_list <- unique(cover_df$site_code)

variance_results <- list()

set.seed(123)
for (current_site in sites_list) {
  
  # 1. Filtrate sites, structure community matrix and collapse years
  sp_wide <- cover_df %>%
    filter(live == 1, site_code == current_site, year_trt > 3) %>%
    
    # Keep only relevant variables
    select(site_code, year_trt, block, plot, trt, Taxon, max_cover) %>% 
    pivot_wider(
      names_from = Taxon,
      values_from = max_cover,
      values_fill = 0,
      values_fn = sum
    )
  
  # Identify treatments (excluding Control)
  trts <- unique(sp_wide$trt[sp_wide$trt != "Control"])
  
  # 2. Nested loop to evaluate each trt against Control
  for (current_trt in trts) {
    
    # Isolate Control and treatment
    sp_subset <- sp_wide %>% filter(trt %in% c("Control", current_trt))
    
    # Separate metadata and community matrix
    metadat <- sp_subset[, 1:5]
    community <- sp_subset[, 6:ncol(sp_subset)]
    
    # Check for enough data for permutations (mínimo exigido)
    if (nrow(metadat) < 6) next 
    
    # Create grouping factor (Ignore year, collapse by block and trt)
    group_factor <- paste(metadat$block, metadat$trt, sep = "_")
    
    # Calculate global distance matrix for the subset
    distance <- vegdist(community, method = "bray")
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # AXIS 1: DIFERENTIATION (PERMANOVA)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # adonis2 checks for centroid difference, controlling for block design
    perm_mod <- adonis2(distance ~ trt, data = metadat,
                        strata = metadat$block, permutations = 999)
    p_permanova <- perm_mod$`Pr(>F)`[1]
    
    # Lógica booleana de diferenciación (alfa = 0.05)
    sig_distance <- ifelse(p_permanova < 0.05, TRUE, FALSE)
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # AXIS 2: DYNAMICS / VARIANCE DIFFERENCE (BETADISPER)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # A. Dispersion by trt
    dispersion_trt <- betadisper(distance, group = metadat$trt, bias.adjust = TRUE, sqrt.dist = TRUE)
    
    # B. Dispersion by block and trt
    dispersion_block_trt <- betadisper(distance, group = group_factor, bias.adjust = TRUE, sqrt.dist = TRUE)
    
    # ---> Extract individual variances based on block+trt
    temp_var_df <- data.frame(
      site_code = current_site,
      block = metadat$block,
      trt = metadat$trt,
      distance_to_centroid = dispersion_block_trt$distances
    )
    
    # Save on variance_results list
    variance_results[[paste(current_site, current_trt, sep = "_")]] <- temp_var_df
    
    
    # Remaining calculations use dispersion_trt:
    disp_test <- permutest(dispersion_trt, permutations = 999)
    p_permutest <- disp_test$tab$`Pr(>F)`[1]
    
    # Lógica booleana de varianza
    sig_variance <- ifelse(p_permutest < 0.05, TRUE, FALSE)
    
    # Calculate directionality of change
    dist_centroid_trt <- dispersion_trt$distances
    mean_var_ctrl <- mean(dist_centroid_trt[metadat$trt == "Control"])
    mean_var_trt  <- mean(dist_centroid_trt[metadat$trt == current_trt])
    
    diff_variance <- mean_var_trt - mean_var_ctrl
    
    # Determine variance difference ("Positive", "Negative", "Null")
    var_difference <- "Null"
    if (sig_variance && diff_variance > 0) var_difference <- "Positive"
    if (sig_variance && diff_variance < 0) var_difference <- "Negative"
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # EXTRACT DISTANCE BETWEEN CENTROIDS (OF TREATMENTS)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    coord_centroid <- dispersion_trt$centroids
    centroid_distance <- as.numeric(dist(coord_centroid[c("Control", current_trt), ]))
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # CLASIFICATION MATRIX
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    category <- NA_character_
    
    if (!sig_distance && var_difference == "Null")     category <- "= Mean, = Variance"
    if (!sig_distance && var_difference == "Positive") category <- "= Mean,  + Variance"
    if (!sig_distance && var_difference == "Negative") category <- "= Mean, - Variance"
    
    if (sig_distance && var_difference == "Null")      category <- "Δ Mean, = Variance"
    if (sig_distance && var_difference == "Positive")  category <- "Δ Mean, + Variance"
    if (sig_distance && var_difference == "Negative")  category <- "Δ Mean, - Variance"
    
    
    category_2 <- NA_character_
    
    if (!sig_distance && var_difference == "Null")     category_2 <- "eqmean_eqvar"
    if (!sig_distance && var_difference == "Positive") category_2 <- "eqmean_divergence"
    if (!sig_distance && var_difference == "Negative") category_2 <- "eqmean_convergence"
    
    if (sig_distance && var_difference == "Null")      category_2 <- "difmean_eqvar"
    if (sig_distance && var_difference == "Positive")  category_2 <- "difmean_divergence"
    if (sig_distance && var_difference == "Negative")  category_2 <- "difmean_convergence"
    
    # Compile iteration results
    iteration_result <- data.frame(
      site_code = current_site,
      trt = current_trt,
      centroid_distance = centroid_distance, # <--- NUEVA COLUMNA
      p_permanova = p_permanova,
      sig_distance = sig_distance,
      diff_variance = diff_variance,
      p_permutest = p_permutest,
      sig_variance = sig_variance,
      var_difference = var_difference,
      category_impact = category,
      category_2 = category_2
    )
    
    # Save on list
    name_list <- paste(current_site, current_trt, sep = "_")
    clasification_results[[name_list]] <- iteration_result
  }
}

# 3. Consolidate final dataframe
df_clasification_final <- bind_rows(clasification_results)

df_variances_final <- do.call(rbind, variance_results)
rownames(df_variances_final) <- NULL # Limpiar nombres de fila para mayor prolijidad


df_variances_final <- df_variances_final %>%
  group_by(site_code,trt,block) %>%
  summarise(variance = mean(distance_to_centroid))


 write.csv(df_clasification_final,
          "./results/centroid-dist_variance-diff.csv")


 write.csv(df_variances_final,"./results/variance_by_trt_block.csv")
