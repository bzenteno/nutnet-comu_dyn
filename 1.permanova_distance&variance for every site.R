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


# Inicializar lista para acumular resultados finales
clasification_results <- list()

# Supuestos del entorno: 
# cover_df es tu base de datos larga.
# sites_list es tu vector con los códigos de los sitios.

sites_list <- unique(cover_df$site_code)

variance_results <- list()

for (current_site in sites_list) {
  
  # 1. Filtrar sitio, estructurar matriz comunitaria y colapsar años
  sp_wide <- cover_df %>%
    filter(live == 1, site_code == current_site, year_trt > 3) %>%
    # Asegurar que solo entren las variables relevantes
    select(site_code, year_trt, block, plot, trt, Taxon, max_cover) %>% 
    pivot_wider(
      names_from = Taxon,
      values_from = max_cover,
      values_fill = 0,
      values_fn = sum
    )
  
  # Identificar tratamientos a evaluar (excluyendo el Control)
  trts <- unique(sp_wide$trt[sp_wide$trt != "Control"])
  
  # 2. Bucle anidado para evaluar cada tratamiento contra el Control
  for (current_trt in trts) {
    
    # Aislar solo el Control y el tratamiento de interés
    sp_subset <- sp_wide %>% filter(trt %in% c("Control", current_trt))
    
    # Separar metadatos y matriz comunitaria
    metadat <- sp_subset[, 1:5]
    community <- sp_subset[, 6:ncol(sp_subset)]
    
    # Validar que haya suficientes datos para permutaciones (mínimo exigido)
    if (nrow(metadat) < 6) next 
    
    # Crear factor de agrupación estática (Ignora el año, colapsa por bloque y trt)
    group_factor <- paste(metadat$block, metadat$trt, sep = "_")
    
    # Calcular matriz de distancias global para el subconjunto
    distance <- vegdist(community, method = "bray")
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # EJE 1: DIFERENCIACIÓN (PERMANOVA)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # adonis2 evalúa si los centroides difieren, controlando por el diseño de bloques
    perm_mod <- adonis2(distance ~ trt, data = metadat,
                        strata = metadat$block, permutations = 999)
    p_permanova <- perm_mod$`Pr(>F)`[1]
    
    # Lógica booleana de diferenciación (alfa = 0.05)
    sig_distance <- ifelse(p_permanova < 0.05, TRUE, FALSE)
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # EJE 2: DINAMISMO / VARIANZA (BETADISPER)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # A. Dispersión por TRATAMIENTO (Para mantener tus estadísticas y categorías globales)
    dispersion_trt <- betadisper(distance, group = metadat$trt, bias.adjust = TRUE, sqrt.dist = TRUE)
    
    # B. Dispersión por BLOQUE_TRATAMIENTO (Para lo que tú quieres: centroide y varianza por bloque)
    dispersion_block_trt <- betadisper(distance, group = group_factor, bias.adjust = TRUE, sqrt.dist = TRUE)
    
    # ---> EXTRAER LAS VARIANZAS INDIVIDUALES BASADAS EN EL BLOQUE_TRATAMIENTO
    temp_var_df <- data.frame(
      site_code = current_site,
      block = metadat$block,
      trt = metadat$trt,
      distance_to_centroid = dispersion_block_trt$distances # <--- Usamos el objeto de bloques
    )
    
    # Guardar en la lista de resultados de varianzas
    variance_results[[paste(current_site, current_trt, sep = "_")]] <- temp_var_df
    
    
    # El resto de tus cálculos siguen usando 'dispersion_trt' para no romperse:
    disp_test <- permutest(dispersion_trt, permutations = 999)
    p_permutest <- disp_test$tab$`Pr(>F)`[1]
    
    # Lógica booleana de varianza
    sig_variance <- ifelse(p_permutest < 0.05, TRUE, FALSE)
    
    # Calcular la dirección del cambio usando las distancias del tratamiento global
    dist_centroid_trt <- dispersion_trt$distances
    mean_var_ctrl <- mean(dist_centroid_trt[metadat$trt == "Control"])
    mean_var_trt  <- mean(dist_centroid_trt[metadat$trt == current_trt])
    
    diff_variance <- mean_var_trt - mean_var_ctrl
    
    # Determinar el estado de la varianza ("Positiva", "Negativa", "Nula")
    var_difference <- "Null"
    if (sig_variance && diff_variance > 0) var_difference <- "Positive"
    if (sig_variance && diff_variance < 0) var_difference <- "Negative"
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # EXTRACCIÓN DE DISTANCIA ENTRE CENTROIDES (DE TRATAMIENTO)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # Esto ya NO dará error porque dispersion_trt sí contiene "Control" y tu current_trt
    coord_centroid <- dispersion_trt$centroids
    centroid_distance <- as.numeric(dist(coord_centroid[c("Control", current_trt), ]))
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # MATRIZ DE CLASIFICACIÓN ECOLÓGICA
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
    if (!sig_distance && var_difference == "Positive") category_2 <- "eqmean_posdifvar"
    if (!sig_distance && var_difference == "Negative") category_2 <- "eqmean_negdifvar"
    
    if (sig_distance && var_difference == "Null")      category_2 <- "difmean_eqvar"
    if (sig_distance && var_difference == "Positive")  category_2 <- "difmean_posdifvar"
    if (sig_distance && var_difference == "Negative")  category_2 <- "difmean_negdifvar"
    
    # Compilar los resultados de esta iteración
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
    
    # Almacenar en la lista
    name_list <- paste(current_site, current_trt, sep = "_")
    clasification_results[[name_list]] <- iteration_result
  }
}

# 3. Consolidar el dataframe final
df_clasification_final <- bind_rows(clasification_results)

df_variances_final <- do.call(rbind, variance_results)
rownames(df_variances_final) <- NULL # Limpiar nombres de fila para mayor prolijidad


df_variances_final <- df_variances_final %>%
  group_by(site_code,trt,block) %>%
  summarise(variance = mean(distance_to_centroid))


# write.csv(df_clasification_final,
  #        "./analysis/results/centroid-dist_variance-diff.csv")


# write.csv(df_variances_final,"./analysis/results/variance_by_trt_block.csv")
