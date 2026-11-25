# ==============================================================================
# ANÁLISIS DE REDUNDANCIA (RDA)
# Relación entre el síndrome de respuesta, las condiciones ambientales y bióticas
# basales y los procesos comunitarios (diferencias NPK − Control) en sitios NutNet
# ==============================================================================


# 0. LIBRERÍAS -----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(permute)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(scales)


# 1. Parameters ---------------------------------------------------
# decision filters

anio_min      <- 3     # se usan años de tratamiento > anio_min
bloques_max   <- 3     # se usan los bloques 1 a bloques_max de cada sitio
# (diseño balanceado, necesario para permutar sitios completos)
sitios_excl   <- c("sevi.us", 
                   "hero.uk", "kiny.au", "marc.ar"
)  # [documentar motivo]
sindrome_excl <- "eqmean_eqvar"  # sitios sin cambio composicional ni de dispersión
n_perm        <- 9999

# Variables de respuesta (matriz Y)
vars_Y <- c("d_richness_change", "d_evenness_change", "d_rank_change",
            "d_gains", "d_losses",
            "richness_diff", "evenness_diff", "rank_diff", "species_diff")

# Predictores candidatos (deben coincidir con los descriptos en métodos)
vars_X <- c(
  # Clima medio y variabilidad intraanual (WorldClim)
  "MAP_v2", "MAT_v2", "MAP_VAR_v2",
  "TEMP_VAR_v2",
  "AI",
  # Variabilidad interanual (NASA POWER)
  "avg_cvT",            # [agregar la variabilidad interanual de precipitación]
  # Suelo
  "pct_N", "ppm_P", "ppm_K", "pH",
  # Comunidad basal
  "richness", "SimpsonEvenness", "Evar","site_richness"
)


# 2. CARGA DE DATOS ------------------------------------------------------------
selected_sites <- read.csv("my_data/selected_sites.csv") %>% select(site_code)

metadata     <- read.csv("nutnet_data/comb-by-plot-2025-12-09.csv", sep = ";")
metadat_clim <- read.csv("nutnet_data/comb-by-plot-clim-soil-2025-12-09.csv",
                         na.strings = c("NULL", "", "N/A"))

rac_change <- read.csv("my_data/RAC_change_results.csv")[, -1]
rac_diff   <- read.csv("my_data/RAC_diff_results.csv")

sindromes   <- read.csv("results/centroid-dist_variance-diff.csv")[, -1]
suelo       <- read.csv("my_data/soil-data_summarised.csv")
comunidad   <- read.csv("my_data/richness_evennes.csv")
clima_inter <- read.csv("my_data/temp_anom.csv")


# 3. MATRIZ DE RESPUESTA (Y) ---------------------------------------------------

# 3.1 Diferencias NPK − Control en los cambios interanuales, promediadas por bloque
parcelas <- selected_sites %>%
  left_join(metadata %>% distinct(site_code, block, plot, trt), by = "site_code") %>%
  filter(trt %in% c("Control", "NPK"))

# Nota: year_trt2 es el segundo año de cada par; con year_trt2 > 3 se incluye
# el cambio entre los años 3 y 4. Si el criterio es "pares posteriores al año 3",
# el filtro debería ser year_trt2 > anio_min + 1. [verificar]
d_cambio <- parcelas %>%
  left_join(rac_change, by = c("site_code", "plot", "trt")) %>%
  filter(year_trt2 > anio_min+1) %>%
  select(site_code, block, year_trt2, trt,
         richness_change, evenness_change, rank_change, gains, losses) %>%
  pivot_wider(names_from  = trt,
              values_from = c(richness_change, evenness_change, rank_change,
                              gains, losses),
              values_fn   = mean) %>%
  transmute(site_code, block, year_trt2,
            d_richness_change = richness_change_NPK - richness_change_Control,
            d_evenness_change = evenness_change_NPK - evenness_change_Control,
            d_rank_change     = rank_change_NPK     - rank_change_Control,
            d_gains           = gains_NPK           - gains_Control,
            d_losses          = losses_NPK          - losses_Control) %>%
  group_by(site_code, block) %>%
  summarise(across(starts_with("d_"), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# 3.2 Diferencias NPK − Control en riqueza, equitatividad, rangos y especies,
#     promediadas por bloque
d_dif <- rac_diff %>%
  filter(year_trt > anio_min) %>%
  group_by(site_code, block) %>%
  summarise(across(c(richness_diff, evenness_diff, rank_diff, species_diff),
                   ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")


# 4. PREDICTORES A NIVEL DE SITIO (X) ------------------------------------------
predictores <- selected_sites %>%
  left_join(sindromes %>% select(site_code, category_2), by = "site_code") %>%
  left_join(metadat_clim %>%
              select(site_code, MAP_v2, MAT_v2, MAP_VAR_v2, 
                     MAT_RANGE_v2,site_richness,TEMP_VAR_v2,
                     AI) %>%
              distinct(),
            by = "site_code") %>%
  left_join(suelo %>%
              filter(trt == "Control") %>%
              select(site_code, pct_N, ppm_P, ppm_K, pH),
            by = "site_code") %>%
  left_join(comunidad   %>% select(-any_of("X")), by = "site_code") %>%
  left_join(clima_inter %>% select(-any_of("X")), by = "site_code") %>%
  # Transformaciones para reducir asimetría
  mutate(MAP_v2   = log(MAP_v2),
         pct_N    = log(pct_N),
         ppm_P    = log(ppm_P),
         ppm_K    = log(ppm_K),
         richness = log(richness),
         site_richness = log(site_richness))


# 5. UNIÓN, FILTRADO Y ESTANDARIZACIÓN -----------------------------------------
# Y y X se arman desde una única tabla para garantizar que las filas coincidan.
datos <- d_cambio %>%
  inner_join(d_dif,       by = c("site_code", "block")) %>%
  inner_join(predictores, by = "site_code") %>%
  filter(block <= bloques_max,
         !site_code %in% sitios_excl,
         category_2 != sindrome_excl
  ) %>%
  drop_na(all_of(c(vars_Y, vars_X, "category_2"))) %>%
  mutate(sample_id  = paste(site_code, block, sep = "_"),
         category_2 = droplevels(factor(category_2))) %>%
  arrange(site_code, block)

cat("Sitios:", n_distinct(datos$site_code), "| Filas (sitio × bloque):", nrow(datos), "\n")
print(table(datos$category_2) / bloques_max)  # n de sitios por síndrome

# Matriz Y estandarizada (media 0, varianza 1)
Y_sc <- datos %>%
  column_to_rownames("sample_id") %>%
  select(all_of(vars_Y)) %>%
  as.matrix() %>%
  scale()

# Matriz X: síndrome + predictores estandarizados
X <- datos %>%
  column_to_rownames("sample_id") %>%
  select(site_code, category_2, all_of(vars_X)) %>%
  mutate(across(all_of(vars_X), ~ as.numeric(scale(.x))))

# Controles
stopifnot(identical(rownames(Y_sc), rownames(X)))
stopifnot(length(unique(table(X$site_code))) == 1)  # mismo n de bloques por sitio


# 6. DISEÑO DE PERMUTACIONES ---------------------------------------------------
# Los predictores y el síndrome se miden a nivel de sitio: se permutan sitios
# completos, sin permutar bloques dentro de cada sitio.
ctrl <- how(plots  = Plots(strata = X$site_code, type = "free"),
            within = Within(type = "none"),
            nperm  = n_perm)


# Last RDA


bio   <- "SimpsonEvenness + Evar + richness + site_richness"
clima <- "MAP_v2 + MAT_v2 + MAP_VAR_v2 + TEMP_VAR_v2 + avg_cvT"
suelo <- "pct_N + ppm_P + ppm_K + pH"

test_bloque <- function(bloque, condicion) {
  set.seed(123)
  anova(rda(as.formula(paste("Y_sc ~", bloque, "+ Condition(", condicion, ")")), data = X),
        permutations = ctrl)
}

test_bloque(bio,   paste("category_2 +", clima, "+", suelo))
test_bloque(clima, paste("category_2 +", bio,   "+", suelo))
test_bloque(suelo, paste("category_2 +", bio,   "+", clima))
# Fracción única del síndrome (falta)
test_bloque("category_2", paste(bio, "+", clima, "+", suelo))

# Modelo de la figura: ajuste global, R² y colinealidad (falta)
rda_bio <- rda(Y_sc ~ category_2 + SimpsonEvenness + Evar + richness + site_richness, data = X)
RsquareAdj(rda_bio)
set.seed(123)
anova(rda_bio, permutations = ctrl)
vif.cca(rda_bio)

varpart(Y_sc, ~ category_2, as.formula(paste("~", bio)),
        as.formula(paste("~", clima)), as.formula(paste("~", suelo)), data = X)


set.seed(123)
anova(rda(Y_sc ~ category_2 + SimpsonEvenness + Evar + richness + site_richness, data = X),
      by = "margin", permutations = ctrl)

anova(rda(Y_sc ~ category_2 + SimpsonEvenness + Evar + richness + site_richness, data = X),
      by = "axis", permutations = ctrl)


bio_terms <- "SimpsonEvenness + Evar + richness + site_richness"
pares <- combn(levels(X$category_2), 2, simplify = FALSE)

comparaciones <- do.call(rbind, lapply(pares, function(par) {
  idx   <- X$category_2 %in% par
  X_par <- droplevels(X[idx, ])
  Y_par <- Y_sc[idx, ]
  ctrl_p <- how(plots  = Plots(strata = X_par$site_code, type = "free"),
                within = Within(type = "none"), nperm = 9999)
  set.seed(123)
  prueba <- anova(rda(as.formula(paste("Y_par ~ category_2 + Condition(", bio_terms, ")")),
                      data = X_par), permutations = ctrl_p)
  data.frame(comparacion = paste(par, collapse = " vs. "),
             n_sitios = n_distinct(X_par$site_code),
             F = prueba$F[1], p = prueba$`Pr(>F)`[1])
}))
comparaciones$p_holm <- p.adjust(comparaciones$p, method = "holm")
comparaciones

# Homogeneidad de la dispersión entre síndromes (se había borrado)
Y_sitio    <- aggregate(Y_sc, by = list(site_code = X$site_code), FUN = mean)
sind_sitio <- X %>% distinct(site_code, category_2) %>%
  slice(match(Y_sitio$site_code, site_code)) %>% pull(category_2)
disp_sind  <- betadisper(dist(Y_sitio[, -1]), group = sind_sitio)
tapply(disp_sind$distances, sind_sitio, mean)
set.seed(123)
test_disp <- permutest(disp_sind, permutations = n_perm, pairwise = TRUE)
test_disp
p.adjust(test_disp$pairwise$permuted, method = "holm")



# %%%%%----
# new plot -----
# %%%%%----

# ==============================================================================
# FIGURE 6 — RDA OF COMMUNITY PROCESSES (RDA1 x RDA2)
#
#   a) Scores of the response variables (processes) on RDA1
#   b) Scores of the response variables (processes) on RDA2
#      (bars coloured by r², on a colour scale shared by a and b)
#   c) Triplot RDA1 x RDA2: sites (mean of blocks) coloured by response
#      syndrome, 1-SD ellipses, syndrome centroids and arrows of the biotic
#      predictors (significance from marginal permutation tests)
#
# Model: syndrome + biotic block (Simpson evenness, E_var, local richness,
#        regional richness), defined a priori (see block tests).
# Requires objects from rda_analisis.R: Y_sc, X, datos, ctrl
# Figure text is defined in LABELS and exported in every language in FIG_LANGS.
# ==============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(vegan)
library(permute)


# 0. SETTINGS ------------------------------------------------------------------

FIG_LANGS   <- c("es", "en")
OUT_DIR     <- "./results"
N_PERM      <- 9999             # permutations for axis and marginal tests
SCALING     <- 3                # symmetric scaling, used in all panels
SIGN_AXES   <- c(-1, 1)         # sign of RDA1 and RDA2 (arbitrary; applied to all panels)
ARROW_FILL  <- 0.85             # predictor arrows span this fraction of the site cloud
ELLIPSE_LVL <- 1 - exp(-1/2)    # 1-SD ellipse of a bivariate normal
SHOW_R2     <- TRUE             # print r² values next to the bars

PREDICTORS <- c("SimpsonEvenness", "Evar", "richness", "site_richness")
VARS_DIFF  <- c("species_diff", "rank_diff", "evenness_diff", "richness_diff")
VARS_CHG   <- c("d_rank_change", "d_gains", "d_losses", "d_evenness_change",
                "d_richness_change")
SYN_LEVELS <- c("difmean_convergence", "difmean_eqvar", "difmean_divergence")

SYN_COLORS <- c(difmean_convergence = "#30123B",
                difmean_eqvar     = "#20EAAB",
                difmean_divergence = "#FEA632")


# 1. LABELS (ES / EN) ----------------------------------------------------------

LABELS <- list(
  es = list(
    syndrome = c(difmean_convergence = "Convergencia temporal",
                 difmean_eqvar     = "Cambio paralelo",
                 difmean_divergence = "Divergencia temporal"),
    process = c(richness_diff = "Riqueza", evenness_diff = "Equitatividad",
                rank_diff = "Rangos", species_diff = "Especies",
                d_richness_change = "Riqueza", d_evenness_change = "Equitatividad",
                d_rank_change = "Reordenamiento", d_gains = "Ganancias", d_losses = "Pérdidas"),
    # plotmath expressions: quoted text, E[var] as subscript
    predictor = c(SimpsonEvenness = "'Equitatividad de Simpson'", Evar = "E[var]",
                  richness = "'Riqueza local'", site_richness = "'Riqueza regional'"),
    type  = c(diff = "Diferencias\nNPK − Control", chg = "Cambios\ninteranuales"),
    score = "Score en RDA%d (%s%%)",
    r2    = "r²",
    decimal = ","
  ),
  en = list(
    syndrome = c(difmean_convergence = "Temporal convergence",
                 difmean_eqvar     = "Parallel change",
                 difmean_divergence = "Temporal divergence"),
    process = c(richness_diff = "Richness", evenness_diff = "Evenness",
                rank_diff = "Ranks", species_diff = "Species",
                d_richness_change = "Richness", d_evenness_change = "Evenness",
                d_rank_change = "Reordering", d_gains = "Gains", d_losses = "Losses"),
    predictor = c(SimpsonEvenness = "'Simpson evenness'", Evar = "E[var]",
                  richness = "'Local richness'", site_richness = "'Regional richness'"),
    type  = c(diff = "NPK − Control\ndifferences", chg = "Interannual\nchanges"),
    score = "RDA%d score (%s%%)",
    r2    = "r²",
    decimal = "."
  )
)

stars <- function(p) as.character(cut(p, c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                                      c("***", "**", "*", "·", ""), right = FALSE))

# Axis labels with the decimal mark of each language (base R, no extra packages)
dec_labels <- function(mark) function(x) format(x, decimal.mark = mark, trim = TRUE,
                                                drop0trailing = TRUE)
fmt_num <- function(x, d, mark) formatC(x, format = "f", digits = d, decimal.mark = mark)


# 2. MODEL AND TESTS -----------------------------------------------------------

form    <- as.formula(paste("Y_sc ~ category_2 +", paste(PREDICTORS, collapse = " + ")))
rda_bio <- rda(form, data = X)

ctrl_fig <- ctrl
setNperm(ctrl_fig) <- N_PERM

set.seed(123)
test_axes   <- anova(rda_bio, by = "axis",   permutations = ctrl_fig)
test_margin <- anova(rda_bio, by = "margin", permutations = ctrl_fig)
print(test_axes); print(test_margin)

pct_axis <- 100 * rda_bio$CCA$eig[1:2] / rda_bio$tot.chi
p_axis   <- test_axes[["Pr(>F)"]][1:2]


# 3. SCORES --------------------------------------------------------------------

flip <- function(df) {
  df[, 1] <- SIGN_AXES[1] * df[, 1]
  df[, 2] <- SIGN_AXES[2] * df[, 2]
  df
}

# Processes: scores on RDA1 and RDA2 + r² (squared correlation with LC scores)
proc_sc <- flip(as.data.frame(scores(rda_bio, display = "species", choices = 1:2,
                                     scaling = SCALING)))
proc_r2 <- as.data.frame(cor(Y_sc, scores(rda_bio, display = "lc", choices = 1:2))^2)

proc <- bind_rows(
  tibble(var = rownames(proc_sc), axis = 1, score = proc_sc[, 1], r2 = proc_r2[rownames(proc_sc), 1]),
  tibble(var = rownames(proc_sc), axis = 2, score = proc_sc[, 2], r2 = proc_r2[rownames(proc_sc), 2])
) %>%
  mutate(type = ifelse(var %in% VARS_DIFF, "diff", "chg"),
         var  = factor(var, levels = rev(c(VARS_DIFF, VARS_CHG))))

r2_max <- max(proc$r2)

# Sites: weighted-average scores, averaged over the blocks of each site
site_sc <- flip(as.data.frame(scores(rda_bio, display = "sites", choices = 1:2,
                                     scaling = SCALING)))
stopifnot(identical(rownames(site_sc), datos$sample_id))

sites <- site_sc %>%
  mutate(site_code = datos$site_code, category = as.character(datos$category_2)) %>%
  group_by(site_code, category) %>%
  summarise(RDA1 = mean(RDA1), RDA2 = mean(RDA2), .groups = "drop") %>%
  mutate(category = factor(category, levels = SYN_LEVELS))

centroids <- sites %>%
  group_by(category) %>%
  summarise(RDA1 = mean(RDA1), RDA2 = mean(RDA2), .groups = "drop")

# Biotic predictors (syndrome indicator variables are shown as centroids instead)
pred <- flip(as.data.frame(scores(rda_bio, display = "bp", choices = 1:2,
                                  scaling = SCALING))) %>%
  rownames_to_column("var") %>%
  filter(var %in% PREDICTORS) %>%
  mutate(p = test_margin[var, "Pr(>F)"])

mult <- ARROW_FILL * max(abs(c(sites$RDA1, sites$RDA2))) / max(sqrt(pred$RDA1^2 + pred$RDA2^2))
pred <- pred %>% mutate(x = RDA1 * mult, y = RDA2 * mult)

print(proc, n = Inf); print(pred)


# 4. PLOTS ---------------------------------------------------------------------

base_theme <- theme_classic(base_size = 11) +
  theme(plot.tag         = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text.y     = element_text(angle = 0, hjust = 0, face = "italic"),
        legend.position  = "bottom")

axis_title <- function(L, k) {
  paste(sprintf(L$score, k, fmt_num(pct_axis[k], 1, L$decimal)), stars(p_axis[k]))
}

# (a, b) process scores on one axis
bar_panel <- function(L, k, tag) {
  d <- filter(proc, axis == k) %>%
    mutate(type = factor(L$type[type], levels = L$type))
  lim <- max(abs(d$score)) * 1.35
  
  p <- ggplot(d, aes(score, var, fill = r2)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    facet_grid(type ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = L$process) +
    scale_fill_gradient(low = "grey88", high = "grey15", limits = c(0, r2_max),
                        name = L$r2,
                        labels = dec_labels(L$decimal),
                        guide = guide_colorbar(barwidth = unit(3, "cm"), title.vjust = 0.8)) +
    scale_x_continuous(labels = dec_labels(L$decimal)) +
    coord_cartesian(xlim = c(-lim, lim)) +
    labs(x = axis_title(L, k), y = NULL, tag = tag) +
    base_theme
  
  if (SHOW_R2) {
    p <- p + geom_text(aes(label = fmt_num(r2, 2, L$decimal),
                           hjust = ifelse(score >= 0, -0.15, 1.15)),
                       size = 2.8, colour = "grey20")
  }
  p
}

# (c) triplot with biotic predictors
triplot_panel <- function(L) {
  lims <- range(c(sites$RDA1, sites$RDA2, pred$x, pred$y)) * 1.15
  
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    stat_ellipse(data = sites, aes(RDA1, RDA2, colour = category, fill = category),
                 geom = "polygon", type = "norm", level = ELLIPSE_LVL,
                 alpha = 0.12, linewidth = 0.4) +
    geom_point(data = sites, aes(RDA1, RDA2, colour = category), size = 2.4,
               alpha = 0.6) +
    geom_point(data = centroids, aes(RDA1, RDA2, fill = category),
               shape = 23, size = 2.6, colour = "black", show.legend = FALSE) +
    geom_segment(data = pred, aes(0, 0, xend = x, yend = y),
                 arrow = arrow(length = unit(0.18, "cm")),
                 colour = "firebrick", linewidth = 0.7) +
    geom_text_repel(data = pred,
                    aes(x, y, label = paste0(L$predictor[var], "~'",
                                             stars(p), "'")),
                    parse = TRUE, size = 3.2, colour = "firebrick4",
                    box.padding = 0.5, max.overlaps = Inf,
                    min.segment.length = 0, seed = 1) +
    scale_colour_manual(values = SYN_COLORS, breaks = SYN_LEVELS, labels = L$syndrome,
                        name = NULL) +
    scale_fill_manual(values = SYN_COLORS, breaks = SYN_LEVELS, labels = L$syndrome,
                      name = NULL) +
    scale_x_continuous(labels = dec_labels(L$decimal)) +
    scale_y_continuous(labels = dec_labels(L$decimal)) +
    coord_equal(xlim = lims, ylim = lims) +
    labs(x = axis_title(L, 1), y = axis_title(L, 2), tag = "c") +
    base_theme
}

make_figure <- function(lang) {
  L <- LABELS[[lang]]
  ((bar_panel(L, 1, "a") / bar_panel(L, 2, "b")) | triplot_panel(L)) +
    plot_layout(widths = c(1, 1.25), guides = "collect") &
    theme(legend.position = "bottom", legend.box = "vertical")
}


# 5. EXPORT --------------------------------------------------------------------

for (lang in FIG_LANGS) {
  fig <- make_figure(lang)
  #print(fig)
  ggsave(file.path(OUT_DIR, paste0("fig6_rda_", lang, ".png")),
         fig, width = 8, height = 6, dpi = 300)
}
