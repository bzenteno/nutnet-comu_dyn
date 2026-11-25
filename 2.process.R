# ==============================================================================
# COMMUNITY PROCESSES UNDERLYING THE RESPONSE SYNDROMES
#
#   Figure 4: NPK − Control differences in community structure
#             (codyn::RAC_difference; one value per block x year)
#   Figure 5: NPK effects on interannual community change
#             (codyn::RAC_change; one value per plot x pair of consecutive years)
#
# Both pipelines exclude sites without response (syndrome "a_eqmean_eqvar")
# and use only years after the initial transient phase (year_trt > MIN_YEAR).
#
# Figure text (axis titles, facet labels, legend, decimal mark) is defined in
# the LABELS dictionary (section 2) and exported in every language listed in
# FIG_LANGS.
# ==============================================================================


# 0. SETTINGS ------------------------------------------------------------------

FIG_LANGS          <- c("en", "es")   # languages of the exported figures
RECOMPUTE_RAC_DIFF <- FALSE           # TRUE: recompute RAC_difference from raw cover data
CHECK_RESIDUALS    <- TRUE           # TRUE: show DHARMa residual diagnostics
RUN_SENSITIVITY    <- FALSE           # TRUE: fit alternative model specifications
ADJUST_METHOD      <- "tukey"         # multiplicity adjustment, among-syndrome comparisons
MIN_YEAR           <- 3               # analyses use year_trt > MIN_YEAR
OUT_DIR            <- "./results"

FIG_WIDTH  <- 8      # inches
FIG_HEIGHT <- 4.5    # inches
FIG_DPI    <- 300


# 1. PACKAGES ------------------------------------------------------------------

library(multcomp)    # cld() generic; loaded before tidyverse so dplyr::select() is not masked
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(emmeans)     # cld() method for emmGrid objects also requires 'multcompView'
library(codyn)
library(ggpubr)
library(scales)


# 2. FIGURE LABELS (EN / ES) ---------------------------------------------------
# Edit here to change any text shown in the figures.

LABELS <- list(
  en = list(
    syndrome = c(
      b_difmean_convergence = "Temporal convergence\n(\u2193 dispersion)",
      c_difmean_eqvar     = "Parallel change\n(= dispersion)",
      d_difmean_divergence = "Temporal divergence\n(\u2191 dispersion)"
    ),
    fig4_facets = c(
      a_richness_diff = "Richness",
      b_evenness_diff = "Evenness",
      c_rank_diff     = "Rank\nabundance",
      d_species_diff  = "Species\ncomposition"
    ),
    fig4_y = "NPK \u2212 Control difference",
    fig5_facets = c(
      a_richness_change = "Richness\nchange",
      b_evenness_change = "Evenness\nchange",
      c_rank_change     = "Rank\nchange",
      d_gains           = "Species\ngains",
      e_losses          = "Species\nlosses"
    ),
    fig5_y       = "NPK effect (SD units)",
    decimal_mark = "."
  ),
  es = list(
    syndrome = c(
      b_difmean_convergence = "Convergencia temporal\n(\u2193 dispersi\u00f3n)",
      c_difmean_eqvar     = "Cambio paralelo\n(= dispersi\u00f3n)",
      d_difmean_divergence = "Divergencia temporal\n(\u2191 dispersi\u00f3n)"
    ),
    fig4_facets = c(
      a_richness_diff = "Riqueza",
      b_evenness_diff = "Equitatividad",
      c_rank_diff     = "Rangos de\nabundancia",
      d_species_diff  = "Composici\u00f3n\nde especies"
    ),
    fig4_y = "Diferencia NPK \u2212 Control",
    fig5_facets = c(
      a_richness_change = "Cambio de\nriqueza",
      b_evenness_change = "Cambio de\nequitatividad",
      c_rank_change     = "Reordenamiento",
      d_gains           = "Ganancias de\nespecies",
      e_losses          = "P\u00e9rdidas de\nespecies"
    ),
    fig5_y       = "Efecto de NPK (unidades de DE)",
    decimal_mark = ","
  )
)

SYNDROME_COLORS <- c(
  b_difmean_convergence = "#30123B",
  c_difmean_eqvar     = "#20EAAB",
  d_difmean_divergence = "#FEA632"
)


# 3. HELPER FUNCTIONS ----------------------------------------------------------

# Syndrome codes with an alphabetical prefix that sets the plotting order
recode_syndrome <- function(x) {
  dplyr::recode(x,
                eqmean_eqvar      = "a_eqmean_eqvar",
                difmean_convergence = "b_difmean_convergence",
                difmean_eqvar     = "c_difmean_eqvar",
                difmean_divergence = "d_difmean_divergence")
}

# Harmonise confidence-limit column names returned by emmeans
get_ci <- function(df) {
  lower <- intersect(names(df), c("asymp.LCL", "lower.CL"))[1]
  upper <- intersect(names(df), c("asymp.UCL", "upper.CL"))[1]
  df$ci_lower <- df[[lower]]
  df$ci_upper <- df[[upper]]
  df
}

p_stars <- function(p) {
  dplyr::case_when(p < 0.001 ~ "***",
                   p < 0.01  ~ "**",
                   p < 0.05  ~ "*",
                   TRUE      ~ "")
}

check_residuals <- function(model) {
  if (CHECK_RESIDUALS) plot(simulateResiduals(model))
  invisible(model)
}


# 4. SHARED DATA ---------------------------------------------------------------

selected_sites <- read.csv("my_data/selected_sites.csv") %>%
  dplyr::select(site_code)

plot_data <- read.csv("./nutnet_data/comb-by-plot-2025-12-09.csv", sep = ";")

plot_meta <- plot_data %>%
  dplyr::distinct(site_code, block, plot)

# Sites with at least 10 years of treatment (excluding sites flagged for removal)
long_sites <- selected_sites
  #plot_data %>%
  #group_by(site_code) %>%
  #summarise(max_year = max(year_trt), .groups = "drop") %>%
  #filter(max_year > 9,
  #       !site_code %in% c("doane.us", "gilb.za", "cereep.fr")) %>%
  #dplyr::select(site_code)

# Response syndrome of each site (one row per site)
syndromes <- read.csv("./results/centroid-dist_variance-diff.csv")[, -1] %>%
  dplyr::distinct(site_code, category = category_2) %>%
  mutate(category = recode_syndrome(category))

stopifnot(!any(duplicated(syndromes$site_code)))


# ==============================================================================
# PART A — FIGURE 4: NPK − CONTROL DIFFERENCES (RAC_difference)
# ==============================================================================

# A.1 RAC_difference (optional recomputation) ----------------------------------
# codyn::RAC_difference requires exactly one NPK and one Control plot per
# block x year. One plot per block x year x treatment is kept (lowest plot ID),
# and only block x year combinations containing both treatments are used.

if (RECOMPUTE_RAC_DIFF) {
  
  cover_data <- read.csv("./nutnet_data/full-cover-2025-12-09.csv", sep = ";")
  
  cover_one_plot <- selected_sites %>%
    left_join(cover_data, by = "site_code") %>%
    filter(trt %in% c("Control", "NPK"), live == 1) %>%
    group_by(site_code, year_trt, block, trt) %>%
    filter(plot == min(plot, na.rm = TRUE)) %>%
    ungroup()
  
  paired_blocks <- cover_one_plot %>%
    dplyr::distinct(site_code, block, year_trt, trt) %>%
    count(site_code, block, year_trt, trt) %>%
    pivot_wider(names_from = trt, values_from = n, values_fill = 0) %>%
    filter(Control == 1, NPK == 1) %>%
    dplyr::select(site_code, block, year_trt)
  
  rac_diff <- paired_blocks %>%
    inner_join(cover_one_plot, by = c("site_code", "block", "year_trt")) %>%
    group_by(site_code, year_trt) %>%
    group_modify(~ RAC_difference(df                  = .x,
                                  species.var         = "Taxon",
                                  abundance.var       = "max_cover",
                                  replicate.var       = "plot",
                                  treatment.var       = "trt",
                                  block.var           = "block",
                                  reference.treatment = "Control")) %>%
    ungroup()
  
  write.csv(rac_diff, "./my_data/RAC_diff_results.csv", row.names = FALSE)
}

rac_diff <- read.csv("./my_data/RAC_diff_results.csv")


# A.2 Analysis data ------------------------------------------------------------

df_diff <- long_sites %>%
  inner_join(rac_diff, by = "site_code") %>%
  filter(year_trt > MIN_YEAR) %>%
  inner_join(syndromes, by = "site_code") %>%
  filter(category != "a_eqmean_eqvar") %>%
  mutate(
    # Unique block ID: block labels are repeated across sites
    block_uid = factor(paste(site_code, block, sep = "_")),
    year_f    = factor(year_trt),
    # Rank difference rescaled to 0–1 [justify the 0.5 denominator in Methods]
    rank_diff_rescaled = rank_diff / 0.5
  )


# A.3 Models -------------------------------------------------------------------
# Response ~ syndrome (no intercept), site random effect, AR(1) among years
# within block. Family and dispersion structure chosen from residual diagnostics.

mod_richness_diff <- glmmTMB(
  richness_diff ~ 0 + category + ar1(year_f + 0 | block_uid) + (1 | site_code),
  dispformula = ~ category,
  data = df_diff, na.action = na.omit
)
check_residuals(mod_richness_diff)

mod_evenness_diff <- glmmTMB(
  evenness_diff ~ 0 + category + ar1(year_f + 0 | block_uid) + (1 | site_code),
  dispformula = ~ site_code,
  family = t_family(),
  data = df_diff, na.action = na.omit
)
check_residuals(mod_evenness_diff)

mod_rank_diff <- glmmTMB(
  rank_diff_rescaled ~ 0 + category + ar1(year_f + 0 | block_uid) + (1 | site_code),
  dispformula = ~ site_code,
  family = t_family(),
  data = df_diff, na.action = na.omit
)
check_residuals(mod_rank_diff)

mod_species_diff <- glmmTMB(
  species_diff ~ 0 + category + ar1(year_f + 0 | block_uid) + (1 | site_code),
  dispformula = ~ category,
  data = df_diff, na.action = na.omit
)
check_residuals(mod_species_diff)


# A.4 Estimated means, tests against zero and among-syndrome comparisons -------

summarise_diff <- function(model, process) {
  emm <- emmeans(model, ~ category)
  
  estimates <- as.data.frame(summary(emm, infer = c(TRUE, TRUE))) %>%
    get_ci()
  
  groups <- as.data.frame(cld(emm, adjust = ADJUST_METHOD, Letters = letters)) %>%
    transmute(category, group = trimws(.group))
  
  estimates %>%
    left_join(groups, by = "category") %>%
    transmute(process = .env$process, category,
              estimate = emmean, SE, ci_lower, ci_upper,
              p_value = p.value, group)
}

fig4_data <- bind_rows(
  summarise_diff(mod_richness_diff, "a_richness_diff"),
  summarise_diff(mod_evenness_diff, "b_evenness_diff"),
  summarise_diff(mod_rank_diff,     "c_rank_diff"),
  summarise_diff(mod_species_diff,  "d_species_diff")
)

fig4_data
write.csv(fig4_data, file.path(OUT_DIR, "fig4_rac_difference_estimates.csv"),
          row.names = FALSE)


# ==============================================================================
# PART B — FIGURE 5: NPK EFFECTS ON INTERANNUAL CHANGE (RAC_change)
# ==============================================================================

# B.1 Analysis data ------------------------------------------------------------
# One row per plot x pair of consecutive years (year_trt -> year_trt2).

rac_change <- read.csv("./my_data/RAC_change_results.csv") %>%
  dplyr::select(-any_of("X")) %>%
  filter(trt %in% c("Control", "NPK"),
         year_trt > MIN_YEAR)

df_change <- long_sites %>%
  inner_join(rac_change, by = "site_code") %>%
  inner_join(syndromes,  by = "site_code") %>%
  filter(category != "a_eqmean_eqvar") %>%
  left_join(plot_meta, by = c("site_code", "plot")) %>%
  mutate(
    # Unique plot ID: plot labels are repeated across sites
    plot_uid = factor(paste(site_code, block, plot, sep = "_")),
    year_f   = factor(year_trt2)
  )


# B.2 Models -------------------------------------------------------------------
# Response ~ treatment x syndrome. Richness, evenness and rank change are
# standardised (scale) so that NPK effects are expressed in SD units; gains and
# losses are modelled with a Tweedie distribution (log link) and their effects
# are standardised afterwards (section B.3). Random-effect, AR(1) and dispersion
# structures were simplified when needed to reach convergence.

# AR(1) term omitted: the model did not converge with it
mod_richness_change <- glmmTMB(
  scale(richness_change) ~ trt * category + (1 | site_code/block/plot),
  dispformula = ~ trt + category,
  data = df_change, na.action = na.omit
)
check_residuals(mod_richness_change)

mod_evenness_change <- glmmTMB(
  scale(evenness_change) ~ trt * category +
    ar1(year_f + 0 | plot_uid) + (1 | site_code/block/plot),
  dispformula = ~ site_code + trt,
  family = t_family(),
  data = df_change, na.action = na.omit
)
check_residuals(mod_evenness_change)

# Block and plot random effects removed (overparameterisation)
mod_rank_change <- glmmTMB(
  scale(rank_change) ~ trt * category +
    ar1(year_f + 0 | plot_uid) + (1 | site_code),
  dispformula = ~ site_code + trt,
  data = df_change, na.action = na.omit
)
check_residuals(mod_rank_change)

mod_gains <- glmmTMB(
  gains ~ trt * category +
    ar1(year_f + 0 | plot_uid) + (1 | site_code/block/plot),
  dispformula = ~ site_code + trt,
  family = tweedie(),
  data = df_change, na.action = na.omit
)
check_residuals(mod_gains)

# Block and plot random effects removed (convergence)
mod_losses <- glmmTMB(
  losses ~ trt * category +
    ar1(year_f + 0 | plot_uid) + (1 | site_code),
  dispformula = ~ site_code * trt,
  family = tweedie(),
  data = df_change, na.action = na.omit
)
check_residuals(mod_losses)


# B.3 NPK − Control contrasts within each syndrome -----------------------------
# Tweedie models: contrasts are computed on the response scale (regrid) and
# divided by the raw SD of the variable, so all effects are in SD units.

summarise_change <- function(model, process, tweedie = FALSE, raw_sd = 1) {
  emm <- emmeans(model, ~ trt | category)
  if (tweedie) emm <- regrid(emm, transform = "response")
  
  npk_effect <- contrast(emm, method = "revpairwise")   # NPK − Control
  
  estimates <- as.data.frame(summary(npk_effect, infer = c(TRUE, TRUE))) %>%
    get_ci()
  
  # Among-syndrome comparison of the NPK effect
  groups <- as.data.frame(cld(npk_effect, by = NULL, adjust = ADJUST_METHOD,
                              Letters = letters)) %>%
    transmute(category, group = trimws(.group))
  
  estimates %>%
    left_join(groups, by = "category") %>%
    transmute(process  = .env$process, category,
              estimate = estimate / raw_sd,
              SE       = SE / raw_sd,
              ci_lower = ci_lower / raw_sd,
              ci_upper = ci_upper / raw_sd,
              p_value  = p.value,
              group)
}

sd_gains  <- sd(df_change$gains,  na.rm = TRUE)
sd_losses <- sd(df_change$losses, na.rm = TRUE)

fig5_data <- bind_rows(
  summarise_change(mod_richness_change, "a_richness_change"),
  summarise_change(mod_evenness_change, "b_evenness_change"),
  summarise_change(mod_rank_change,     "c_rank_change"),
  summarise_change(mod_gains,  "d_gains",  tweedie = TRUE, raw_sd = sd_gains),
  summarise_change(mod_losses, "e_losses", tweedie = TRUE, raw_sd = sd_losses)
)

fig5_data
write.csv(fig5_data, file.path(OUT_DIR, "fig5_rac_change_contrasts.csv"),
          row.names = FALSE)


# B.4 Sensitivity checks (optional) --------------------------------------------
# Not used for the figures.

if (RUN_SENSITIVITY) {
  
  # Rank change as a proportion: beta family with zero inflation
  # (~2.7% exact zeros, not admitted by beta_family)
  mod_rank_change_beta <- glmmTMB(
    rank_change ~ trt * category +
      ar1(year_f + 0 | plot_uid) + (1 | site_code/block),
    dispformula = ~ site_code + trt,
    ziformula   = ~ 1,
    family      = beta_family(),
    data = df_change, na.action = na.omit
  )
  check_residuals(mod_rank_change_beta)
  summary(mod_rank_change_beta)
  
  # Losses with the full nested random-effect structure
  mod_losses_full_re <- glmmTMB(
    losses ~ trt * category +
      ar1(year_f + 0 | plot_uid) + (1 | site_code/block/plot),
    dispformula = ~ category + trt,
    family      = tweedie(),
    data = df_change, na.action = na.omit
  )
  check_residuals(mod_losses_full_re)
  summary(mod_losses_full_re)
}


# ==============================================================================
# PART C — FIGURES
# ==============================================================================

# Bars: estimate ± 95% CI per syndrome.
# Asterisks: test against zero. Letters: among-syndrome comparisons
# (hidden in panels where all syndromes share the same letter).

plot_effects <- function(df, lang, facet_key, y_key, free_y = FALSE) {
  
  L <- LABELS[[lang]]
  
  d <- df %>%
    mutate(process  = factor(process, levels = names(L[[facet_key]])),
           category = factor(category, levels = names(L$syndrome))) %>%
    group_by(process) %>%
    mutate(
      span      = max(ci_upper, 0) - min(ci_lower, 0),
      side      = if_else(estimate >= 0, 1, -1),
      edge      = if_else(estimate >= 0, ci_upper, ci_lower),
      stars     = p_stars(p_value),
      stars_y   = edge + side * 0.06 * span,
      letters_y = edge + side * 0.16 * span,
      group     = if (n_distinct(group) == 1) "" else group
    ) %>%
    ungroup()
  
  ggplot(d, aes(x = category, y = estimate, fill = category)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
    geom_col(width = 1, colour = "black") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                  width = 0.2, linewidth = 0.5) +
    geom_text(aes(y = stars_y,   label = stars), size = 5,   vjust = 0.5) +
    geom_text(aes(y = letters_y, label = group), size = 3.8, vjust = 0.5) +
    facet_wrap(~ process, nrow = 1,
               scales   = if (free_y) "free_y" else "fixed",
               labeller = as_labeller(L[[facet_key]])) +
    scale_fill_manual(values = SYNDROME_COLORS,
                      breaks = names(L$syndrome),
                      labels = L$syndrome) +
    scale_y_continuous(labels = label_number(decimal.mark = L$decimal_mark)) +
    labs(x = NULL, y = L[[y_key]]) +
    theme_pubr() +
    theme(legend.position = "bottom",
          legend.title    = element_blank(),
          axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank())
}

for (lang in FIG_LANGS) {
  
  fig4 <- plot_effects(fig4_data, lang, facet_key = "fig4_facets", y_key = "fig4_y")
  fig5 <- plot_effects(fig5_data, lang, facet_key = "fig5_facets", y_key = "fig5_y")
  
  print(fig4)
  print(fig5)
  
  ggsave(file.path(OUT_DIR, paste0("fig4_rac_difference_", lang, ".png")),
         fig4, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI)
  ggsave(file.path(OUT_DIR, paste0("fig5_rac_change_", lang, ".png")),
         fig5, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI)
}
