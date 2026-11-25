set.seed(123)
# data ----

library(tidyverse)

glmulti_difvariance_result<-read.csv("./analysis/results/glmulti_difvariance_result.csv") # no need to run the glmulti model
glmulti_centroiddist_result<-read.csv("./analysis/results/glmulti_centrdoiddist_result.csv") # no need to run the glmulti model

metadat_<-read.csv("./nutnet_data/comb-by-plot-clim-soil-2025-12-09.csv", sep = ",",
                   na.strings = c("NULL", "", "N/A")) 


sites <- read.csv("./my_data/selected_sites.csv")
sites<-sites%>%dplyr::select(site_code)

#variance_trt_bloc_dat <- read.csv("./analysis/results/variance_by_trt_block.csv",
#sep=",")[,-1]

dif_variance_dat <- read.csv("./results/centroid-dist_variance-diff.csv")[,-1]

#dif_variance_df <- sites %>%
 # left_join(dif_variance_dat,by="site_code")

amb_metadat <- sites %>%
  left_join(metadat_  %>%
              dplyr::select(burned,site_code,elevation,latitude,
                            RAIN_PET,MAT_v2,MAT_RANGE_v2,ISO_v2,TEMP_VAR_v2,
                            MIN_TEMP_v2,MAX_TEMP_v2,site_richness,MAT_RANGE_v2,
                            MAP_v2,MAP_VAR_v2,AI,PET,N_Dep,elevation,ISO_v2,
                            TEMP_WET_Q_v2) %>% distinct(),
            by="site_code")

soil_dat <-read.csv("./my_data/soil-data_summarised.csv",sep=",") %>%
  filter(trt=="Control") %>%
  dplyr::select(!X) %>% dplyr::select(!trt)

rich_dat <- read.csv("./my_data/richness_evennes.csv")

sync_dat <- read.csv("./my_data/synchrony_data.csv")

sp_stab_dat <- read.csv("./my_data/control_max_stability.csv")


df <- dif_variance_df%>%select(site_code,centroid_distance,diff_variance,category_2) %>%
  left_join(amb_metadat %>%
              select(site_code,#max_year,
                     burned,
                     MAT_v2,TEMP_VAR_v2,MAP_v2,MAP_VAR_v2,ISO_v2,elevation,
                     TEMP_WET_Q_v2,site_richness,MAT_RANGE_v2,
                     AI),by="site_code")%>%
  left_join(soil_dat%>%
              select(site_code,pct_C,pct_N,ppm_P,ppm_K,pH),
            by="site_code") %>%
  left_join(rich_dat,by="site_code")%>%
  left_join(sp_stab_dat%>%
              select(site_code,Control_stability),by="site_code")%>%
  left_join(sync_dat,by="site_code")

# libraries ----
library(viridis)
library(ggpubr)
library(glmulti)

# GLMULTI ----

# check correlations first!
  library(GGally)
  
  ggcorr(df)



predictors2 <- c(
  
  "TEMP_VAR_v2",
  #"TEMP_WET_Q_v2", #NEW!
  "log(MAP_v2)",
  "log(richness)",
  "Evar",
  "log(site_richness)",
  "SimpsonEvenness", # no need for log
  #"log(Control_stability)", 
  "log(MAP_VAR_v2)", 
  #"AI", # no need for log
  "MAT_v2",
  #"log(Control_sync)",
  "pH",
  "log(pct_N)",
  "log(ppm_P)",
  "log(ppm_K)"
)

df_glm <- df %>%
  filter(category_2 != "eqmean_eqvar") #exclude no-response sites (n=2)

maxsize_glm <- 3
n_models    <- sum(choose(length(predictors2), 0:maxsize_glm))   # all candidate models

set.seed(123)
glmulti_difvariance2 <- glmulti(
  y           = "diff_variance",
  xr          = predictors2,
  data        = df_glm,
  fitfunction = "glm",
  family      = gaussian(),
  level       = 1,
  maxsize     = maxsize_glm,
  crit        = "aicc",
  method      = "h",
  confsetsize = n_models,
  plotty      = FALSE,
  report      = FALSE
)

set.seed(123)
glmulti_centroiddist2 <- glmulti(
  y           = "centroid_distance",
  xr          = predictors2,
  data        = df_glm,
  fitfunction = "glm",
  family      = gaussian(),
  level       = 1,
  maxsize     = maxsize_glm,
  crit        = "aicc",
  method      = "h",
  confsetsize = n_models,
  plotty      = FALSE,
  report      = FALSE
)

# Check: the intercept should appear in all candidate models
n_models

coef(glmulti_difvariance2)
coef(glmulti_centroiddist2)

plot(glmulti_difvariance2,  type = "s")
plot(glmulti_centroiddist2, type = "s")

## extract models output ----

### diff variance ----

mmi2 <- as.data.frame(coef(glmulti_difvariance2, varweighting="Johnson"))
mmi2 <- data.frame(Estimate=mmi2$Est, SE=sqrt(mmi2$Uncond), Importance=mmi2$Importance, row.names=row.names(mmi2))
mmi2$z <- mmi2$Estimate / mmi2$SE
mmi2$p <- 2*pnorm(abs(mmi2$z), lower.tail=FALSE)
names(mmi2) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi2$ci.lb <- mmi2[[1]] - qnorm(.975) * mmi2[[2]]
mmi2$ci.ub <- mmi2[[1]] + qnorm(.975) * mmi2[[2]]
mmi2 <- mmi2[order(mmi2$Importance, decreasing=TRUE), c(1,2,4:7,3)]
glmulti_difvariance_result2 <-data.frame(round(mmi2, 4))



### centroid distance ----

mmi.3 <- as.data.frame(coef(glmulti_centroiddist2, varweighting="Johnson"))
mmi.3 <- data.frame(Estimate=mmi.3$Est, SE=sqrt(mmi.3$Uncond), Importance=mmi.3$Importance, row.names=row.names(mmi.3))
mmi.3$z <- mmi.3$Estimate / mmi.3$SE
mmi.3$p <- 2*pnorm(abs(mmi.3$z), lower.tail=FALSE)
names(mmi.3) <- c("Estimate", "Std. Error", "Importance", "z value", "Pr(>|z|)")
mmi.3$ci.lb <- mmi.3[[1]] - qnorm(.975) * mmi.3[[2]]
mmi.3$ci.ub <- mmi.3[[1]] + qnorm(.975) * mmi.3[[2]]
mmi.3 <- mmi.3[order(mmi.3$Importance, decreasing=TRUE), c(1,2,4:7,3)]

glmulti_centroiddist_result2 <-data.frame(round(mmi.3, 4))


 write.csv(glmulti_difvariance_result2,"./results/glmulti_difvariance_results.csv")
 write.csv(glmulti_centroiddist_result2,"./results/glmulti_centrdoiddist_results.csv")

# plots ----

glmulti_centd_results_df <- read.csv("./results/glmulti_centrdoiddist_results.csv")


glmulti_vardif_results_df <- read.csv("./results/glmulti_difvariance_results.csv")


ggplot(glmulti_centd_results_df[-1,],
       aes(x = reorder(X, Importance), y = Importance,
           fill = Estimate)) +
  geom_col()+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1))+
  geom_hline(yintercept = 0.8,linetype="dashed")+
  scale_fill_viridis_b(option = "turbo",
                       begin = 0.1,
                       end = 0.45,
                       direction = -1)+
  labs(fill="Slope",
       title = "Centroid distance")+
  scale_x_discrete(labels =c("(log) % TN",
                             "MAT", "Simpson Evenness",
                             "(log) ppm K",
                             "PH","(log) Sp. Synchrony",
                             "(log) MAP",
                             "AI","Temp. var.",
                             "(log) MAP var.",
                             "(log) ppm P",
                             "(log) sp. Richness"))+
  coord_flip()+
  theme_pubr()+theme(legend.position = "right",
                     axis.title.y = element_blank())


ggplot(glmulti_vardif_results_df[-1,],
       aes(x = reorder(X, Importance), y = Importance,
           fill = Estimate)) +
  geom_col()+
  geom_hline(yintercept = 0.8,linetype="dashed")+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1))+
  scale_fill_viridis_b(option = "turbo",
                       begin = 0.3,
                       end = 0.8,
                       direction = -1)+#,
#                       limits = c(-0.5, 0.031), oob = scales::squish)+
  scale_x_discrete(labels =c("(log) Sp. Synchrony",
                             "(log) ppm K",
                             "(log) MAP var.",
                             "(log) ppm P",
                             "(log) % TN",
                             "(log) sp. Richness",
                             "PH",
                             "(log) MAP",
                             "AI",
                             "MAT",
                             "Simpson Evenness",
                             "Temp. var."
                             
                             ))+
  labs(fill="Slope",
       title = "Dispersion")+
  coord_flip()+
  theme_pubr()+theme(legend.position = "right",
                     axis.title.y = element_blank())

