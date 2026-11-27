library(tidyverse)
library(ggpubr)
library(viridis)

dat_disp <- read.csv("analysis/results/variance_by_trt_block.csv")

dat_synd <- read.csv("analysis/results/centroid-dist_variance-diff.csv")

sites <- read.csv("my_data/selected_sites.csv")

dat_rich <- read.csv("my_data/richness_evennes.csv")

metadat <- read.csv("nutnet_data/comb-by-plot-clim-soil-2025-12-09.csv")

df <- sites %>%
  left_join(dat_synd%>%dplyr::select(site_code,category_impact),
            by="site_code") %>%
  left_join(dat_disp %>% dplyr::select(trt,block,variance,site_code),
            by="site_code") #%>%
  #left_join(dat_rich,
   #         by="site_code") %>%
  #left_join(metadat%>%
#              filter(trt%in%c("Control","NPK"))%>%
 #             distinct(),by="site_code")


ggplot(df,aes(category_impact,variance,fill = trt))+geom_boxplot()

library(glmmTMB)
library(DHARMa)


m_1 <- glmmTMB(
  variance ~ trt*category_impact+(1|site_code/block),
  dispformula = ~ category_impact,
  data = df %>% filter(category_impact != "= Mean, = Variance")
)

plot(simulateResiduals(m_1))

summary(m_1)


library(emmeans)
library(multcomp)
library(multcompView)


df_res<-data.frame(cld(emmeans(m_1,~category_impact|trt),Letters = letters))


ggplot(df_res,
       aes(category_impact,emmean,fill=category_impact))+
  
  geom_col()+
  geom_errorbar(aes(ymax=upper.CL,ymin=lower.CL),
                width = 0.2)+
  facet_grid(~trt)+
  ylab("Dispersion")+
  ylim(0,0.7)+
  geom_text(aes(label = .group,x = category_impact, y = emmean),
            position = position_dodge(width = 2), vjust = -2.5)+
  theme_pubr()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank())



ggplot(df,aes(as.numeric(MAP_VAR_v2),RAIN_PET,colour = category_impact))+
  geom_point(size=4)
