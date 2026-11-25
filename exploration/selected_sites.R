# Select sites

library(tidyverse)

dat <- read.csv("./nutnet_data/full-cover-2025-12-09.csv",sep=";")


# 1. paired blocks (Control & NPK - same block, same year)

paired <- dat %>%
  filter(trt %in% c("Control", "NPK")) %>%
  distinct(site_code, year_trt, block, trt) %>%
  group_by(site_code, year_trt, block) %>%
  filter(n_distinct(trt) == 2) %>%
  ungroup() %>%
  distinct(site_code, year_trt, block)

# 2. Valid years: at least 3 paired blocks

years_ok <- dat %>%
  filter(trt %in% c("Control", "NPK")) %>%
  count(site_code, year_trt, name = "n_blocks") %>%
  filter(n_blocks >= 3)

# 3. Site criteria:
#    - 10+ years of trt (year_trt max >= 10)
#    - no more than 2 missing years

keep_sites <- dat%>%
  filter(trt %in% c("Control", "NPK")) %>%
  group_by(site_code) %>%
  summarise(last_yr  = max(year_trt),
            n_years  = n(),
            max_run_missing = {
              r <- rle(0:max(year_trt) %in% year_trt)
              max(c(0, r$lengths[!r$values]))
            },
            .groups = "drop") %>%
  filter(last_yr >= 10,
         max_run_missing <= 2)

# 4. Final selection

dat_final <- dat %>%
  filter(trt %in% c("Control", "NPK")) %>%
  semi_join(keep_sites, by = "site_code") %>%
  semi_join(years_ok,   by = c("site_code", "year_trt")) %>%
  semi_join(paired,     by = c("site_code", "year_trt", "block")) %>%
  
  #keep only site_code names
  distinct(site_code)%>%
  
  # Doane keeps only 1 block after year_trt 1
  filter(site_code!="doane.us")

write.csv(dat_final,"./my_data/selected_sites.csv")

