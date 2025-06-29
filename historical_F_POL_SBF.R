# get F time series for pollock and sablefish
# we use model fits from Cole Monnahan (pollock) and Matt Cheng (sablefish) for the 2024 assessments
# Atlantis will use the HCRs for both stocks, but we need historical F for some scenarios (including other GOACLIM models)

library(tidyverse)

rm(list = ls())


# Sablefish ---------------------------------------------------------------
sbf <- readRDS("data/sablefish_assessment.RDS")

# Fmort is dimensioned by region, years, and fleet
# Region 1 and 2 are BS and AI
# Region 3, 4, and 5 are WGOA, CGOA, and EGOA
# Years 1 - 62 = 1960 - 2021
# Fleet 1 = Longline, Fleet 2 = Trawl
# Dan Goethel:  I'd probably say that a total F of 0.05 - 0.10 would be in the ballpark for sablefish in recent years (as you can see in table 3.10)
# Approach: get F for WGOA, CGOA, EGOA
# Sum across fleets
# Average over areas to get GOA-wide F
# This is approximate but in the ballpark
f_sbf <- sbf$Fmort
f_sbf <- apply(f_sbf[3:5,,], c(1,2), sum)
f_sbf <- f_sbf %>% 
  t() %>% 
  as.data.frame() %>% 
  set_names(c("WGOA","CGOA","EGOA")) %>% 
  mutate(year = 1960:2021) %>%
  pivot_longer(-year, names_to = "area", values_to = "f")

# pull selected biomass to weigh these F's between the areas
# weigh by total biomass estimated for each area. Exploitable biomass would have been better
# Region 1 and 2 are BS and AI
# Region 3, 4, and 5 are WGOA, CGOA, and EGOA
# Years 1 - 62 = 1960 - 2021
biom_sbf <- sbf$Total_Biom[3:5,] %>%
  t() %>% 
  as.data.frame() %>% 
  set_names(c("WGOA","CGOA","EGOA")) %>% 
  mutate(year = 1960:2021) %>%
  pivot_longer(-year, names_to = "area", values_to = "biom")

f_sbf <- f_sbf %>%
  left_join(biom_sbf, by = c("year","area")) %>%
  group_by(year) %>%
  summarise(mean_f = weighted.mean(f,biom)) %>%
  pull(mean_f)

# Pollock -----------------------------------------------------------------
pol <- readRDS("data/pollock_assessment.RDS")
pol_f <- pol$sd %>% filter(name=='F') %>% select(est,year)

# how does this compare to Stock SMART?
pol_f_ss <- f_ss %>% filter(grepl("Walleye pollock",`Stock Name`)) %>% select(Year,Value) %>% mutate(Year = as.numeric(Year), Value = as.numeric(Value))
pol_f %>% 
  left_join(pol_f_ss, by = c("year"="Year")) %>%
  mutate(comp = est / Value) %>%
  pull(comp) %>%
  summary() # fairly close overall

# How does this compare to the exloitation rates in the assessment?
pol_f %>% filter(year >= 1977) %>% mutate(mu = 1 - exp(-est)) %>% pull(mu) # pretty poorly

# TODO: figure out the difference between F and exploitation rate in the assessments
# For now use F as it seems to be more consistent with the rest of the assessments
pol_f %>%
  filter(year > 1976) %>%
  pull(est)



