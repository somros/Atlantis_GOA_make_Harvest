# Alberto Rovellini
# 4/5/2024
# This code gets proportion of the catch for each species by fleet
# It is based on Adam's metier analysis
# The idea is to use these proportions to break down F (or mFC), whichever mFC we will use
# Some applications:
# break down background mFC into mFC by fleet
# break down FMSY, Ftarg, Fref, etc. by fleet

library(tidyverse)

# Read in data ------------------------------------------------------------

# read in RDS object from Adam
# read fleets and fleet keys
fleets <- readRDS("fleets/fleet_total_catch_atl.RDS")
fleet_key <- read.csv("fleets/fleet_Atlantis.csv")
grps <- read.csv("data/GOA_Groups.csv")

# remove space (model-wide F)
fleets_tot <- fleets %>%
  group_by(fleet, year, spp) %>%
  summarise(mt = sum(weight_mton))

# expand with all missing combinations
template <- expand.grid(fleet = unique(fleets_tot$fleet), 
                        year = unique(fleets_tot$year),
                        spp = unique(fleets_tot$spp))
fleets_tot <- merge(template, fleets_tot, by = c("fleet","year","spp"), all.x = T)
fleets_tot$mt[is.na(fleets_tot$mt)] <- 0

# for each species, get proportion of total catch by fleet
fleets_prop <- fleets_tot %>%
  group_by(year, spp) %>%
  mutate(tot_mt = sum(mt)) %>%
  ungroup() %>%
  mutate(prop = mt / tot_mt)

# is this summing up to 1?
# check <- fleets_prop %>%
#   group_by(year, spp) %>%
#   summarise(check = sum(prop))
# yes

# pick a point in time:
# say average 2015-2020
# This assumes that proportional effects of different fleets is time-invariant: pot caught x% of cod in 2020, so it will in 2080, so it did in 1995.
fleets_prop_term <- fleets_prop %>%
  filter(year >= 2015, year <= 2020) %>%
  group_by(spp, fleet) %>%
  summarise(prop = mean(prop, na.rm = T))

# check
# check <- fleets_prop_term %>%
#   group_by(spp) %>%
#   summarise(check = sum(prop))
# okay

# view
fleets_prop_term %>%
  ggplot(aes(x = spp, y = prop, fill = fleet))+
  geom_bar(stat = "identity", position = "stack")+
  theme_bw()

# write this out
write.csv(fleets_prop_term, "data/mFC_prop_by_fleet.csv", row.names = F)

# Apply to background mFC -------------------------------------------------

# open the harvest.prm from the Base model
# Only for the species that appear here, break down mFC across fleets
# make sure that the order of the fleets is consistent with the fleets.csv input file

file_path <- "C:/Users/Alberto Rovellini/Documents/GOA/Atlantis_GOA_OY_MS/"
bg_harvest_file <- paste0(file_path, "GOA_harvest_background.prm")
goa_fleets_file <- paste0(file_path, "GOA_fisheries.csv")

# read in
bg_harvest <- readLines(bg_harvest_file)
goa_fisheries <- read.csv(goa_fleets_file)
