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

# for each species, get proportion of total catch by fleet
fleets_prop <- fleets_tot %>%
  group_by(year, spp) %>%
  mutate(tot_mt = sum(mt)) %>%
  ungroup() %>%
  mutate(prop = mt / tot_mt)

# is this summing up to 1?
check <- fleets_prop %>%
  group_by(year, spp) %>%
  summarise(check = sum(prop))
# yes

# pick a point in time:
# say average 2015-2020
# This assumes that proportional effects of different fleets is time-invariant: pot caught x% of cod in 2020, so it will in 2080, so it did in 1995.

