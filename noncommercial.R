# Alberto Rovellini
# 3/28/2024
# Document to explore the non-commercial catch data set from AKFIN Answers:
# A collection of non-commercial fishery catch sources to include subsistence, sport, and survey data as compiled by NMFS Alaska Region: Sourced by the AKR.V_NONCOMMERCIAL_FISHERY_CATCH table.

library(tidyverse)

dat <- read.csv("../../Catch_data/data/AKFIN/Non_Commercial/Noncommercial Fishery Catch.csv", fileEncoding = "UTF-8-BOM")
gear_key <- read.csv("../../Catch_data/data/AKFIN/Fish_Tickets/Target_codes_BLEND.csv", fileEncoding = "UTF-8-BOM")
sp_key <- read.csv("../../Catch_data/data/eLandings_Atlantis_lookup.csv")

glimpse(dat)

# let's have a look at some of the fields
sort(unique(dat$Collection.Agency)) # "ADFG" "IPHC" "NMFS"
sort(unique(dat$Collection.Name))
# a lot of the collection names are from surveys of all kinds and from different agencies throughout the GOA
# The most promising-sounding of these are: 
# "Sport Fishery" "Subsistence Fishery"
# Should probably filter the df to those two only
dat_fishery <- dat %>%
  filter(Collection.Name %in% c("Sport Fishery", "Subsistence Fishery"))

# it's a rather small data set
sort(unique(dat_fishery$Collection.Year))
sort(unique(dat_fishery$OBS.Gear.Code)) # empty
sort(unique(dat_fishery$Reporting.Area)) # 610 620 630 640 649 650 659
sort(unique(dat_fishery$ADFG.Statistical.Area)) # empty
sort(unique(dat_fishery$Species.Group.Name)) # 21 species

# there is no gear nor ADFG stat area information in this data
# The former is not too important but the latter means that we cannot allocate these catches in space other than using the NMFS areas

# let's do some exploratory plotting
dat_fishery %>%
  group_by(Collection.Year, Collection.Name, Species.Group.Name) %>%
  summarize(catch = sum(Weight)) %>%
  ggplot(aes(x = Collection.Year, y = catch, fill = Species.Group.Name))+
  geom_area()+
  theme_bw()+
  facet_wrap(~Collection.Name)

# the problem is apparent: "Subsistence" fishery landings end in 2005; sport fishery landings begin in 2010
# I am unsure why this happens but it is not particularly useful data

