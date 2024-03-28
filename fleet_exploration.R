# Alberto Rovellini
# 03/26/2024
# This document explores the fleet structure that Adam developed from metier analysis of GOA fish tickets
# it also identifies possible gaps in the information
# Things to look for:
# Fleet names
# Flet target species
# Annual catches per fleet?
# Links with fish ticket / catch data (gear, target, port, etc...). Likely will need to do a new pull from AKFIN
# Spatial information - can we derive a footprint from this data?

library(tidyverse)
library(rbgm)
library(sf)
library(maps)
library(mapdata)
library(lubridate)
library(data.table)
library(viridis)

select <- dplyr::select


# Read in data ------------------------------------------------------------

# read in RDS object from Adam
# read fleets and fleet keys
fleets <- readRDS("fleets/fleet_total_catch_atl.RDS")
fleet_key <- read.csv("fleets/fleet_Atlantis.csv")
grps <- read.csv("data/GOA_Groups.csv")

# read in spatial data
goa_bgm <- read_bgm("data/GOA_WGS84_V4_final.bgm")
goa_sf <- goa_bgm %>% box_sf()
st_crs(goa_sf) <- st_crs(attr(goa_sf$geometry, "crs")$proj)
boundary_boxes <- goa_sf %>% filter(boundary == TRUE) %>% pull(box_id) # get boundary boxes

# coast shapefile
coast <- maps::map(database = 'worldHires', regions = c('USA','Canada'), plot = FALSE, fill=TRUE)

coast_sf <- coast %>% 
  st_as_sf(crs = 4326) %>% 
  st_transform(crs = st_crs(goa_sf)) %>% 
  st_combine() %>%
  st_crop(goa_sf %>% st_bbox())

# join 
fleets <- fleets %>%
  left_join(fleet_key, by = c("fleet"="Fleet")) %>%
  left_join(grps %>% select(Code, LongName), by = c("spp"="Code"))


# Summary table of the 21 fleets ------------------------------------------

# make a summary table for the fleets, ordered by recent catch
fleets_summary <- fleets %>%
  filter(year > 2015) %>%
  group_by(year, fleet, Primary.Spp, Gear, Fishing.Area, Landing.Area) %>%
  summarize(catch_mt = sum(weight_mton)) %>%
  group_by(fleet, Primary.Spp, Gear, Fishing.Area, Landing.Area) %>%
  summarize(catch_mt = mean(catch_mt)) %>%
  ungroup() %>%
  select(Primary.Spp, Gear, Fishing.Area, Landing.Area, catch_mt, fleet) %>%
  arrange(-catch_mt)

colnames(fleets_summary) <- c("Primary target", "Gear", "Fishing area", "Landing area", "Mean catch 2016-2021 (mt)")

# write.csv(fleets_summary, "fleets/fleets_summary_by_recent_catch.csv", row.names = F)


# Total catch by fleet of each species, no areas -------------------------------
# for this to be informative, let's break it down into one plot per fleet, with only the needed species
fleet_names <- unique(fleets$fleet)

for(i in 1:length(fleet_names)){
  
  # pick out a fleet
  this_fleet_name <- fleet_names[i]
  
  # filter df to that fleet
  this_fleet <- fleets %>%
    filter(fleet == this_fleet_name, year > 2007)
  
  # define gear, target, fishing area, landing area
  this_target <- this_fleet %>% pull(Primary.Spp) %>% unique()
  this_gear <- this_fleet %>% pull(Gear) %>% unique()
  this_fishing_area <- this_fleet %>% pull(Fishing.Area) %>% unique()
  this_landing_area <- this_fleet %>% pull(Landing.Area) %>% unique()
  
  area_chart_by_fleet <-  this_fleet %>%
    group_by(year, LongName) %>%
    summarize(catch_mt = sum(weight_mton)) %>%
    ggplot(aes(x = year, y = catch_mt, fill = LongName))+
    geom_area()+
    theme_bw()+
    labs(title = paste0("Target = ", this_target, "\n",
                       "Gear = ", this_gear, "\n",
                       "Fishing area = ", this_fishing_area, "\n",
                       "Landing area = ", this_landing_area))+
    guides(fill = guide_legend(ncol = 2))
  
  # save plot
  ggsave(paste0("fleets/area_chart_by_fleet/", this_fleet_name, ".png"), area_chart_by_fleet, width = 8, height = 6)
  
}

# Total catch by species irrespective of fleets and areas -----------------

fleets_total <- fleets %>%
  group_by(year, LongName) %>%
  summarize(catch_mt = sum(weight_mton))
  
fleets_total_plot <- fleets_total %>%
  ggplot(aes(x = year, y = catch_mt))+
  geom_line()+
  geom_point()+
  theme_bw()+
  facet_wrap(~LongName, scales = "free")

# Spatial analysis --------------------------------------------------------

# map catches in space for the top key fleets
# make a summary table for the fleets, ordered by recent catch
fleets_spatial <- fleets %>%
  filter(year > 2015) %>%
  group_by(year, fleet, Primary.Spp, Gear, Fishing.Area, Landing.Area, box_id) %>%
  summarize(catch_mt = sum(weight_mton)) %>%
  group_by(fleet, Primary.Spp, Gear, Fishing.Area, Landing.Area, box_id) %>%
  summarize(catch_mt = mean(catch_mt)) %>%
  ungroup()

# make it spatial
fleets_spatial <- goa_sf %>%
  select(box_id) %>%
  left_join(fleets_spatial, by = "box_id")

for(i in 1:length(fleet_names)){
  
  # pick out a fleet
  this_fleet_name <- fleet_names[i]
  
  # filter df to that fleet
  this_fleet <- fleets_spatial %>%
    filter(fleet == this_fleet_name)
  
  # define gear, target, fishing area, landing area
  this_target <- this_fleet %>% pull(Primary.Spp) %>% unique()
  this_gear <- this_fleet %>% pull(Gear) %>% unique()
  this_fishing_area <- this_fleet %>% pull(Fishing.Area) %>% unique()
  this_landing_area <- this_fleet %>% pull(Landing.Area) %>% unique()
  
  map_by_fleet <-  ggplot()+
    geom_sf(data = this_fleet, aes(fill = catch_mt))+
    geom_sf(data = coast_sf)+
    scale_fill_viridis()+
    theme_bw()+
    labs(title = paste0("Target = ", this_target, "\n",
                        "Gear = ", this_gear, "\n",
                        "Fishing area = ", this_fishing_area, "\n",
                        "Landing area = ", this_landing_area))
  
  # save plot
  ggsave(paste0("fleets/maps_by_fleet_recent_catch/", this_fleet_name, ".png"), map_by_fleet, width = 8, height = 4)
  
}

# allocation bled into the boundary boxes, either because I did not explain that they should have been ignored or because
# the ADF&G statistical areas (squares) overlap with them. 
# Either way, we will need to handle those either by just dropping catches in those boxes, or by re-allocating it to neighboring boxes

# Port allocation ---------------------------------------------------------

# The data from Adam has port information. 
unique(fleets$port)

# Can we, for each fleet, allocate what went to which port?
fleets_by_port <- fleets %>%
  filter(year > 2015) %>%
  group_by(year, fleet, Primary.Spp, Gear, Fishing.Area, Landing.Area, port, LongName) %>%
  summarize(catch_mt = sum(weight_mton)) %>%
  group_by(fleet, Primary.Spp, Gear, Fishing.Area, Landing.Area, port, LongName) %>%
  summarize(catch_mt = mean(catch_mt)) %>%
  ungroup()

# plot
for(i in 1:length(fleet_names)){
  
  # pick out a fleet
  this_fleet_name <- fleet_names[i]
  
  # filter df to that fleet
  this_fleet <- fleets_by_port %>%
    filter(fleet == this_fleet_name)
  
  # define gear, target, fishing area, landing area
  this_target <- this_fleet %>% pull(Primary.Spp) %>% unique()
  this_gear <- this_fleet %>% pull(Gear) %>% unique()
  this_fishing_area <- this_fleet %>% pull(Fishing.Area) %>% unique()
  this_landing_area <- this_fleet %>% pull(Landing.Area) %>% unique()
  
  bar_ports <-  this_fleet %>%
    ggplot()+
    geom_bar(aes(x = port, y = catch_mt, fill = LongName), stat = "identity", position = "stack")+
    theme_bw()+
    labs(title = paste0("Target = ", this_target, "\n",
                        "Gear = ", this_gear, "\n",
                        "Fishing area = ", this_fishing_area, "\n",
                        "Landing area = ", this_landing_area))+
    guides(fill = guide_legend(ncol = 2))
  
  # save plot
  ggsave(paste0("fleets/catch_by_port/", this_fleet_name, ".png"), bar_ports, width = 8, height = 6)
  
}

# Comparison with catch reconstruction for catch.ts -----------------------

# From Adam:
# The spatial catch data represented by fish tickets and catch accounting data are statistical 
# areas defined by the ADFG, which differ from the spatial boxes used in the Atlantis model. 
# Catch in ADFG statistical areas are converted to Atlantis box by calculating the overlap 
# between the ADFG area and each Atlantis box. The catch is apportioned to the Atlantis boxes 
# according to the proportion of the statistical area the box overlaps with; e.g., if a given 
# Atlantis box overlaps 50% with an ADFG statistical area, then 50% of the catch in that statistical 
# area will be attributed to that Atlantis box. This may create issues for boundary boxes where little 
# fishing occurs if an ADFG statistical area with significant fishing overlaps with it, creating an erroneous spatial imputation of catch. 

# This suggests that the numbers in the RDS object from Adam should be total catch by fleet accounting for both tickets and catch accounting data
# But, to get a somewhat realistic F, we'd need to ensure that this captures global catch to satisfation
# do some vieweing and compare to:
# 1. stock assessment catch reconstructions
# 2. catch time series we did

# We will need to turn these into F's. If they are in line with the catch.ts data and/or the stock assessment catch reconstruction,
# we'll need to get best estimates of stock status (we have those) and get F based on catch (remember selex here)
# how do we parameterize selectivity for these?

# Adam recommends that the evaluation of F is done with data after 2008
# read in time series

# Atlantis groups for column names
all_fg <- grps %>% filter(IsImpacted > 0) %>% pull(Code)

# list forcing files
details_ak <- file.info(list.files('../../Catch_data/output/AKFIN/', full.names = T))
details_ak <- details_ak[with(details_ak, order(as.POSIXct(mtime))), ]
files_ak <- rownames(details_ak)

details_bc <- file.info(list.files('../../Catch_data/output/DFO/', full.names = T))
details_bc <- details_bc[with(details_bc, order(as.POSIXct(mtime))), ]
files_bc <- rownames(details_bc)

# read all files
all_boxes <- 0:108

# separate AK boxes from BC boxes
ak_boxes <- all_boxes[all_boxes < 92]
bc_boxes <- all_boxes[all_boxes > 91]

# Alaska first

ak_data <- list()

for(b in 1:length(ak_boxes)){
  this_box <- ak_boxes[b]
  this_tab <- read.table(files_ak[b], skip = 309)
  colnames(this_tab) <- c('Time', all_fg) # add column names
  this_tab <- this_tab %>% mutate(box_id = this_box) # add a column for the box number
  ak_data[[b]] <- this_tab
} 

all_ak <- bind_rows(ak_data)

# not doing canada because Adam's data is from AK

# together
all_catch <- all_ak

# transform time steps into dates
this_origin <- as.Date('1991-01-01', tz = 'UTC')

# change format, keep box at first
all_catch_long_box <- all_catch %>%
  mutate(Date = as.Date(Time, "%Y-%m-%d", origin = this_origin, tz = 'UTC')) %>%
  select(Date, box_id, KWT:BIV) %>%
  pivot_longer(-c(Date, box_id), names_to = 'Species', values_to = 'Catch_mgs') %>%
  mutate(Catch_mt_day = Catch_mgs * 60 * 60 * 24 * 20 * 5.7 / 1e9,
         Catch_mt_month = Catch_mt_day * 30, # THIS IS AN APPROXIMATION FOR VISUALISATION PURPOSES
         Year = year(Date)) %>% 
  select(Year, Date, box_id, Species, Catch_mt_month) %>%
  distinct() %>%
  group_by(Year, box_id, Species) %>%
  mutate(Catch_mt_year = sum(Catch_mt_month)) %>% # sum over months for annual catch
  ungroup() 

all_catch_long <- all_catch_long_box %>%
  select(Year, box_id, Species, Catch_mt_year) %>%
  distinct() %>%
  group_by(Year, Species) %>%
  summarise(Catch_mt_year = sum(Catch_mt_year)) %>% # sum across boxes
  ungroup() 

# join the two sets
# reconstruction from Adam
fleets_total <- fleets_total %>%
  mutate(set = "adam")

catch_total <- all_catch_long %>%
  left_join(grps %>% select(Code, LongName), by = c("Species"="Code")) %>%
  select(Year, LongName, Catch_mt_year) %>%
  rename(year = Year, catch_mt = Catch_mt_year) %>%
  mutate(set = "ts")

# join
comp <- rbind(fleets_total, catch_total)

# view
comp_plot <- comp %>%
  filter(year >= 2008) %>%
  ggplot(aes(x = year, y = catch_mt, color = set))+
  geom_point()+
  geom_line()+
  theme_bw()+
  facet_wrap(~LongName, scales = "free")

# for most groups, the catch reconstruction to produce the catch.ts files seems to have higher catches
# there are several spots where this may have occurred
# For some species the old catch reconstruction used all data that we had available, so I would assume that this would be a better estimate of removals
# view for some specific groups

comp %>%
  filter(LongName %in% (grps %>% 
                          select(Code, LongName) %>% 
                          filter(Code %in% c("POL","ATF","COD","FFS","FHS","FFD","REX","POP","SBF","RFS","RFP","HAL")) %>% 
                          pull(LongName))) %>%
  filter(year >= 2008) %>%
  ggplot(aes(x = year, y = catch_mt, color = set))+
  geom_point()+
  geom_line()+
  theme_bw()+
  facet_wrap(~LongName, scales = "free")

# groundfish looks OK (T3 in particular, remembering which stocks are important here)
# salmon is more problematic
comp %>%
  filter(LongName %in% (grps %>% 
                          select(Code, LongName) %>% 
                          filter(Code %in% c("SCH","SCM","SCO","SSO","SPI")) %>% 
                          pull(LongName))) %>%
  filter(year >= 2008) %>%
  ggplot(aes(x = year, y = catch_mt, color = set))+
  geom_point()+
  geom_line()+
  theme_bw()+
  facet_wrap(~LongName, scales = "free")

# chinook is close, the others are not, especially chum and sockeye
# part of the problem may be how we mapped salmon catch to the GOA
# in the original reconstruction I mapped ADFG areas KMLAEH, but maybe Adam did not account for Cook Inlet and PWS (and the southeast)
# probably here it makes sense to stick with what Adam got, for consistency

# overall, it likely does not matter so much, because this F will require calibration anyway, until it yields catch that is in the ballpark of the catch for the desired period
# my catch reconstruction only goes to 2020
# 
