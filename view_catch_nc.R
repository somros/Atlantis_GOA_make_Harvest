# Alberto Rovellini
# 6/7/2024
# this script extracts and views catch in space from netcdf files
# allow for comparisons between runs (e.g., Base no MPA, MPA fully open, in the future partially closed MPAs, etc)
# allows for comparisons between a run and Adam's catch streams (aka "the truth")
# gets residuals between compared runs
# TODO: add to the function an option of getting a time series instead of the terminal average
# In fact, that should just be the default, and then any averaging could be done afterwards

library(tidync)
library(ncdf4)
library(tidyverse)
library(sf)
library(rbgm)
library(viridis)

select <- dplyr::select

# nc files
# base
catch_nc_file_base <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1517/outputGOA01517_testCATCH.nc"
bio_nc_file_base <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1517/outputGOA01517_test.nc"
# mpa
catch_nc_file_mpa <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1553/outputGOA01553_testCATCH.nc"
bio_nc_file_mpa <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1553/outputGOA01553_test.nc"
# groups
grps <- read.csv("data/GOA_Groups.csv")
all_fg <- grps %>% filter(GroupType %in% c("FISH", "SHARK")) %>% pull(Name)

# geometry
goa_bgm <- read_bgm("data/GOA_WGS84_V4_final.bgm")
goa_sf <- goa_bgm %>% box_sf()

# catch reconstruction
fleets <- readRDS("fleets/fleet_total_catch_atl.RDS")
fleet_key <- read.csv("data/GOA_fisheries.csv")

# process fishery data

# load function to extract catch from netcdf
source("build_catch_output.R")

# Compare between Atlantis runs -------------------------------------------
# Use v2
catch_nc_base <- build_catch_output_v2(catch_nc = catch_nc_file_base, 
                                    fleet_struc = F,
                                    relative = T,
                                    run = 1517,
                                    key = fleet_key)

catch_nc_mpa <- build_catch_output_v2(catch_nc = catch_nc_file_mpa, 
                                   fleet_struc = F,
                                   relative = T,
                                   run = 1555,
                                   key = fleet_key)

# get residuals (it gets hazy for proportions - what's a big residual and how do you translate that to catch?)
catch_diff <- catch_nc_base %>%
  left_join(catch_nc_mpa %>% st_set_geometry(NULL), by = c("box_id", "Name")) %>%
  mutate(residual_mt = mt_tot.x - mt_tot.y,
         residual_prop = prop.x - prop.y)
  
# view
catch_diff %>%
  filter(Name == "Pollock") %>%
  ggplot()+
  geom_sf(aes(fill = residual_prop))+
  scale_fill_viridis()+
  theme_bw()+
  facet_wrap(~Name)

# These plots are difficult to interpret

# Compare between Atlantis run and data -----------------------------------

catch_atlantis <- build_catch_output_v2(catch_nc = catch_nc_file_base, 
                                      fleet_struc = F,
                                      relative = T,
                                      run = 1517,
                                      key = fleet_key)

# average of end of the run
catch_atlantis_end <- catch_atlantis %>%
  filter(ts > (max(ts)-5)) %>%
  group_by(box_id, Name) %>%
  summarize(mt = mean(mt)) %>%
  ungroup()

# first no fleets
catch_data <- fleets %>%
  left_join(grps %>% select(Code,Name), by = c("spp"="Code")) %>%
  select(year, box_id, Name, fleet, weight_mton) %>%
  group_by(year, box_id, Name) %>%
  summarise(mt = sum(weight_mton)) %>% # sum up across fleets
  ungroup() %>%
  filter(year > (max(year)-5)) %>%
  group_by(box_id, Name) %>%
  summarize(mt = mean(mt)) %>%
  ungroup()

# check differences
catch_diff <- catch_data %>%
  left_join(catch_atlantis_end, by = c("box_id", "Name")) %>%
  mutate(residual = mt.x - mt.y) # observed - predicted

# add space
catch_diff <- goa_sf %>%
  select(box_id) %>%
  left_join(catch_diff, by = "box_id")

# view
# do only verts
catch_diff %>%
  filter(Name %in% all_fg) %>%
  ggplot()+
  geom_sf(aes(fill = residual))+
  scale_fill_viridis()+
  theme_bw()+
  facet_wrap(~Name)

# TODO: move this forward - need to add 0 boxes to the data, reassign island catches (check bc we have done it)
# Now it will be way out of sorts, but the final goal WHEN WE WORK WITH REAL F's will be to:
# 1. Look at catches in space in the reconstruction
# 2. Get a scalar of the ratio between what is caught and what should be caught
# 3. Rescale MPAYYY entries so that they force more / less catch
# I find this a better hack (in theory) than modifying mFC - thogh I am sure in practice it won't be as easy
