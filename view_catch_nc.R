# Alberto Rovellini
# 6/7/2024
# this script extracts and views catch in space from netcdf files
# allow for comparisons between runs (e.g., Base no MPA, MPA fully open, in the future partially closed MPAs, etc)
# allows for comparisons between a run and Adam's catch streams (aka "the truth")
# gets residuals between compared runs

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
catch_nc_file_mpa <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1555/outputGOA01555_testCATCH.nc"
bio_nc_file_mpa <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1555/outputGOA01555_test.nc"
# groups
grps <- read.csv("data/GOA_Groups.csv")
all_fg <- grps %>% filter(GroupType %in% c("FISH", "SHARK")) %>% pull(Name)

# geometry
goa_bgm <- read_bgm("data/GOA_WGS84_V4_final.bgm")
goa_sf <- goa_bgm %>% box_sf()

# fishery data

# process fishery data

# load function to extract catch from netcdf
source("build_catch_output.R")

# Compare between Atlantis runs -------------------------------------------
catch_nc_base <- build_catch_output(catch_nc = catch_nc_file_base, 
                                    bio_nc = bio_nc_file_base, 
                                    age_struc = F,
                                    relative = T,
                                    run = 1517)

catch_nc_mpa <- build_catch_output(catch_nc = catch_nc_file_mpa, 
                                   bio_nc = bio_nc_file_mpa,
                                   age_struc = F,
                                   relative = T,
                                   run = 1555)

# get residuals (it gets hazy for proportions - what's a big residual and how do you translate that to catch?)
catch_diff <- catch_nc_base %>%
  left_join(catch_nc_mpa %>% st_set_geometry(NULL), by = c("box_id", "Name")) %>%
  mutate(residual_mt = mt_tot.x - mt_tot.y,
         residual_prop = prop.x - prop.y)
  
# view
catch_diff %>%
  filter(Name == "Salmon_pink") %>%
  ggplot()+
  geom_sf(aes(fill = residual_prop))+
  scale_fill_viridis()+
  theme_bw()+
  facet_wrap(~Name)

# When the MPA cells are fully open, the spatial distribution seems to be similar
# However, the total catch is much lower, so why is that the case?


# Compare between Atlantis run and data -----------------------------------





# Diagnostic with CSV files -----------------------------------------------
# catch seems much lower in the file with MPAs
# This can also be seen in the CSV files
# catch_flat_base <- read.table("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1553/outputGOA01553_testCatch.txt", sep = " ", header = T)
# catch_flat_mpa <- read.csv("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1558/outputGOA01558_testCatch.txt", sep = " ", header = T)
# 
# catch_flat_base_long <- catch_flat_base %>%
#   pivot_longer(-Time, names_to = "Code", values_to = "mt") %>%
#   mutate(run = "base")
# 
# catch_flat_mpa_long <- catch_flat_mpa %>%
#   pivot_longer(-Time, names_to = "Code", values_to = "mt") %>%
#   mutate(run = "mpa")
# 
# catch_flat <- rbind(catch_flat_base_long, catch_flat_mpa_long)
# 
# # view
# catch_flat %>%
#   filter(Code %in% (grps %>% filter(Name %in% all_fg) %>% pull(Code)))  %>%
#   ggplot(aes(x = Time, y = mt, color = run))+
#   geom_line()+
#   theme_bw()+
#   facet_wrap(~Code, scales = "free")
