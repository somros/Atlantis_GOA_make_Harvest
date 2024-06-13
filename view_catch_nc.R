# Alberto Rovellini
# 6/7/2024
# this script extracts and views catch in space from netcdf files
# allow for comparisons between runs (e.g., Base no MPA, MPA fully open, in the future partially closed MPAs, etc)
# allows for comparisons between a run and Adam's catch streams (aka "the truth")
# gets residuals between compared runs
# TODO: add to the function an option of getting a time series instead of the terminal average
# In fact, that should just be the default, and then any averaging could be done afterwards

# Important caveat:
# Adam allocated catches to boundary boxes and island boxes, likely as a result of a spatial join between statistical areas and Atlantis boxes
# For island boxes, catch may come from any number of the boxes around an island
# For boundary boxes, catch may be a misallocation due to the chunky box shape or it could legitimately be catch occurring off the shelf
# At the moment (6/12/2024) it seems unlikely that we will be able to fully re-run his analysis, so we need to make assumtpions for both cases
# Island boxes: split catch among neighboring boxes based on the existing proportions
# boundary boxes: assume that this is catch that occurred outside the Atlantis domain and drop it from the reconstruction

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
catch_nc_file_mpa <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1562/outputGOA01562_testCATCH.nc"
bio_nc_file_mpa <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1562/outputGOA01562_test.nc"
# groups
grps <- read.csv("data/GOA_Groups.csv")
all_fg <- grps %>% filter(GroupType %in% c("FISH", "SHARK")) %>% pull(Name)

# geometry
goa_bgm <- read_bgm("data/GOA_WGS84_V4_final.bgm")
goa_sf <- goa_bgm %>% box_sf()

# catch data from reconstruction
source("handle_fleet_catch_reconstruction.R")

# load function to extract catch from netcdf
source("read_catch_nc_functions.R")

# Compare between Atlantis runs -------------------------------------------
# This is to compare the catch in space between two Atlantis runs
# Use v2
catch_nc_base <- build_catch_output_v2(catch_nc = catch_nc_file_base, 
                                       fleet_struc = F,
                                       relative = T,
                                       run = 1517,
                                       key = fleet_key)

catch_nc_mpa <- build_catch_output_v2(catch_nc = catch_nc_file_mpa, 
                                      fleet_struc = F,
                                      relative = T,
                                      run = 1562,
                                      key = fleet_key)

# get residuals (it gets hazy for proportions - what's a big residual and how do you translate that to catch?)
catch_diff <- catch_nc_base %>%
  left_join(catch_nc_mpa, by = c("ts", "box_id", "Name")) %>%
  mutate(residual_mt = mt.x - mt.y,
         residual_prop = prop.x - prop.y)

# add space
catch_diff <- goa_sf %>%
  select(box_id) %>%
  left_join(catch_diff, by = "box_id")

# view

to_plot <- c("Pollock","Cod","Arrowtooth_flounder","Flatfish_shallow","Flatfish_deep","Rex_sole","Flathead_sole","Pacific_ocean_perch","Sablefish","Halibut")

catch_diff %>%
  filter(ts == 15, Name %in% to_plot) %>%
  ggplot()+
  geom_sf(aes(fill = residual_mt))+
  scale_fill_viridis()+
  theme_bw()+
  facet_wrap(~Name)

# other view
catch_diff %>%
  filter(ts == 15, Name %in% to_plot) %>%
  ggplot(aes(x = mt.x, y = mt.y, color = box_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+  # Add 1:1 line
  labs(x = "Catch in non-spatial run", y = "Catch in spatial run")+
  facet_wrap(~Name, scales = "free")

# This step shows you the difference between the normal mFC Atlantis run and the run with the new MPA setup
# This offers some interesting avenues for calibration
# For example, one could derive box-specific indices of the catch ratio between base and MPA runs
# Problem here is that the MPAYYY vectors are fleet-specific, not species-specific

# Compare spatial patterns between Atlantis run and data -----------------------------------
# This only makes sense either in relative terms or when the model operates under realistic F

catch_atlantis <- build_catch_output_v2(catch_nc = catch_nc_file_mpa, 
                                        fleet_struc = F,
                                        relative = T,
                                        run = 1562,
                                        key = fleet_key)

# average of end of the run
catch_atlantis_end <- catch_atlantis %>%
  filter(ts > (max(ts)-5)) %>%
  group_by(box_id, Name) %>%
  summarize(mt = mean(mt)) %>%
  ungroup() %>%
  group_by(Name) %>%
  mutate(tot_mt = sum(mt)) %>%
  ungroup() %>%
  mutate(prop = mt / tot_mt)

# start wih no fleets here - space only
catch_data <- fleets %>%
  ungroup() %>%
  left_join(grps %>% select(Code,Name), by = c("spp"="Code")) %>%
  select(year, box_id, Name, fleet, weight_mton) %>%
  mutate(weight_mton = replace_na(weight_mton, 0)) %>%
  group_by(year, box_id, Name) %>%
  summarise(mt = sum(weight_mton)) %>% # sum up across fleets
  ungroup() %>%
  filter(year > (max(year)-5)) %>%
  group_by(box_id, Name) %>%
  summarize(mt = mean(mt)) %>%
  ungroup()

# add missing combinations
dummy_df <- expand.grid("box_id" = 0:108, "Name" = unique(catch_data$Name))
catch_data <- catch_data %>%
  full_join(dummy_df) %>%
  mutate(mt = replace_na(mt, 0))

# get relative catch
catch_data <- catch_data %>%
  group_by(Name) %>%
  mutate(tot_mt = sum(mt, na.rm = T)) %>%
  ungroup() %>%
  mutate(prop = mt / tot_mt)

# check differences
catch_diff <- catch_data %>%
  left_join(catch_atlantis_end, by = c("box_id", "Name"))# %>%
# mutate(residual = mt.x - mt.y,
#        ratio = mt.x / mt.y) # observed - predicted

# view
catch_diff <- catch_diff %>%
  select(box_id, Name, prop.x, prop.y) %>%
  rename(data = prop.x, model = prop.y) %>%
  pivot_longer(-c(box_id, Name), names_to = "Type", values_to = "Prop")

# add space
catch_diff <- goa_sf %>%
  select(box_id, boundary) %>%
  left_join(catch_diff, by = "box_id")

for(i in 1:length(all_fg)){
  this_fg <- all_fg[i]
  
  if(this_fg %in% unique(catch_diff$Name)){
    print(paste("doing", this_fg))
    
    p <- catch_diff %>%
      filter(Name == this_fg) %>%
      ggplot()+
      geom_sf(aes(fill = Prop))+
      scale_fill_viridis()+
      theme_bw()+
      facet_grid(Name~Type)
    
    ggsave(paste0("fleets/data_vs_catch/by_area/", this_fg, ".png"), p, width = 9, height = 3)
  }
  
}

# Now it will be way out of sorts, but the final goal WHEN WE WORK WITH REAL F's will be to:
# 1. Look at catches in space in the reconstruction
# 2. Get a scalar of the ratio between what is caught and what should be caught
# 3. Rescale MPAYYY entries so that they force more / less catch
# I find this a better hack (in theory) than modifying mFC - thogh I am sure in practice it won't be as easy

# Before we get to look at absolute catch, we can use this information to rescale catch in each box?
# We need to break this down by fleet

# Spatial patterns by fleet -----------------------------------------------

catch_atlantis_spatial_fleets <- build_catch_output_v2(catch_nc = catch_nc_file_mpa, 
                                                       fleet_struc = T,
                                                       relative = T,
                                                       run = 1562,
                                                       key = fleet_key)

# first - look at total catch across the model and species
# Is the split across fleets the same in the data and the model?
# deal with codes for the fleets being different
rewrite_codes <- function(original_string){
  # Split the string into words based on '_'
  words <- unlist(strsplit(original_string, "_"))
  
  # Convert each word to Title Case
  title_case_words <- sapply(words, function(word) {
    paste0(toupper(substr(word, 1, 1)), tolower(substr(word, 2, nchar(word))))
  })
  
  # Concatenate the words back together
  final_string <- paste0(title_case_words, collapse = "")
  
  return(final_string)
}

tot_atlantis <- catch_atlantis_spatial_fleets %>%
  filter(fleet != "background", fleet != "Canada") %>%
  filter(ts > (max(ts)-5)) %>%
  group_by(box_id, fleet) %>%
  summarize(mt = mean(mt)) %>%
  ungroup() %>%
  group_by(fleet) %>%
  summarize(mt = sum(mt)) %>%
  ungroup() %>%
  mutate(tot = sum(mt)) %>%
  mutate(prop = mt / tot) %>%
  mutate(Type = "model") %>%
  select(Type, fleet, prop)

bg <- catch_atlantis_spatial_fleets %>%
  filter(fleet == "background")

tot_data <- fleets %>%
  ungroup() %>%
  left_join(grps %>% select(Code,Name), by = c("spp"="Code")) %>%
  select(year, box_id, Name, fleet, weight_mton) %>%
  mutate(weight_mton = replace_na(weight_mton, 0)) %>%
  group_by(year, fleet) %>%
  summarise(mt = sum(weight_mton)) %>% # sum up across fleets
  ungroup() %>%
  filter(year > (max(year)-5)) %>%
  group_by(fleet) %>%
  summarize(mt = mean(mt)) %>%
  ungroup() %>%
  mutate(tot = sum(mt),
         prop = mt / tot) %>%
  rowwise() %>%
  mutate(fleet = rewrite_codes(fleet))%>%
  mutate(Type = "data") %>%
  select(Type, fleet, prop)

catch_diff <- rbind(tot_atlantis, tot_data) %>%
  left_join(fleet_key %>% select(Code, Name), by = c("fleet"="Code"))

# view
# drop empty fleets and background F
to_drop <- c("Canada",
             "background",
             unique(catch_diff$fleet)[grepl("dummy", unique(catch_diff$fleet))])
to_keep <- setdiff(unique(catch_diff$fleet), to_drop)

# view
catch_diff %>%
  filter(fleet %in% to_keep) %>%
  ggplot(aes(x = Type, y = prop, fill = Name))+
  geom_bar(stat = "identity", position = "stack")

# Now produce spatial plots but fleet-by-fleet
# model output
catch_atlantis_spatial_fleets_end <- catch_atlantis_spatial_fleets %>%
  filter(ts > (max(ts)-5)) %>% # keep last 5 years of run
  group_by(box_id, Name, fleet) %>%
  summarize(mt = mean(mt)) %>% # means across last 5 years
  group_by(box_id, fleet) %>%
  summarise(mt = sum(mt)) %>% # sum across species
  group_by(fleet) %>%
  mutate(tot_mt = sum(mt)) %>% # get total annual catch from a fleet
  ungroup() %>% 
  mutate(Prop = mt / tot_mt) %>% # get proportions by box
  mutate(Type = "model") %>%
  select(box_id, fleet, Prop, Type)

# data
catch_data_fleets <- fleets %>%
  ungroup() %>%
  left_join(grps %>% select(Code,Name), by = c("spp"="Code")) %>%
  select(year, box_id, Name, fleet, weight_mton) %>%
  mutate(weight_mton = replace_na(weight_mton, 0)) %>%
  filter(year >= (max(year)-5)) %>%
  group_by(box_id, Name, fleet) %>%
  summarize(mt = mean(weight_mton)) %>% # averages over time
  group_by(box_id, fleet) %>%
  summarise(mt = sum(mt)) %>% # sum across species
  ungroup()

# add missing combinations
dummy_df <- expand.grid("box_id" = 0:108, "fleet" = unique(catch_data_fleets$fleet))
catch_data_fleets <- catch_data_fleets %>%
  full_join(dummy_df) %>%
  mutate(mt = replace_na(mt, 0))

# get relative catch
catch_data_fleets <- catch_data_fleets %>%
  group_by(fleet) %>%
  mutate(tot_mt = sum(mt)) %>% # get total annual catch from a fleet
  ungroup() %>%
  mutate(Prop = mt / tot_mt) %>% # get proportions by box
  mutate(Type = "data") %>%
  select(box_id, fleet, Prop, Type) %>%
  rowwise() %>%
  mutate(fleet = rewrite_codes(fleet)) %>%
  ungroup()

# produce plots - these will be a lot to unpack
catch_diff <- rbind(catch_atlantis_spatial_fleets_end, catch_data_fleets) %>%
  left_join(fleet_key %>% select(Code, Name), by = c("fleet" = "Code"))

# add space
catch_diff_sf <- goa_sf %>%
  select(box_id) %>%
  left_join(catch_diff, by = "box_id")

# plots
for(i in 1:length(to_keep)){
  this_fleet <- to_keep[i]
  this_fleet_name <- catch_diff_sf %>%
    filter(fleet == this_fleet) %>%
    pull(Name) %>%
    unique()
  
  print(paste("doing", this_fleet))
  
  p <- catch_diff_sf %>%
    filter(fleet == this_fleet) %>%
    ggplot()+
    geom_sf(aes(fill = Prop))+
    scale_fill_viridis()+
    theme_bw()+
    facet_grid(~Type)+
    labs(title = this_fleet_name)
  
  ggsave(paste0("fleets/data_vs_catch/by_area_and_fleet/", this_fleet, ".png"), p, width = 9, height = 3)
  
}

# Fleet makeup of a species' catch --------------------------------------------------
# Check that catch split by fleet is similar to the data
# Because the data is only for AK, drop BC from the catch

catch_atlantis_species <- build_catch_output_v2(catch_nc = catch_nc_file_mpa, 
                                                fleet_struc = T,
                                                relative = T,
                                                run = 1562,
                                                key = fleet_key)

# average of end of the run
catch_atlantis_species_end <- catch_atlantis_species %>%
  filter(ts > (max(ts)-5)) %>%
  filter(box_id < 92) %>% # keep AK only
  group_by(fleet, Name) %>%
  summarize(mt = mean(mt)) %>%
  ungroup() %>%
  group_by(Name) %>% # group by species
  mutate(tot_mt = sum(mt)) %>%
  ungroup() %>%
  mutate(prop = mt / tot_mt)

# no space - fleets only
catch_data <- fleets %>%
  ungroup() %>%
  left_join(grps %>% select(Code,Name), by = c("spp"="Code")) %>%
  select(year, box_id, Name, fleet, weight_mton) %>%
  mutate(weight_mton = replace_na(weight_mton, 0)) %>%
  group_by(year, fleet, Name) %>%
  summarise(mt = sum(weight_mton)) %>% # sum up across boxes
  ungroup() %>%
  filter(year > (max(year)-5)) %>%
  group_by(fleet, Name) %>%
  summarize(mt = mean(mt)) %>%
  ungroup()

# add missing combinations
dummy_df <- expand.grid("fleet" = unique(catch_data$fleet), "Name" = unique(catch_data$Name))
catch_data <- catch_data %>%
  full_join(dummy_df) %>%
  mutate(mt = replace_na(mt, 0))

# get relative catch
catch_data <- catch_data %>%
  group_by(Name) %>%
  mutate(tot_mt = sum(mt, na.rm = T)) %>%
  ungroup() %>%
  mutate(prop = mt / tot_mt)

# deal with code names
catch_data <- catch_data %>% 
  rowwise() %>%
  mutate(fleet = rewrite_codes(as.character(fleet))) %>%
  ungroup()

# stack
catch_data <- catch_data %>% mutate(Type = "data")
catch_atlantis_species_end <- catch_atlantis_species_end %>% mutate(Type = "model")
catch_diff <- rbind(catch_data, catch_atlantis_species_end)
catch_diff <- catch_diff %>%
  left_join(fleet_key %>% select(Code, Name), by = c("fleet" = "Code"))

# view
# drop empty fleets and background F, we do not need to plot those
to_drop <- c("Canada",
             "background",
             unique(catch_diff$fleet)[grepl("dummy", unique(catch_diff$fleet))])
to_keep <- setdiff(unique(catch_diff$fleet), to_drop)

# make a plot per species
for(i in 1:length(all_fg)){
  this_fg <- all_fg[i]
  
  if(this_fg %in% unique(catch_data$Name)){
    
    print(paste("Doing", this_fg))
    
    p <- catch_diff %>%
      filter(Name.x == this_fg) %>%
      filter(fleet %in% to_keep) %>%
      filter(prop > 0) %>%
      ggplot(aes(x = Type, y = prop, fill = Name.y))+
      geom_bar(stat = "identity", position = "stack")+
      theme_bw()+
      facet_wrap(~Name.x)
    
    ggsave(paste0("fleets/data_vs_catch/by_species/", this_fg, ".png"), p, width = 8, height = 4.5)
    
  }
  
}

# the proportion of total catch by fleet works well
# Though this is the easy part - you may have a species that is barely caught at all being caught in the right proportions by different fleets
# Something that may be informative in that sense is the species makeup of the catch from a fleet, in relative terms
# I.e., looking at fleet Y, how does the species composition of the catch change?


# Species makeup of a fleet's catch ----------------------------------------
# simlar approach to above
# TODO: these should all become functions for easy deployment every time you have a new run to look at
# Also need to find a way to version the plots

# average of end of the run
catch_atlantis_fleets_end <- catch_atlantis_species %>%
  filter(ts > (max(ts)-5)) %>%
  filter(box_id < 92) %>% # keep AK only
  group_by(fleet, Name) %>%
  summarize(mt = mean(mt)) %>%
  ungroup() %>%
  group_by(fleet) %>% # group by fleets
  mutate(tot_mt = sum(mt)) %>%
  ungroup() %>%
  mutate(prop = mt / tot_mt)

# no space - fleets only
catch_data <- fleets %>%
  ungroup() %>%
  left_join(grps %>% select(Code,Name), by = c("spp"="Code")) %>%
  select(year, box_id, Name, fleet, weight_mton) %>%
  mutate(weight_mton = replace_na(weight_mton, 0)) %>%
  group_by(year, fleet, Name) %>%
  summarise(mt = sum(weight_mton)) %>% # sum up across boxes
  ungroup() %>%
  filter(year > (max(year)-5)) %>%
  group_by(fleet, Name) %>%
  summarize(mt = mean(mt)) %>%
  ungroup()

# add missing combinations
dummy_df <- expand.grid("fleet" = unique(catch_data$fleet), "Name" = unique(catch_data$Name))
catch_data <- catch_data %>%
  full_join(dummy_df) %>%
  mutate(mt = replace_na(mt, 0))

# get relative catch
catch_data <- catch_data %>%
  group_by(fleet) %>%
  mutate(tot_mt = sum(mt, na.rm = T)) %>%
  ungroup() %>%
  mutate(prop = mt / tot_mt)

# deal with codes for the fleets being different
catch_data <- catch_data %>% 
  rowwise() %>%
  mutate(fleet = rewrite_codes(as.character(fleet))) %>%
  ungroup()

# stack
catch_data <- catch_data %>% mutate(Type = "data")
catch_atlantis_fleets_end <- catch_atlantis_fleets_end %>% mutate(Type = "model")
catch_diff <- rbind(catch_data, catch_atlantis_fleets_end)
catch_diff <- catch_diff %>%
  left_join(fleet_key %>% select(Code, Name), by = c("fleet" = "Code"))

# make a plot per species
for(i in 1:length(to_keep)){
  this_fleet <- to_keep[i]
  
  print(paste("Doing", this_fleet))
  
  p <- catch_diff %>%
    filter(fleet == this_fleet) %>%
    filter(Name.x %in% all_fg) %>% # do vertebrates only
    filter(prop > 0) %>%
    ggplot(aes(x = Type, y = prop, fill = Name.x))+
    geom_bar(stat = "identity", position = "stack")+
    theme_bw()+
    facet_wrap(~Name.y)
  
  ggsave(paste0("fleets/data_vs_catch/by_fleet/", this_fleet, ".png"), p, width = 8, height = 4.5)
  
  
}

# Most fleets look somewhat similar to the data in their catch compositions, but there are some clear problems
# of availability.
# A great example is pollock. 99% of it is caught by a trawl fleet (same as data). 
# Looking at this trawl fleet however, we see that in the model it catches much more arrowtooth (or less pollock) than it should
# So, is this fleet catching more arrowtooth or less pollock? 

# There is a need to calibrate this. A few ways to do it:
# 1. Manipulate mFC. This gives you control over the species and the fleet, but not over space
# 2. Manipulate MPAYYY (it takes entries >1). This gives you control over a fleet and the areas, but not the species

# When the spatial element of the setup is what causes it to underharvest, it seems sensible to use approach 2.
# Problem with that is - you amp up a fleet and you may end up overharvesting some of the species in the fleet to get the main target right

# Perhaps step 1 should be: at model-level, can we catch the same amount of a species that we get in the non-spatial setup?
# Then we can look at species and fleet makeup and compare it to the data.

# Should we change the proportions of the mFC vector? This would make us drift away from the fleet makeup of the catch of a species, which now is on point
# So perhaps you shouls start at species level: how far off is the final catch of ATF between the two setup?
# Apply that scalar to mFC - ACROSS all entries. This way, proportions of catch among fleets will be the same(ish?), but the total catch will be higher.
# Spatial distribution of catch should not change much (or should it?)
# Species makeup of each fleet WILL change (for the better or for the worse?)


# Compare total catch by species between runs -----------------------------
catch_nc_base <- build_catch_output_v2(catch_nc = catch_nc_file_base, 
                                       fleet_struc = F,
                                       relative = T,
                                       run = 1517,
                                       key = fleet_key)

catch_nc_mpa <- build_catch_output_v2(catch_nc = catch_nc_file_mpa, 
                                      fleet_struc = F,
                                      relative = T,
                                      run = 1562,
                                      key = fleet_key)

# get residuals (it gets hazy for proportions - what's a big residual and how do you translate that to catch?)
catch_diff <- catch_nc_base %>%
  left_join(catch_nc_mpa, by = c("ts", "box_id", "Name")) %>%
  select(ts, Name, mt_goa.x, mt_goa.y) %>%
  distinct() %>%
  rename(base = mt_goa.x, mpa = mt_goa.y) 

# for plotting
catch_diff_long <- catch_diff %>%
  pivot_longer(-c(ts,Name), names_to = "run", values_to = "mt")

# view
catch_diff_long %>%
  filter(Name %in% to_plot) %>%
  ggplot()+
  geom_line(aes(x = ts, y = mt, color = run), linewidth = 1.5)+
  facet_wrap(~Name, scales = "free")

# One key question here is: How short a run can we get away with?
# The plot above indicates that, while there can still be some divergence at the end of the run,
# for most groups it is fairly apparent already at the start of the run that the fishery can't hit the non-spatial quota
# So in the interest of ironing out the big kinks fast, let's do some bried (1-2 years) runs

# for ease of visualization, add a plot for the first time step only (you will need it when you have 1-yr long runs)
catch_diff_long %>%
  filter(Name %in% all_fg) %>%
  filter(ts == 1) %>%
  ggplot()+
  geom_bar(aes(x = run, y = mt), stat = "identity")+
  facet_wrap(~Name, scales = "free")

# Get mFC scalars and apply them to the harvest.nc file -------------------

# So, based on the difference in catch between base and mpa setup, work out scalars to crank up mFC in the mpa setup
# Do so at species level, but know that this will have implications for the species makeup of a fleet's catch
# Also be mindful that: this may lead to grossly incorrect spatial patterns; and it may have no effect if a fishery fails in a box because its target isn't there
# Save these scalars because you will need to apply them when you change mFC (e.g. for Ftarg, FMSY, and others)
# use t=1 as first coarse pass. There will be remaining divergences 
scalars <- catch_diff %>%
  filter(ts == 1) %>%
  mutate(scalar = base / mpa) %>%
  select(Name, scalar) %>%
  left_join(grps %>% select(Name, Code))

# which harvest.prm file do you want to work on
harvest_old <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1559/GOA_harvest_fleets_mpa_v2.prm"
harvest_new <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1559/GOA_harvest_fleets_mpa_v3.prm"

harvest_tab <- readLines(harvest_old)

for(i in 1:nrow(scalars)){
  
  this_code <- scalars[i,]$Code
  this_scalar <- scalars[i,]$scalar
  
  # old mFC
  mfc_old <- harvest_tab[grep(paste0("mFC_", this_code, " 33"), harvest_tab)+1]
  mfc_old_vec <- as.numeric(unlist(strsplit(mfc_old, " ")))
  
  # new vec
  mfc_new_vec <- mfc_old_vec * this_scalar
  mfc_new <- paste(as.character(mfc_new_vec), collapse = " ")
  
  # replace relevant line
  harvest_tab[grep(paste0("mFC_", this_code, " 33"), harvest_tab)+1] <- mfc_new
}

# write out
writeLines(harvest_tab, con = harvest_new)
