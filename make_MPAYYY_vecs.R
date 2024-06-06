# Alberto Rovellini
# 04/24/2024
# This code prepares the MPAYYY vectors for the fishing fleets based on Adam's catch reconstruction
# the goal is, for each fleet:
# MPAYYY 109
# m0 m1 m2 m3 ... m108
# as a first pass, let's open fully those boxes where some fishing has occurred
# we also need to: set catches in the boudary boxes to 0, and to allocate catch from Island boxes to nearby boxes
# let's use the last 5 years of data for each fleet for the spatial coverage. Why? because we want to project this.

library(tidyverse)
library(rbgm)
library(sf)
library(viridis)

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
island_boxes <- goa_sf %>% filter(botz == 0) %>% pull(box_id)

fleets_spatial <- fleets %>%
  left_join(fleet_key, by = c("fleet"="Fleet")) %>%
  filter(year > 2015) %>% # last 5 years
  group_by(year, fleet, Primary.Spp, Gear, Fishing.Area, Landing.Area, box_id) %>%
  summarize(catch_mt = sum(weight_mton)) %>%
  group_by(fleet, Primary.Spp, Gear, Fishing.Area, Landing.Area, box_id) %>%
  summarize(catch_mt = mean(catch_mt)) %>%
  ungroup()

all_fleets <- unique(fleets_spatial$fleet)
# first, expand to all possible combinations fleets / box
combinations <- expand.grid("fleet" = all_fleets, "box_id" = 0:108) %>%
  left_join(fleet_key, by = c("fleet"="Fleet"))

fleets_spatial <- combinations %>%
  full_join(fleets_spatial)

# turn NA's to 0
fleets_spatial$catch_mt[is.na(fleets_spatial$catch_mt)] <- 0

# handle boundary boxes
# Here I assume that Adam was not made aware that we do not include those and mapped actual catch from stat areas to the respective boxes
# This as opposed to catch in the BB being an artefact of spillover from adjacent boxes and inaccurate spatial mismatch
# Note that this can reduce the catched of some species substantially (e.g. sablefish), which has implications for things like calculating F
fleets_spatial <- fleets_spatial %>%
  mutate(catch_mt = ifelse(box_id %in% boundary_boxes, 0, catch_mt))

# handle island boxes
# do it as a loop - for each fleet for each island box
# make a key first
n_ls <- list()
for(i in 1:length(island_boxes)){
  # get neighboring boxes
  neighbors <- (goa_sf %>% filter(box_id == island_boxes[i]) %>% st_touches(goa_sf) %>% unlist()) - 1
  # get the areas of the neighbors as weights
  neighbor_area_props <- goa_sf %>% filter(box_id %in% neighbors) %>% st_area() %>% as.numeric() %>% prop.table()
  # bind
  n_ls[[i]] <- data.frame("box_id" = island_boxes[i],
                          "neighbor" = neighbors,
                          "props" = neighbor_area_props)
}
n_df <- bind_rows(n_ls)

# now loop through fleets and then through island boxes and then through neighbors
# for now, which species is getting caught does not matter, it only matters that the fleet is operating here
for(f in 1:length(all_fleets)){
  this_f <- all_fleets[f]
  
  for (b in 1:length(island_boxes)){
    this_b <- island_boxes[b]
    this_key <- n_df %>% filter(box_id == this_b)
    neighbors <- this_key %>% pull(neighbor)
    neighbor_ap <- this_key %>% pull(props)
    this_catch <- fleets_spatial %>% 
      filter(fleet == this_f, box_id == this_b) %>%
      pull(catch_mt)
    
    # allocate it to the neighbors
    for(n in 1:length(neighbors)){
      this_n <- neighbors[n]
      this_p <- neighbor_ap[n]
      
      this_n_catch <- fleets_spatial[fleets_spatial$fleet == this_f & fleets_spatial$box_id == this_n,]$catch_mt
      
      # add catch from island
      fleets_spatial[fleets_spatial$fleet == this_f & fleets_spatial$box_id == this_n,]$catch_mt <- this_n_catch + (this_catch * this_p)
      
    }
    
    # set it to 0 in the island box itself
    fleets_spatial[fleets_spatial$fleet == this_f & fleets_spatial$box_id == this_b,]$catch_mt <- 0
    
  }
  
}

# view this:
# make it spatial
# fleets_spatial_plot <- goa_sf %>%
#   select(box_id) %>%
#   left_join(fleets_spatial, by = "box_id") %>%
#   mutate(catch_mt = ifelse(catch_mt == 0, NA, catch_mt))
# 
# for(i in 1:length(all_fleets)){
#   
#   # pick out a fleet
#   this_fleet_name <- all_fleets[i]
#   
#   # filter df to that fleet
#   this_fleet <- fleets_spatial_plot %>%
#     filter(fleet == this_fleet_name)
#   
#   # define gear, target, fishing area, landing area
#   this_target <- this_fleet %>% pull(Primary.Spp) %>% unique()
#   this_gear <- this_fleet %>% pull(Gear) %>% unique()
#   this_fishing_area <- this_fleet %>% pull(Fishing.Area) %>% unique()
#   this_landing_area <- this_fleet %>% pull(Landing.Area) %>% unique()
#   
#   map_by_fleet <-  ggplot()+
#     geom_sf(data = this_fleet, aes(fill = catch_mt))+
#     scale_fill_viridis()+
#     theme_bw()+
#     labs(title = paste0("Target = ", this_target, "\n",
#                         "Gear = ", this_gear, "\n",
#                         "Fishing area = ", this_fishing_area, "\n",
#                         "Landing area = ", this_landing_area))
#   
#   # save plot
#   ggsave(paste0("fleets/maps_by_fleet_clean/", this_fleet_name, ".png"), map_by_fleet, width = 8, height = 4)
#   
# }

# now produce vectors
# start from a prm we have
file_path <- "data/"
bg_harvest_file <- paste0(file_path, "GOA_harvest_fleets_v2.prm")
goa_fleets_file <- paste0(file_path, "GOA_fisheries.csv")

# read in
bg_harvest <- readLines(bg_harvest_file)
goa_fisheries <- read.csv(goa_fleets_file)

# fix typo
goa_fisheries$Code <- gsub("CgOthSpiKo", "CgOthSpiKi", goa_fisheries$Code)

# get codes, and rewrite fleet codes to match format
fleets_prop_term <- read.csv("data/mFC_prop_by_fleet.csv")
fg_codes <- unique(fleets_prop_term$spp)
fleet_codes_csv <- goa_fisheries %>% pull(Code)
fleet_codes <- fleets_spatial %>% select(fleet) %>% distinct()

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

fleets_spatial <- fleets_spatial %>% 
  rowwise() %>%
  mutate(fleet_code = rewrite_codes(as.character(fleet))) %>%
  ungroup()

fleet_codes <- unique(fleets_spatial$fleet_code)
for(i in 1:length(fleet_codes)){
  this_f <- fleet_codes[i]
  this_df <- fleets_spatial %>% filter(fleet_code == this_f)
  idx <- grep(paste0('MPA', this_f, ' 109'), bg_harvest)
  
  # turn numbers to 0-1
  this_vec <- this_df %>% 
    mutate(is_open = ifelse(catch_mt == 0, 0, 1)) %>% 
    pull(is_open) %>%
    as.character() %>%
    paste(collapse = " ")
  
  bg_harvest[idx + 1] <- this_vec
  # other relevant parameters:
  idx <- grep(paste0(this_f, '_flagmpa'), bg_harvest)
  lab <- gsub("^.*?(##)", "\\1", bg_harvest[idx])
  bg_harvest[idx] <- paste0(this_f, "_flagmpa 1", " ", lab)
  
}

# make a vector for Canada
# assume that, within Canada, all boxes are open
idx <- grep('MPACanada 109', bg_harvest)

# turn numbers to 0-1
# make a vector of 1 for each box in Canada (boxes 93 and up)
this_vec <- c(rep(0,92), rep(1, (109-92))) %>%
  as.character() %>%
  paste(collapse = " ")

bg_harvest[idx + 1] <- this_vec
# other relevant parameters:
idx <- grep('Canada_flagmpa', bg_harvest)
lab <- gsub("^.*?(##)", "\\1", bg_harvest[idx])
bg_harvest[idx] <- paste0(this_f, "_flagmpa 1", " ", lab)

# write out
writeLines(bg_harvest, con = "data/GOA_harvest_fleets_mpa_v2.prm")

# Make a dummy MPA harvest file where the MPA vectors are all 1 -----------
# This has the purpose of testing if the MPA machinery works in Atlantis

# read in
bg_harvest <- readLines("data/GOA_harvest_fleets_mpa_v2.prm")

for(i in 1:length(fleet_codes)){
  this_f <- fleet_codes[i]
  idx <- grep(paste0('MPA', this_f, ' 109'), bg_harvest)
  
  # make a vector of 1 for each box
  this_vec <- rep(1, 109) %>%
    as.character() %>%
    paste(collapse = " ")
  
  bg_harvest[idx + 1] <- this_vec
  
}

writeLines(bg_harvest, con = "data/GOA_harvest_fleets_mpa_check.prm")
