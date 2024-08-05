# Alberto Rovellini
# 04/24/2024
# This code prepares the MPAYYY vectors for the fishing fleets based on Adam's catch reconstruction
# the goal is, for each fleet:
# MPAYYY 109
# m0 m1 m2 m3 ... m108
# options: open the whole box or open a proportion of it
# Use a broad enough period of time that we do not limit the fishery to a small area where it operated in recent years
# So, use at least the last 10 years
# we also need to: set catches in the boudary boxes to 0, and to allocate catch from Island boxes to nearby boxes
# Do so based on the proportion of catch of that fleet in the nearby boxes

library(tidyverse)
library(rbgm)
library(sf)
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

# at this stage, drop species from the data frame - we only care about total catch by fleet
# doing this now should help us avoid some issues with reallocating catch from island boxes
fleets <- fleets %>%
  group_by(year, fleet, box_id) %>%
  summarise(weight_mton = sum(weight_mton, na.rm = T))

# some catch bleeds into box 92, which is in BC. Set that to 0
fleets <- fleets %>%
  rowwise() %>%
  mutate(weight_mton = ifelse(box_id >= 92, 0, weight_mton)) %>%
  ungroup()

# handle boundary boxes
# Here I assume that Adam was not made aware that we do not include those and mapped actual catch from stat areas to the respective boxes
# This as opposed to catch in the BB being an artefact of spillover from adjacent boxes and inaccurate spatial mismatch
# Note that this can reduce the catches of some species substantially (e.g. sablefish), which has implications for things like calculating F
fleets <- fleets %>%
  mutate(weight_mton = ifelse(box_id %in% boundary_boxes, 0, weight_mton))

# handle island boxes
# allocate catch based on proportion of the catch from that fleet in each neighbor
# 1. identify island boxes
island_boxes <- goa_sf %>% filter(botz == 0) %>% pull(box_id) # 21, 44, 99 (Shumagin, Kodiak, Haida Gwaii)

# 2. identify their neighbors
neighbors <- st_touches(goa_sf) # identify THE INDEX of the touching geometries
neighbors <- neighbors[island_boxes+1] # keep islands only

# Initialize an empty data frame
neighbors_key <- data.frame(island_box = integer(0), neighbor = integer(0))

# Loop through the list and create the data frame
for (i in seq_along(island_boxes)) {
  this_island_box <- island_boxes[i]
  these_neighbors <- neighbors[[i]] - 1 # turn from INDEX to box_id
  
  # Create a temporary data frame for the current focal geometry
  temp_df <- data.frame(island_box = this_island_box, neighbor = these_neighbors)
  
  # Append the temporary data frame to the result data frame
  neighbors_key <- rbind(neighbors_key, temp_df)
}

# nest this
neighbors_key <- nest(neighbors_key, .by = neighbor, .key = "island_box")

# 3. for each island box, year, fleet, species, port, get the catch value from the neighboring boxes
# First subset the data to the boxes of interest - neighbors and their anchor islands
# Then, for each island box and by accounting for grouping vars, look at the split of the catch between the neighbors
fleets_islands <- fleets %>%
  left_join(neighbors_key, by = c("box_id" = "neighbor")) %>%
  unnest(cols = c(island_box)) 

# 4. Turn said catch values into proportions
fleets_islands_prop <- fleets_islands %>%
  group_by(year, island_box, fleet) %>%
  mutate(total_weight = sum(weight_mton, na.rm = TRUE)) %>% # total catch among the neighbors by year by fleet
  ungroup() %>%
  mutate(prop = weight_mton / total_weight) %>%
  select(year, island_box, fleet, box_id, prop) %>%
  rename(neighbors = box_id) %>%
  rename(box_id = island_box)

# check that these all add up to 1
# fleets_islands_prop %>%
#   group_by(year, box_id, fleet) %>%
#   summarise(check = sum(prop)) %>%
#   pull(check) %>%
#   summary() # OK

# 5. split catch from the island box into chunks based on calculated proportions
fleets_islands_reallocated <- fleets_islands_prop %>%
  left_join(fleets, by = c("year", "box_id", "fleet")) %>%
  mutate(weight_mton_2 = weight_mton * prop) %>% # this is the new weight to add to the neighbors
  drop_na()

# 6. Add the new catch to the respective neighbors
fleets_islands_reallocated <- fleets_islands_reallocated %>%
  select(year, neighbors, fleet, weight_mton_2)

fleets_2 <- fleets %>%
  left_join(fleets_islands_reallocated, by = c("year","box_id" = "neighbors", "fleet")) %>%
  mutate(weight_mton_2 = replace_na(weight_mton_2, 0)) %>% # new weight to add is 0 if you are not the neighbor of an island
  mutate(weight_mton_3 = weight_mton + weight_mton_2) %>%
  select(-weight_mton, -weight_mton_2) %>%
  rename(weight_mton = weight_mton_3)

# 7. set catch in all island boxes to 0
fleets_2 <- fleets_2 %>%
  rowwise() %>%
  mutate(weight_mton = ifelse(box_id %in% island_boxes, 0, weight_mton)) %>%
  ungroup()

# average across years - i.e. only keep the spatial element
fleets_spatial <- fleets_2 %>%
  left_join(fleet_key, by = c("fleet"="Fleet")) %>%
  filter(year > 2010) %>% # last 10 years
  group_by(fleet, Primary.Spp, Gear, Fishing.Area, Landing.Area, box_id) %>% # mean across years
  summarize(catch_mt = mean(weight_mton)) %>%
  ungroup()

# pad the data frame with missing combinatiosn of fleets and boxes so that they show up as 0's
all_fleets <- unique(fleets_spatial$fleet)
# first, expand to all possible combinations fleets / box
combinations <- expand.grid("fleet" = all_fleets, "box_id" = 0:108) %>%
  left_join(fleet_key, by = c("fleet"="Fleet"))

fleets_spatial <- combinations %>%
  full_join(fleets_spatial)

# turn NA's to 0
fleets_spatial$catch_mt[is.na(fleets_spatial$catch_mt)] <- 0

# view this:
# make it spatial
fleets_spatial_plot <- goa_sf %>%
  select(box_id) %>%
  left_join(fleets_spatial, by = "box_id") %>%
  mutate(catch_mt = ifelse(catch_mt == 0, NA, catch_mt))

for(i in 1:length(all_fleets)){

  # pick out a fleet
  this_fleet_name <- all_fleets[i]

  # filter df to that fleet
  this_fleet <- fleets_spatial_plot %>%
    filter(fleet == this_fleet_name)

  # define gear, target, fishing area, landing area
  this_target <- this_fleet %>% pull(Primary.Spp) %>% unique()
  this_gear <- this_fleet %>% pull(Gear) %>% unique()
  this_fishing_area <- this_fleet %>% pull(Fishing.Area) %>% unique()
  this_landing_area <- this_fleet %>% pull(Landing.Area) %>% unique()

  map_by_fleet <-  ggplot()+
    geom_sf(data = this_fleet, aes(fill = catch_mt))+
    scale_fill_viridis()+
    theme_bw()+
    labs(title = paste0("Target = ", this_target, "\n",
                        "Gear = ", this_gear, "\n",
                        "Fishing area = ", this_fishing_area, "\n",
                        "Landing area = ", this_landing_area))

  # save plot
  ggsave(paste0("fleets/maps_by_fleet_clean/", this_fleet_name, ".png"), map_by_fleet, width = 8, height = 4)

}

# Produce MPAYY vectors ---------------------------------------------------
# start from a prm we have
file_path <- "data/"
bg_harvest_file <- paste0(file_path, "GOA_harvest_fleets_v2_10YR.prm") # mFC is apportioned to fleets but MPAYYY are still undefined
goa_fleets_file <- paste0(file_path, "GOA_fisheries.csv")

# read in
bg_harvest <- readLines(bg_harvest_file)
goa_fisheries <- read.csv(goa_fleets_file)

# fix typo in CSV
goa_fisheries$Code <- gsub("CgOthSpiKo", "CgOthSpiKi", goa_fisheries$Code)

# rewrite fleet codes to match format
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

# rewrite codes
fleets_spatial <- fleets_spatial %>% 
  rowwise() %>%
  mutate(fleet_code = rewrite_codes(as.character(fleet))) %>%
  ungroup()

# define actual fleet codes
fleet_codes <- c(unique(fleets_spatial$fleet_code), "Canada")

# define what level of access the fleet has to each box:
# access_case:
# 1: if the fleet has fished there in the reference period, the box is fully open (1)
# 2: assign "openness" to boxes based on the proportion of the fleet's catch coming from that box
# 3: determines that the box with the highest CPUE is fully open, and the other boxes are open based on the proportion of their CPUE compared to the highest (most arbitrary)
# In either case, the realized catch out of that box will be a function of not only MPAYYY but also the spatial distributions of the fished species
# This means that every case will need calibration of mFC

access_case <- 2

for(i in 1:length(fleet_codes)){
  this_f <- fleet_codes[i]
  this_df <- fleets_spatial %>% filter(fleet_code == this_f)
  idx <- grep(paste0('MPA', this_f, ' 109'), bg_harvest)
  
  if(access_case == 1) {
    
    # turn catches to 0-1
    this_vec <- this_df %>% 
      rowwise() %>%
      mutate(is_open = ifelse(catch_mt == 0, 0, 1)) %>% 
      ungroup() %>%
      pull(is_open) %>%
      as.character() %>%
      paste(collapse = " ")
    
  } else if (access_case == 2) {
    
    this_vec <- this_df %>% 
      mutate(tot_catch = sum(catch_mt, na.rm = T),
             prop = catch_mt / tot_catch) %>% 
      pull(prop) %>%
      format(digits = 7) %>%
      as.character() %>%
      paste(collapse = " ")
    
  } else if (access_case == 3){
    
    this_vec <- this_df %>%
      left_join(goa_sf %>% # bring in areas
                  st_set_geometry(NULL) %>%
                  select(box_id,area),
                by = "box_id") %>%
      mutate(cpue = catch_mt / area, # get cpue by box
             max_cpue = max(cpue, na.rm = F), # get max cpue
             rel_cpue = cpue / max_cpue) %>% # scale other cpues to max value
      pull(rel_cpue) %>%
      format(digits = 7) %>%
      as.character() %>%
      paste(collapse = " ")
    
  } else {
    
    stop("This access_level does not exist")
    
  }
  
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

# general parameter flagmpa
idx <- grep('\\bflagmpa', bg_harvest)
lab <- gsub("^.*?(##)", "\\1", bg_harvest[idx])
bg_harvest[idx] <- paste0("flagmpa 1", " ", lab)

# flagYYYday to 2 so that the fishery operates day and night
for(i in 1:length(fleet_codes)){
  
  idx <- grep(paste0('flag', fleet_codes[i], "day"), bg_harvest)
  bg_harvest[idx] <-  gsub(" 0$", " 2", bg_harvest[idx])
  
}

# write out
writeLines(bg_harvest, con = "data/GOA_harvest_fleets_mpa_CASE2_10YR.prm")

# Make a dummy MPA harvest file where the MPA vectors are all 1 -----------
# This has the purpose of testing if the MPA machinery works in Atlantis

# # read in
# bg_harvest <- readLines("data/GOA_harvest_fleets_mpa_v2.prm")
# 
# for(i in 1:length(fleet_codes)){
#   this_f <- fleet_codes[i]
#   idx <- grep(paste0('MPA', this_f, ' 109'), bg_harvest)
#   
#   # make a vector 
#   if(this_f != "Canada"){
#     this_vec <- c(rep(1,92), rep(0, (109-92)))
#   } else {
#     this_vec <- c(rep(0,92), rep(1, (109-92)))
#   }
#   
#   this_vec <- this_vec %>%
#     as.character() %>%
#     paste(collapse = " ")
#   
#   bg_harvest[idx + 1] <- this_vec
#   
# }
# 
# writeLines(bg_harvest, con = "data/GOA_harvest_fleets_mpa_check_v2.prm")
