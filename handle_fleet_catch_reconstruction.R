# This scrips is called to perform some operations on the catch reconstruction from Adam
# Main tasks are:
# Turns catch in boundary boxes to 0
# Reallocates catch from island boxes to neighboring boxes based on the proportion of the catch in said naighbors

# catch reconstruction
fleets <- readRDS("fleets/fleet_total_catch_atl.RDS")
fleet_key <- read.csv("data/GOA_fisheries.csv")

# handle islands and boundary boxes
# zero-out BBs (TODO: thinks if this makes sense)
bboxes <- goa_sf %>% filter(boundary == TRUE) %>% pull(box_id)
fleets <- fleets %>%
  rowwise() %>%
  mutate(weight_mton = ifelse(box_id %in% bboxes,0,weight_mton)) %>%
  ungroup()

# handle islands: attribute catch from an island box to its neighbors based on the existing catch props
# 1. identify island boxes
# 2. identify their neighbors
# 3. for each island box, year, fleet, species, port, get the catch value from the neighboring boxes
# 4. Turn said catch values into proportions
# 5. split catch from the island box into chunks based on calculated proportions
# 6. Add the new catch to the respective neighbors
# 7. set catch in all island boxes to 0

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
  group_by(year, island_box, fleet, spp, port) %>%
  mutate(total_weight = sum(weight_mton, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(prop = weight_mton / total_weight) %>%
  select(year, island_box, fleet, spp, port, box_id, prop) %>%
  rename(neighbors = box_id) %>%
  rename(box_id = island_box)

# check that these all add up to 1
# fleets_islands_prop %>%
#   group_by(year, island_box, fleet, spp, port) %>%
#   summarise(check = sum(prop)) %>%
#   pull(check) %>%
#   summary() # OK

# 5. split catch from the island box into chunks based on calculated proportions
fleets_islands_reallocated <- fleets_islands_prop %>%
  left_join(fleets, by = c("year", "box_id", "fleet", "spp", "port")) %>%
  mutate(weight_mton_2 = weight_mton * prop) %>% # this is the new weight to add to the neighbors
  drop_na()

# now let's do some checks. The total catch in column weight_mton_2 summed by year needs to be the same as the raw catch
# if it is not, we lost some of it on the way (for instance what if a metier/species.port combo fished in the island box but nowhere else? Then it would not appear here)
# Expect differences
# check1 <- fleets %>%
#   filter(box_id %in% island_boxes) %>%
#   group_by(year, fleet) %>%
#   summarise(tot = sum(weight_mton))
# 
# check2 <- fleets_islands_reallocated %>%
#   group_by(year, fleet) %>%
#   summarise(tot = sum(weight_mton_2))
# 
# check3 <- left_join(check1, check2, by = c("year", "fleet")) %>%
#   mutate(lost = tot.x - tot.y,
#          lost_percent = (tot.x - tot.y) / tot.x * 100)

# we lose some catch this way 
# the worst cases are: 
# 75% of CG_OTH_SMX_ROA in 1998
# 65% of CG_OTH_SMX_ROA in 1997
# 21% of CG_OTH_SSO_KI in 2019
# 20% of CG_OTH_SMX_ROA in 1995
# 17% of CG_OTH_SSO_KI in 2021
# 4.4% of CG_OTH_SMX_ROA in 2006

# 6. Add the new catch to the respective neighbors
fleets_islands_reallocated <- fleets_islands_reallocated %>%
  select(year, neighbors, fleet, spp, port, weight_mton_2)

fleets_2 <- fleets %>%
  left_join(fleets_islands_reallocated, by = c("year","box_id" = "neighbors", "fleet", "spp", "port")) %>%
  mutate(weight_mton_2 = replace_na(weight_mton_2, 0)) %>% # new weight to add is 0 if you are not the neighbor of an island
  mutate(weight_mton_3 = weight_mton + weight_mton_2) %>%
  select(-weight_mton, -weight_mton_2) %>%
  rename(weight_mton = weight_mton_3)

# 7. set catch in all island boxes to 0
fleets_2 <- fleets_2 %>%
  rowwise() %>%
  mutate(weight_mton = ifelse(box_id %in% island_boxes, 0, weight_mton)) %>%
  ungroup()

# now check
# check4 <- fleets %>%
#   group_by(year, fleet) %>%
#   summarise(tot  = sum(weight_mton))
# 
# check5 <- fleets_2 %>%
#   group_by(year, fleet) %>%
#   summarise(tot  = sum(weight_mton))
# 
# check6 <- left_join(check4,check5, by = c("year","fleet")) %>%
#   mutate(lost = tot.x - tot.y,
#          lost_percent = (tot.x - tot.y) / tot.x * 100)
# small amount lost, revisit this in the future

# replace the original fleet
fleets <- fleets_2