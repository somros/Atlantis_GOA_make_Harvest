# Alberto Rovellini
# 08/05/2024
# Make a box-port matching key, by species, based on the catch reconstruction from Adam
# Desired end product:
# a matching key that maps proportion of catch in any given box to the ports in this list
# by year, by species
# then, we can average over the last 10 years, but not before looking at some time trends (what melissa mentioned)

pacman::p_load(tidyverse, sf, rbgm, viridis, here, RColorBrewer)

select <- dplyr::select
dir <- here()

# Read data ---------------------------------------------------------------

# Atlantis groups file
grps <- read.csv(here(dir, "data/GOA_Groups.csv"))
all_fg <- grps %>% filter(GroupType %in% c("FISH", "SHARK")) %>% pull(Name)

# geometry
goa_bgm <- read_bgm(here(dir, "data/GOA_WGS84_V4_final.bgm"))
goa_sf <- goa_bgm %>% box_sf()

# catch reconstruction from data
catch_dat <- readRDS(here(dir, "fleets/fleet_total_catch_atl.RDS"))
fleet_key <- read.csv(here(dir, "data/GOA_fisheries.csv"))
port_names <- read.csv("data/port_codes_and_names.csv")

# load functions
source("catch_compare_functions.R")

# View catch by port ------------------------------------------------------

# handle the catch reconstruction
catch_data <- prepare_catch_data(catch_dat = catch_dat, do_boundary = TRUE, do_islands = TRUE)

# keep only 2011-2020
catch_data <- catch_data %>% filter(year > 2010, year <= 2020)

# make key by species and year. Aggregate across fleets
catch_data <- catch_data %>%
  group_by(year, box_id, spp, port) %>% # aggregate across fleets
  summarise(mt_box = sum(weight_mton, na.rm = T))

# drop box 95, an artefact of the cleaning but empty and unwanted here
catch_data <- catch_data %>% filter(box_id < 95)

# to start, visualize this a bit
catch_data %>%
  group_by(year, port, spp) %>%
  summarise(mt = sum(mt_box)) %>%
  filter(spp == "POL") %>%
  ggplot()+
  geom_bar(aes(x = factor(year), y = mt, fill = port), stat = "identity", position = "dodge")

# moving on
# for this to work, we need to expand the data frame with all possible combinations (it will be large)
# that way we avoid NAs and ensure that averaging goes smooth at all steps
catch_data_template <- expand.grid(
  year = sort(unique(catch_data$year)),
  box_id = sort(unique(catch_data$box_id)),
  spp = sort(unique(catch_data$spp)),
  port = sort(unique(catch_data$port))
)

# merge in data
catch_data_full <- merge(catch_data_template, catch_data, by = c("year","box_id","spp","port"), all.x = T)

# now replace NA's with 0's
catch_data_full$mt_box[is.na(catch_data_full$mt_box)] <- 0

# calculate the proportion going to each port
catch_data_port <-  catch_data_full %>%
  group_by(year, box_id, spp) %>%
  mutate(tot_mt_box = sum(mt_box, na.rm = T)) %>% # total catch by box and species across ports
  ungroup() %>%
  mutate(prop = mt_box / tot_mt_box) # proportion of catch of this species in this port going to each port

# now average across years
box_port_key <- catch_data_port %>%
  group_by(box_id, spp, port) %>%
  summarise(mean_prop = mean(prop, na.rm = T),
            mean_mt = mean(mt_box, na.rm = T))

# replace many NaN's with 0 (i.e. no catch going to that port)
box_port_key$mean_prop[is.nan(box_port_key$mean_prop)] <- 0

# do these add up to one?
# box_port_key %>%
#   group_by(box_id, spp) %>%
#   summarise(check = sum(mean_prop)) %>%
#   pull(check) %>%
#   summary()
# yes

# add port long names
box_port_key <- box_port_key %>% 
  left_join(port_names, by = c('port'='Port.Code'))

# prepare for writing out
box_port_key <- box_port_key %>%
  select(box_id, spp, port, Port.Name, mean_prop, mean_mt) %>%
  rename(Code = spp,
         Port.Code = port)

saveRDS(box_port_key, file = "data/box_to_port_key.RDS")

# view some examples
# composition
# pollock
box_port_key %>%
  filter(Code == "POL") %>%
  filter(mean_prop > 0) %>%
  ggplot()+
  geom_bar(aes(x = box_id, y = mean_prop, fill = Port.Name), stat = "identity", position = "stack")

# cod
box_port_key %>%
  filter(Code == "COD") %>%
  filter(mean_prop > 0) %>%
  ggplot()+
  geom_bar(aes(x = box_id, y = mean_prop, fill = Port.Name), stat = "identity", position = "stack")

# we also need proportion of the catch by port
# this will be used to evaluate Atlantis skill
# here the catch by port is from the data, in the run output it will be based off this key
# total catch by port
catch_by_port <- box_port_key %>%
  group_by(Code, Port.Code, Port.Name) %>%
  summarize(mean_mt_port = sum(mean_mt, na.rm = T)) %>% # drop boxes
  group_by(Code) %>%
  mutate(mt_tot = sum(mean_mt_port, na.rm = T)) %>%
  ungroup() %>%
  mutate(prop = mean_mt_port / mt_tot)

catch_by_port$prop[is.nan(catch_by_port$prop)] <- 0 # NaN to 0

# view
catch_by_port %>%
  filter(Code == "POL", prop > 0.01) %>%
  mutate(Port.Name = fct_reorder(Port.Name, prop, .desc = T)) %>%
  ggplot()+
  geom_bar(aes(x = Port.Name, y = prop), stat = "identity", position = "stack")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Apportion catch ---------------------------------------------------------
# now let's take a test run and apportion it based on the key
# this will need to become a function that writes out plots
# use the CATCHTOT.nc file here (no txt file exists with spatial catch output)

# box_port_key <- readRDS("data/box_to_port_key.RDS")
# 
# run <- 1574 
# data_dir <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/" # directory with the Atlantis runs
# 
# # nc file
# catch_nc_file_TOT <- paste0(data_dir, "/out_", run, "/outputGOA0", run, "_testTOTCATCH.nc") # full nc catch file
# 
# catch_tot <- build_catch_output_TOT(catch_nc_file_TOT, run = run)
# 
# # drop BC
# catch_tot <- catch_tot %>% filter(box_id < 92)
# 
# # tie in the port key and perform operations
# catch_mapped <- full_join(catch_tot, key %>% select(-mean_mt)) 
# 
# # there are several NA's that appear from these joint. They include:
# # 1. Groups not in the data reconstruction (FOS)
# # 2. Boundary boxes, such as box 2
# # 3. Invertebrates
# # Need to look a bit deeper into this, for now scratch these
# catch_mapped <- catch_mapped[complete.cases(catch_mapped),]
# 
# # tie in the port key and perform operations
# catch_mapped <- catch_mapped %>%
#   mutate(mt_new = mt * mean_prop) %>% # get catch to each port by box and species
#   group_by(ts, Name, Code, Port.Name, Port.Code) %>%
#   summarise(mt_tot_port = sum(mt_new)) %>% # drop boxes
#   group_by(ts, Name, Code) %>%
#   mutate(mt_tot = sum(mt_tot_port, na.rm =T)) %>%
#   ungroup() %>%
#   mutate(prop = mt_tot_port / mt_tot)
# 
# # view
# catch_mapped %>%
#   slice_max(ts) %>%
#   filter(Code == "POL", prop > 0.001) %>%
#   mutate(Port.Name = fct_reorder(Port.Name, prop, .desc = T)) %>%
#   ggplot()+
#   geom_bar(aes(x = Port.Name, y = prop), stat = "identity", position = "stack")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# # set up comparison plot
# tt <- catch_mapped %>%
#   rename(prop_out = prop) %>%
#   select(-mt_tot_port,-mt_tot) %>%
#   left_join(catch_by_port %>%
#               rename(prop_in = prop) %>%
#               select(-mean_mt_port, -mt_tot), 
#             by = c('Code','Port.Name','Port.Code')) %>%
#   pivot_longer(c(prop_in, prop_out), names_to = "type", values_to = "prop")
# 
# tt %>%
#   slice_max(ts) %>%
#   filter(Code == "POL", prop > 0.001) %>%
#   mutate(Port.Name = fct_reorder(Port.Name, prop, .desc = T)) %>%
#   ggplot()+
#   geom_bar(aes(x = Port.Name, y = prop, fill = type), stat = "identity", position = "dodge")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# 
# # Now put these into a function so that it is done for all species in a run
