# Alberto Rovellini
# Jan 5 2022
# Modified from Hem Nalini Morzaria Luna's code for Puget Sound to work on Gulf of Alaska

library(tidyverse)
library(data.table)

# read in catch data for 1990
source('make_dummy_catch.R')
all.goa.fisheries <- fread("data/all_fisheries_goa_100FMSY.csv") %>% 
  as_tibble

# this is catch by fleet
all.goa.fisheries %>% pull(year) %>% unique()
all.goa.fisheries %>% pull(gear_name) %>% unique()
all.goa.fisheries %>% pull(atlantis_fg) %>% unique()
all.goa.fisheries %>% pull(box) %>% unique()
all.goa.fisheries %>% pull(month) %>% unique()

group.list <- read_csv("data/GOA_Groups.csv") # this is the group.csv file

# this is initial biomass, taken from the biomindx.txt file from any one run for t0
goa.biomass <- read_delim("data/GOA_BiomIndx.txt", delim = " ") %>% 
  dplyr::select(Time:DR) %>% 
  gather(Code,biomass_mt, -Time) %>% 
  filter(Time==0) %>% 
  left_join(group.list, by="Code") %>% 
  dplyr::rename(atlantis_fg = Code) %>% 
  dplyr::select(atlantis_fg, biomass_mt)

# use GOA biomass by age as an option
# selex
selex <- read.csv('data/age_at_selex.csv')

goa.biomass.selex <- biom_age <- read.table('data/GOA_AgeBiomIndx.txt', sep = ' ', header = T) %>%
  filter(Time == 0) %>%
  pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass_mt') %>%
  separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
  select(-Time) %>%
  left_join(selex, by = 'Code') %>%
  mutate(idx = as.numeric(Age) - as.numeric(age_class)) %>%
  filter(is.na(idx) | idx >= 0) %>%
  group_by(Code) %>%
  summarise(biomass_mt = sum(biomass_mt)) %>%
  ungroup() %>%
  rename(atlantis_fg = Code)

# this is the list of fleets
# this includes fleets that have 0 catch in the catch time series read in above
fleet.list <- read_csv("data/GOA_fisheries.csv") %>% 
  dplyr::rename(gear_name = Name)

# read in functions
source('at_harvest_functions.R')

# make ancillary prm files used to populate harvest.prm below
source('make_ancillary_prm.R')

#this generates the fishing mortality rates
mfc.tibble <- make_mfc(all.goa.fisheries, fleet.list, group.list, goa.biomass)

#create at_harvest
grp.file  <-  "data/GOA_Groups.csv" # grp file
fsh.file  <-  "data/goa_fleet_values.csv" # same as GOA_fisheries.csv but with a couple more fields
temp  <-  'data/PrmFishTemplate.csv' # this is a spreadsheet with fishing parameters
bgm.file  <-  'data/GOA_WGS84_V4_final.bgm' # geometry
cum.depths  <-  c(0, 30, 100, 200, 500, 1000, 4000) # depths
harvest.file.name <-  "data/GOA_harvest_100FMSY.prm"
run.type  <- "future"#"historical"
this.mfc <- "data/mfc_vector.prm"

make_at_harvest(grp.file, fsh.file, temp, bgm.file, cum.depths, harvest.file.name, run.type, this.mfc)

# need to then add the parameters below for newer versions of the code
file.name <- "data/new_parameters.prm"
file.create(file.name)

# spawn_closure_
newpar <- group.list %>% 
  select(Code) %>% 
  mutate(param = paste0('spawn_closure_', Code, ' 1')) %>% 
  pull(param) 

for(i in 1:length(newpar)){
  cat(newpar[i],file=file.name,append=TRUE, '\n')
  cat(0, file=file.name,append = TRUE, '\n')
}

# _mFC_endage
newpar <- group.list %>% 
  select(Code) %>% 
  mutate(param = paste0(Code, '_mFC_endage', ' 33')) %>% 
  pull(param) 

for(i in 1:length(newpar)){
  cat(newpar[i],file=file.name,append=TRUE, '\n')
  cat(rep(10, 33), file=file.name,append = TRUE, '\n')
}
