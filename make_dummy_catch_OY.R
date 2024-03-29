# Alberto Rovellini
# 1/5/2023
# Make dummy catch file for Atlantis GOA background F based on:
# Initial biomass (1990)
# FMSY from latest stock assessment where available
# M where FMSY is not available
# Consider background F (or a starting point) 25% of FMSY or M
# This assumes that F in 2020 is comparable to F in 1990
# This code also adds catch (0) for all other fleets, to have them as placeholders

# UPDATE 2/15/2023
# Added the option to calculate catch and in turn mFC on selected age classes
# This is based on Atlantis' implementation of age selectivity, i.e. knife-edge 
# For tier 3 we put this at the age class that results in >50% selex; for all others to age at 50% maturity

# UPDATE 12/21/2023
# This was based on old selectivity calculations that have since been updated
# Recoded this to use more recent values of the knife edge 
# Simply pulling startage from the harvest.prm from 1473, which is a recent version with updated knife edges


library(tidyverse)
library(readxl)
library(data.table)
library(lubridate)

# groups
fg <- read.csv('data/GOA_Groups.csv', header = T)
fg_codes <- fg %>% pull(Code)

# fleets
fleets <- read.csv('data/GOA_fisheries.csv')
fleet_codes <- fleets %>% pull(Code)
fleet_names <- fleets %>% pull(Name)
n_fleets <- length(fleet_codes)

# initial biomass from Atlantis output
biom <- read.table('../../Parametrization/output_files/data/out_1473/outputGOA01473_testBiomIndx.txt', sep = ' ', header = T) %>%
  filter(Time == 0) %>%
  pivot_longer(-Time, names_to = 'Code', values_to = 'Biomass') %>%
  filter(Code %in% fg$Code) %>%
  select(-Time)

# initial biomass at age from Atlantis output
biom_age <- read.table('../../Parametrization/output_files/data/out_1473/outputGOA01473_testAgeBiomIndx.txt', sep = ' ', header = T) %>%
  filter(Time == 0) %>%
  pivot_longer(-Time, names_to = 'Code.Age', values_to = 'Biomass') %>%
  separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
  select(-Time)

# selectivity curves
# in this version we are taking these from the harvest.prm
# to avoid circularity (this file is needed to build the harvest.prm file), we pull these from CSV files that have been informed from the prm files
selex <- read.csv("data/age_at_selex_OY.csv", header = T) %>% rename(age_class = age_class_selex)
mat <- read.csv("data/age_at_mat_OY.csv", header = T)

# boxes
bgm  <- readLines('data/GOA_WGS84_V4_final.bgm')
n_box    <- as.numeric(gsub('nbox', '', grep("nbox", bgm, value = TRUE)))

# read in F and M
tier3 <- read_xlsx('data/GOA MSY estimates tables.xlsx', sheet = 1, range = 'A3:J19') %>%
  select(Stock, FOFL) %>%
  set_names(c('Stock', 'FMSY')) 

# fix cod with info from Pete Hulson:
tier3[tier3$Stock == "Pacific cod","FMSY"] <- 0.5100000

# add halibut (M = 0.2)
tier3 <- rbind(tier3, data.frame("Stock"="Pacific halibut","FMSY"= 0.2))

tier4_5 <- read_xlsx('data/GOA MSY estimates tables.xlsx', sheet = 2, range = 'A3:I10') %>%
  select(`Stock/Stock complex`, `M or FMSY`)%>%
  set_names(c('Stock', 'FMSY'))

tier_3_4_5 <- rbind(tier3, tier4_5)

# make key
tier_3_4_5 <- tier_3_4_5 %>%
  mutate(Code = c('POL','COD','SBF','FFS','FFS','FFS','FFS','FFD',
                  'REX','REX','ATF','FHS','POP','RFS','RFS','RFP',
                  'HAL','FFS','RFD','RFD','RFD','THO','DOG','SKB')) %>%
  group_by(Code) %>%
  summarise(FMSY = mean(FMSY)) %>%
  set_names('Code','FMSY')

# extend with SKL and SKO
tier_3_4_5 <- rbind(tier_3_4_5, data.frame("Code"=c("SKL","SKO"), "FMSY"=tier_3_4_5[tier_3_4_5$Code=="SKB",]$"FMSY"))

# get M for other groups from parameters
other_M <- read.csv('data/life_history_parameters.csv') %>%
  select(Code, M_FUNC) %>%
  filter(Code %in% setdiff(fg$Code, tier_3_4_5$Code)) %>%
  set_names('Code','FMSY')

all <- rbind(tier_3_4_5, other_M)

# write selex out to use in the make_ancillary_prm.R script
# write.csv(selex, 'data/age_at_selex.csv', row.names = F)

# Calculate selected virgin biomass
# Because Atlantis will apply mFC as knife-edge, set all classes < age at selex as 0 biomass
biom_age_selected <- biom_age %>%
  left_join(selex, by = 'Code') %>%
  mutate(idx = as.numeric(Age) - as.numeric(age_class)) %>%
  filter(is.na(idx) | idx >= 0) %>%
  group_by(Code) %>%
  summarise(Biomass = sum(Biomass)) %>%
  ungroup()

# get catch based on initial biomass (total or selected) and F 

selected <- TRUE

if(selected){
  dat <- biom_age_selected %>%
    left_join(all, by = 'Code') %>%
    mutate(FMSY = replace_na(FMSY, 0),
           FMSY_25 = FMSY/4,
           FMSY_50 = FMSY/2,
           mu = 1-exp(-FMSY_25), # proportion of exploited population, CHANGE THIS FOR DIFFERENT FRACTIONS OF FMSY
           Catch = Biomass * mu)
  
} else {
  dat <- biom %>%
    left_join(all, by = 'Code') %>%
    mutate(FMSY = replace_na(FMSY, 0),
           FMSY_25 = FMSY/4,
           FMSY_50 = FMSY/2,
           mu = 1-exp(-FMSY_25), # proportion of exploited population, CHANGE THIS FOR DIFFERENT FRACTIONS OF FMSY
           Catch = Biomass * mu)
}

# we need to add invertebrate catch. Since I don't have a sense of what the equivalent of FMSY (or M) should be, let's base this one off of catch time series 1990-2020
# values will be in mg N s-1 for the first day of each month
# go from that to catch mt day, then * 30 to estimate monthly catch, then sum to get annual catch, then take average of annual catch, then sum across all boxes
catch_files <- c(list.files('C:/Users/Alberto Rovellini/Documents/GOA/Catch_data/output/AKFIN/', full.names = T), 
                 list.files('C:/Users/Alberto Rovellini/Documents/GOA/Catch_data/output/DFO/', full.names = T))

# extract column names
columns <- readLines(catch_files[1])
columns <- columns[grepl('long_name',columns)]
columns <- gsub('## COLUMN..long_name ', '', columns)
columns <- gsub('## COLUMN...long_name ', '', columns)

# which of these are inverts?
inverts <- fg %>% filter(NumCohorts == 1) %>% pull(Code)
caught_inverts <- intersect(columns, inverts)

origin <- as.Date('1991-01-01')

catch_by_box <- list()

for(i in 1:length(catch_files)){
 
  this_catch <- read.table(catch_files[i], skip = 309)
  colnames(this_catch) <- columns
  
  this_box <- as.numeric(gsub('.ts', '', gsub('.*catch', '', catch_files[i])))
  
  this_catch_long <- this_catch %>%
    pivot_longer(-Time, values_to = 'Catch_mgN_day', names_to = 'Code') %>%
    filter(Code %in% caught_inverts) %>%
    mutate(Box = this_box,
           Date = origin + Time,
           Year = year(Date),
           Month = month(Date),
           Catch_mt_day = Catch_mgN_day * 20 * 5.7 * 60 * 60 * 24 / 1e9,
           Catch_mt_month = Catch_mt_day * 30) %>%
    select(Box, Year, Month, Code, Catch_mt_month)
  
  catch_by_box[[i]] <- this_catch_long
}

catch_inverts <- do.call('rbind', catch_by_box) %>%
  group_by(Box, Year, Code) %>%
  summarise(Catch_mt = sum(Catch_mt_month)) %>%
  group_by(Box, Code) %>%
  summarise(Catch_mt = mean(Catch_mt)) %>%
  group_by(Code) %>%
  summarise(Catch_mt = sum(Catch_mt)) %>%
  ungroup() %>%
  mutate(Catch_background = Catch_mt / 4) %>% # assuming they have been fished at MSY 1991-2020, probably wrong # proportion of exploited population, CHANGE THIS FOR DIFFERENT FRACTIONS OF FMSY
  filter(Catch_background > 0) %>% # we have no catches for KIN, BFF, COR, PWN, SPG, although they are indicated as impacted in the Groups.csv file
  select(Code, Catch_background) 

# create a dummy file for catch from background F: 1 year, 1 month, 1 box (basically we already aggregate over time and space)
bgF <- dat %>%
  mutate(year = 1990,
         gear_name = 'Background F',
         box = 1,
         month = 1) %>% # this is dummy
  select(year, gear_name, Code, box, month, Catch) %>%
  rename(atlantis_fg = Code, tot_catch_mt = Catch) %>%
  mutate(catch_tons = tot_catch_mt) 

# bind the invertebrate catches, after formatting them like the dummy template
catch_inverts <- catch_inverts %>%
  mutate(year = 1990,
         gear_name = 'Background F',
         box = 1,
         month = 1,
         catch_tons = Catch_background) %>%
  rename(tot_catch_mt = Catch_background, atlantis_fg = Code) %>%
  select(names(bgF))

bgF <- rbind(bgF, catch_inverts)

# remove duplicates for inverts by summarising
bgF <- bgF %>%
  group_by(across(c(year:month))) %>%
  summarise(tot_catch_mt = sum(tot_catch_mt),
            catch_tons = sum(catch_tons)) %>%
  ungroup()

# add in all other fleets, for 1990, month 1, box 1, all catch 0
month <- 1
box <- 1
year <- 1990
all_fleets <- expand.grid(year, fleet_names, fg_codes, box, month)
all_fleets <- all_fleets %>%
  mutate(tot_catch_mt = 0, catch_tons = 0) %>%
  set_names(names(bgF)) %>%
  arrange(factor(gear_name, levels = fleet_names)) %>%
  filter(gear_name != 'Background F')

# bind
fleets_tmp <- rbind(bgF, all_fleets)

write.csv(fleets_tmp, 'data/all_fisheries_goa_25FMSY_OY.csv', row.names = F)
