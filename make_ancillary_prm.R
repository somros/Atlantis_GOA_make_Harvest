# Script to write the ancillary prm files that are used to populate the Harvest.prm file
# Alberto Rovellini
# 01/04/2022
# At first, this will be geared toward building a template / placeholder with all GOA fleets inactive, and only the background F active,
# for the purposes of completing model calibration
# Note: many of these vectors are only needed for dynamic fishing, but we still need them for the model to run
# UPDATE 2/15/2023 Adding XXX_mFC_startage vectors based on age at first selex instead of 0

library(tidyverse)

# groups
grp <- read.csv('data/GOA_Groups.csv')
grp.codes <- grp %>% pull(Code)
n.grp <- length(grp.codes)

# fleets
fleets <- read.csv('data/GOA_fleet_values.csv')
fleet.codes <- fleets %>% pull(Code)
n.fleets <- length(fleet.codes)

# boxes
bgm  <- readLines('data/GOA_WGS84_V4_final.bgm') # bgm file
n.box    <- as.numeric(gsub('nbox', '', grep("nbox", bgm, value = TRUE)))

# selex
selex <- read.csv('data/age_at_selex.csv')

# q_vector.prm
# Catchability. Needed for all groups, as many entries as there are fisheries
file.create('data/q_vector.prm')

for(i in 1:n.grp){
  cat(paste0('q_', grp.codes[i], ' ', n.fleets), file='data/q_vector.prm', append=TRUE,'\n')
  cat(rep(0, n.fleets), file='data/q_vector.prm', append=TRUE, '\n')
}

# target_vector.prm
# seems to only be used for dynamic fishing. Set all to 0. One vector per fleet, with as many entries as there are functional groups
file.create('data/target_vector.prm')

for(i in 1:n.fleets){
  cat(paste0('target_', fleet.codes[i], ' ', n.grp), file='data/target_vector.prm', append=TRUE,'\n')
  cat(rep(0, n.grp), file='data/target_vector.prm', append=TRUE, '\n')
}

# mpa_vector.prm
# Proportion of the box that is open to each fleet. 1 is fully open, 0 is fully closed
file.create('data/mpa_vector.prm')

for(i in 1:n.fleets){
  cat(paste0('MPA', fleet.codes[i], ' ', n.box), file='data/mpa_vector.prm', append=TRUE,'\n')
  cat(rep(1, n.grp), file='data/mpa_vector.prm', append=TRUE, '\n')
}

# startage_vector.prm
# first age class that is affected by fishing
# one vector per species with as many entries as there are fleets
# Option 1: Set all to 0, meaning that all age classes are vulnerable
# Option 2: Set to age at 50% selex or age at 50% maturity
file.create('data/startage_vector.prm')

selectivity <- FALSE # set this

if(selectivity){
  
  # for now this assumes that the fleet we work on is the first
  for(i in 1:n.grp){
    
    this_age_selex <- selex %>% filter(Code == grp.codes[i]) %>% pull(age_class)
    
    if(length(this_age_selex)==0) this_age_selex <- 0
    
    cat(paste0(grp.codes[i], '_mFC_startage', ' ', n.fleets), file='data/startage_vector.prm', append=TRUE,'\n')
    cat(c(this_age_selex, rep(0, n.fleets-1)), 
        file='data/startage_vector.prm', append=TRUE, '\n')
  }
  
} else {
  
  for(i in 1:n.grp){
    cat(paste0(grp.codes[i], '_mFC_startage', ' ', n.fleets), file='data/startage_vector.prm', append=TRUE,'\n')
    cat(rep(0, n.fleets), file='data/startage_vector.prm', append=TRUE, '\n')
  }
  
}

# agedistrib_vector.prm
# This defines the proportion of catch applied to each cohort
# the need to set up to 1
# Only relevant for imposed catch and for age-structured groups it seems
# for now, use dummy values
as.grp <- grp %>% filter(NumCohorts > 1)
as.grp.codes <- as.grp %>% pull(Code)
n.as.grp <- length(as.grp.codes)

file.create('data/agedistrib_vector.prm')

for(i in 1:n.as.grp){
  n.cohorts <- grp %>% filter(Code == as.grp.codes[i]) %>% pull(NumCohorts)
  
  cat(paste0('CatchTS_agedistrib', as.grp.codes[i], ' ', n.cohorts), file='data/agedistrib_vector.prm', append=TRUE,'\n')
  cat(rep(1/n.cohorts, n.cohorts), file='data/agedistrib_vector.prm', append=TRUE, '\n')
}
