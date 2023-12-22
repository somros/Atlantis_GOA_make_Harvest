# 12/21/2023
# Script to set up a new harvest.prm file based on FMP values for T3 stocks
# What happened here is that the selectivity knife edges have been adjusted
# Stock productivity has also changed, but tuning steepness and limiting unselected spawning
# These two are enough reason why the originla calibration of mFC may no longer be valid
# For example, previous misspecification of the selectivity curve may have been the cause of the insufficient fishing on some of the flatfish
# instead of tuning for that by increasing mFC, that has been addressed by fixing the selectivity curve, thereby making previous mFC tuning unjustified
# We will keep other stocks to the default values of the Base model (run 1328)

# This script build a new harvest.prm file where mFC values for the T3 stocks are pulled fresh from the FMP single-species FOFL values
# This will constitute the starting point for some rounds of mFC tuning

library(tidyverse)
library(readxl)

nfleet <- 33 # how many fleets in the harvest.prm file?

# get the default harvest.prm from 1473. It is important that this file has the correct startage setup for selectivity
# we will use this as the basis for all non-tier 3 stocks, because selectivity, steepness, and fecundity have not been changed as much
prm_file <- "../../Parametrization/output_files/data/out_1473/GOA_harvest_background.prm"
prm_vals <- readLines(prm_file)

fofl <- read_xlsx("data/GOA MSY estimates tables.xlsx", sheet = 1, range = 'A3:J19') %>%
  select(Stock, FOFL) %>%
  mutate(Code = c('POL','COD','SBF','FFS','FFS','FFS','FFS','FFD',
                  'REX','REX','ATF','FHS','POP','RFS','RFS','RFP')) %>%
  group_by(Code) %>%
  summarise(FOFL = mean(FOFL)) %>%
  ungroup() 

# fix cod with info from Pete Hulson:
fofl[fofl$Code == "COD","FOFL"] <- 0.5100000

# add halibut (use M = 0.2 as a proxy for FMSY)
fofl <- rbind(fofl, data.frame("Code"="HAL","FOFL"= 0.2))

# now get 1/4 FOFL
fofl <- fofl %>%
  mutate(bg_f = FOFL / 4)

# now turn F to mFC
# Formula for mFC is mfc = 1-exp(-F / 365)
fofl <- fofl %>%
  mutate(mfc = 1 - exp(-bg_f / 365))

# now write out a prm with these new values
# do one stock at a time
stocks <- unique(fofl$Code)

# create file
newfile <-  "mFC_tuning/GOA_harvest_background_OY.prm"

file.create(newfile)

prm_new <- prm_vals
for(i in stocks){

  this_mfc <- fofl %>% filter(Code == i) %>% pull(mfc)

  # find the lines that have parameters for the species of interest
  # mFC_XXX 
  mfc_line <- grep(paste0("mFC_",i), prm_vals) + 1
  new_mFC <- paste(as.character(c(this_mfc, rep(0, nfleet-1))), collapse = ' ')
  
  # replace the line
  prm_new[mfc_line] <- new_mFC 
  
}
# write to file
writeLines(prm_new, con = newfile)

