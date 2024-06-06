# Alberto Rovellini
# 4/5/2024
# This code gets proportion of the catch for each species by fleet
# It is based on Adam's metier analysis
# The idea is to use these proportions to break down F (or mFC), whichever mFC we will use
# Some applications:
# break down background mFC into mFC by fleet
# break down FMSY, Ftarg, Fref, etc. by fleet

library(tidyverse)

# Read in data ------------------------------------------------------------

# read in RDS object from Adam
# read fleets and fleet keys
fleets <- readRDS("fleets/fleet_total_catch_atl.RDS")
fleet_key <- read.csv("fleets/fleet_Atlantis.csv")
grps <- read.csv("data/GOA_Groups.csv")

# remove space (model-wide F)
fleets_tot <- fleets %>%
  group_by(fleet, year, spp) %>%
  summarise(mt = sum(weight_mton))

# expand with all missing combinations
template <- expand.grid(fleet = unique(fleets_tot$fleet), 
                        year = unique(fleets_tot$year),
                        spp = unique(fleets_tot$spp))
fleets_tot <- merge(template, fleets_tot, by = c("fleet","year","spp"), all.x = T)
fleets_tot$mt[is.na(fleets_tot$mt)] <- 0

# for each species, get proportion of total catch by fleet
fleets_prop <- fleets_tot %>%
  group_by(year, spp) %>%
  mutate(tot_mt = sum(mt)) %>%
  ungroup() %>%
  mutate(prop = mt / tot_mt)

# is this summing up to 1?
# check <- fleets_prop %>%
#   group_by(year, spp) %>%
#   summarise(check = sum(prop))
# yes

# pick a point in time:
# say average 2015-2020
# This assumes that proportional effects of different fleets is time-invariant: pot caught x% of cod in 2020, so it will in 2080, so it did in 1995.
fleets_prop_term <- fleets_prop %>%
  filter(year >= 2015, year <= 2020) %>%
  group_by(spp, fleet) %>%
  summarise(prop = mean(prop, na.rm = T)) %>%
  ungroup()

# check
# check <- fleets_prop_term %>%
#   group_by(spp) %>%
#   summarise(check = sum(prop))
# okay

# view
fleets_prop_term %>%
  ggplot(aes(x = spp, y = prop, fill = fleet))+
  geom_bar(stat = "identity", position = "stack")+
  theme_bw()

# write this out
# write.csv(fleets_prop_term, "data/mFC_prop_by_fleet.csv", row.names = F)

# Apply to background mFC -------------------------------------------------

# open the harvest.prm from the Base model
# Only for the species that appear here, break down mFC across fleets
# make sure that the order of the fleets is consistent with the fleets.csv input file

file_path <- "C:/Users/Alberto Rovellini/Documents/GOA/Atlantis_GOA_OY_MS/"
bg_harvest_file <- paste0(file_path, "GOA_harvest_background.prm")
goa_fleets_file <- paste0(file_path, "GOA_fisheries.csv")

# read in
bg_harvest <- readLines(bg_harvest_file)
goa_fisheries <- read.csv(goa_fleets_file)

# fix typo
goa_fisheries$Code <- gsub("CgOthSpiKo", "CgOthSpiKi", goa_fisheries$Code)

# get codes, and rewrite fleet codes to match format
fg_codes <- unique(fleets_prop_term$spp)
fleet_codes_csv <- goa_fisheries %>% pull(Code)
fleet_codes <- fleets_prop_term %>% dplyr::select(fleet) %>% distinct()

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

fleet_codes <- fleet_codes %>% 
  rowwise() %>%
  mutate(fleet_key = rewrite_codes(as.character(fleet))) %>%
  ungroup()

# Add Canada fleet
# For now we are keeping that to mFC
# No other fleets catch in Canada and Canada fleet does not catch anywhere else so this should be fine
fleet_codes <- rbind(fleet_codes,
                     data.frame("fleet" = "Canada", "fleet_key" = "Canada"))

for(i in 1:length(fg_codes)) {
  
  grp <- fg_codes[i]
  idx <- grep(paste0('mFC_', grp, ' 33'), bg_harvest)
  
  parname <- bg_harvest[idx]
  parvals <- bg_harvest[idx + 1] %>% 
    strsplit(' ') %>% 
    unlist() %>% 
    as.numeric() %>% 
    matrix(nrow=1) %>%
    data.frame() %>%
    select_if(~ !any(is.na(.))) %>%
    set_names(fleet_codes_csv)
  
  # bg mfc
  bg_mfc <- parvals[1,1]
  
  # get modifiers by fleet
  scalars <- fleets_prop_term %>%
    filter(spp == grp) %>%
    dplyr::select(-spp)
  
  # add Canada
  scalars <- rbind(scalars, data.frame("fleet" = "Canada", prop = 1))
  
  # replace line in prm file
  parvals_long <- parvals %>%
    pivot_longer(everything(), names_to = "fleet_key", values_to = "mFC") %>%
    left_join(fleet_codes) %>%
    left_join(scalars, by = "fleet") %>%
    distinct() %>%
    mutate(prop = replace_na(prop, 0)) %>%
    rowwise() %>%
    mutate(mFC_new = bg_mfc * prop) %>%
    ungroup()
  
  # check
  check <- (parvals_long %>% pull(mFC_new) %>% sum()) - 2 * bg_mfc # this is 2*mFC because mFC for BC stays the same
  if(check != 0){
    print(paste(grp, ": new mFC - original mFC =", check))
  }
  
  # new mFC vector
  new_parvals <- parvals_long %>% pull(mFC_new) %>% paste(collapse = " ")
  
  # replace the string in the prm
  bg_harvest[idx + 1] <- new_parvals
  
}

# throughout, fix typo
bg_harvest <- gsub("CgOthSpiKo", "CgOthSpiKi", bg_harvest)

# Changing other parameters -----------------------------------------------
# Other parameters to change are:
# flagF_XXX
for(i in 1:length(fg_codes)) {
  
  grp <- fg_codes[i]
  idx <- grep(paste0("flagF_", grp, " 59"), bg_harvest)
  
  parname <- bg_harvest[idx]
  parvals <- bg_harvest[idx + 1] %>% 
    strsplit(' ') %>% 
    unlist() %>% 
    as.numeric()
  
  # all fleets can catch some small proportion of a species
  # open all of them, and let mFC determine how much of a species is caught
  # but, close the first two (bg and imposed catch), to minimize unintended behavior
  # if mFC is 0, catch will be 0

  # new flagF vector
  new_parvals <- c(0, 0, rep(1, (length(fleet_codes_csv)-2))) %>% paste(collapse = " ")
  
  # replace the string in the prm
  bg_harvest[idx + 1] <- new_parvals
  
  # also fix typo in the parmaeter name
  bg_harvest[idx] <- paste0("flagF_", grp, " 33")
  
}

# XXX_mFC_startage
# for now, set them all to the value of the background, then we will need fleet-specific selex
for(i in 1:length(fg_codes)) {
  
  grp <- fg_codes[i]
  idx <- grep(paste0(grp, '_mFC_startage 33'), bg_harvest)
  
  parname <- bg_harvest[idx]
  parvals <- bg_harvest[idx + 1] %>% 
    strsplit(' ') %>% 
    unlist() %>% 
    as.numeric()
  
  # bg mfc_startage
  bg_mfc_startage <- parvals[1]
  
  # new mFC_startage vector
  new_parvals <- rep(bg_mfc_startage, length(fleet_codes_csv)) %>% paste(collapse = " ")
  
  # replace the string in the prm
  bg_harvest[idx + 1] <- new_parvals
  
}

# write out the new prm
writeLines(bg_harvest, con = "data/GOA_harvest_fleets_v2.prm")
