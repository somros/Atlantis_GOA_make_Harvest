# Alberto Rovellini
# 08/02/2024
# Compare flat text to nc files for catch output
# This is needed if we want to use CATCH.nc for the plotting and calibration
# The Problem:
# I want to use CATCH.nc for plots and to calculate calibration factors 
# CATCH.nc is the only file that has catch output by species, fleet, and box
# This is what Isaac thinks about this file:
# I always start with the simple CatchPerFishery.txt (which is not spatial). I trust that the most.  Then I check to see which NC files agree with it.   (I have to admit I have not played with these catch NC files in a while, and I forgot one might in units of numbers). My recollection is that whatever  NC file gave us catch per species per polygon (but not per fleet) was correct; but that no NC file correctly gave catch per species per polygon per fleet.    I can dig around more in old files, but I definitely recommend using the txt files as (spatially aggregated) truth.
# So, if we want to use, we need to:
# 1. Compare catch as FC vars (t) from the nc file to the flat txt file - not spatial
# 2. Compare catch by box between CATCH.nc and CATCHTOT.nc (we have done this before)

run <- 1574 # no fleets

select <- dplyr::select
dir <- here()
data_dir <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/" # directory with the Atlantis runs

# Read data ---------------------------------------------------------------

# nc files
catch_nc_file <- paste0(data_dir, "/out_", run, "/outputGOA0", run, "_testCATCH.nc") # full nc catch file
catch_nc_file_TOT <- paste0(data_dir, "/out_", run, "/outputGOA0", run, "_testTOTCATCH.nc") # full nc catch file
catch_txt_file <- paste0(data_dir, "/out_", run, "/outputGOA0", run, "_testCatchPerFishery.txt") # full nc catch file

# Atlantis groups file
grps <- read.csv(here(dir, "data/GOA_Groups.csv"))
all_fg <- grps %>% filter(GroupType %in% c("FISH", "SHARK")) %>% pull(Name)

# geometry
goa_bgm <- read_bgm(here(dir, "data/GOA_WGS84_V4_final.bgm"))
goa_sf <- goa_bgm %>% box_sf()

# fleet key
fleet_key <- read.csv(here(dir, "data/GOA_fisheries.csv"))

# TXT vs NC by fleet ------------------------------------------------------
# The txt file is not spatially-explicit, but let's see if the FC variables are even reliable at all

# read txt first
catch_txt <- read.table(catch_txt_file, sep = " ", header = T)
catch_txt <- catch_txt %>%
  pivot_longer(-c(Time, Fishery), names_to = "Code", values_to = "mt") %>%
  mutate(year = ceiling(Time / 365)) %>%
  rename(fleet = Fishery) %>%
  select(year, fleet, Code, mt)

# now nc file
catch_nc_ls <- list()

for(n in 1:length(all_fg)){
  fg <- all_fg[n] # this needs to use the "Name" to pull from the NC file
  this_tidync <- tidync(catch_nc_file)
  this_nc <- ncdf4::nc_open(catch_nc_file)
  
  code <- grps %>% filter(Name == fg) %>% pull(Code)
  
  #Extract from the output .nc file the appropriate catch time series variables
  catch_vars <- this_tidync %>%
    activate("D1,D0") %>%
    hyper_vars() %>% # all variables in the .nc file active grid
    filter(grepl("_Catch_FC.",name)) %>% # filter for reserve N
    filter(grepl(code,name)) # filter for specific functional group
  
  catch_vars <- catch_vars %>%
    mutate(Index = as.numeric(gsub(paste0(code, "_Catch_FC"), "", name))) %>%
    left_join(fleet_key %>% dplyr::select(Code, Index))
  
  catch <- purrr::map(catch_vars$name,ncdf4::ncvar_get,nc=this_nc) 
  
  catch_box_ls <- list()
  
  # Loop over each matrix in the collapsed list to fill the data frame
  for(i in 1:length(catch)) {
    # Get the current matrix
    mat <- catch[[i]]
    this_fleet <- catch_vars[i,]$Code
    
    # Convert matrix to a data frame
    mat_df <- as.data.frame((mat))
    
    # add box_id
    mat_df <- mat_df %>% mutate(box_id = 0:108)
    
    # reshape
    mat_df_long <- mat_df %>%
      pivot_longer(-box_id, names_to = "ts", values_to = "mt")
    
    # turn ts column to integer
    mat_df_long <- mat_df_long %>%
      mutate(ts = gsub("V","",ts)) %>%
      mutate(ts = as.numeric(ts)) %>%
      mutate(ts = ts - 1) # start numbering ts from 0
    
    # add fleet
    mat_df_long <- mat_df_long %>%
      mutate(fleet = this_fleet)
    
    catch_box_ls[[i]] <- mat_df_long
    
  }
  
  catch_box_df <- bind_rows(catch_box_ls) # this is in tons
  
  # drop ts = 0
  catch_box_df <- catch_box_df %>%
    filter(ts > 0)
  
  # add species
  catch_box_df <- catch_box_df %>%
    mutate(Name = fg)
  
  catch_nc_ls[[n]] <- catch_box_df
  
}

catch_nc_FC <- bind_rows(catch_nc_ls)

# collapse space for comparison with txt
catch_nc_FC_nospace <- catch_nc_FC %>%
  group_by(ts, fleet, Name) %>%
  summarise(mt = sum(mt, na.rm = T)) %>%
  left_join(grps %>% select(Code, Name)) %>%
  rename(year = ts) %>%
  select(year, fleet, Code, mt)

# compare
catch_comp <- catch_txt %>%
  filter(Code %in% unique(catch_nc_FC_nospace$Code)) %>%
  left_join(catch_nc_FC_nospace, by = c('year','fleet','Code')) %>%
  mutate(ratio = mt.x/mt.y)

summary(catch_comp$ratio)
hist(catch_comp$ratio, breaks = 100)

# at a level of total catch per species per fleet, CATCH.nc and CatchPerFishery.txt are very similar
# the discrepancies are as usual on migrating species, whose recording in the nc files is dodgy (and this will become important)
# so that's good!
# unfortunately, comparisons by box cannot be done (the txt file is spatially aggregated)
# no txt files for catch are by box
# We can compare catch in space between the two nc files

# Compare CATCH.nc to CATCHTOT.nc -----------------------------------------

# now nc file
catch_nc_tot_ls <- list()

for(n in 1:length(all_fg)){
  fg <- all_fg[n] # this needs to use the "Name" to pull from the NC file
  this_tidync <- tidync(catch_nc_file_TOT)
  this_nc <- ncdf4::nc_open(catch_nc_file_TOT)
  
  code <- grps %>% filter(Name == fg) %>% pull(Code)
  
  #Extract from the output .nc file the appropriate catch time series variables
  catch_vars <- this_tidync %>%
    activate("D1,D0") %>%
    hyper_vars() %>% # all variables in the .nc file active grid
    filter(grepl("Tot_.*_Catch",name)) %>% # filter for reserve N
    filter(grepl(code,name)) # filter for specific functional group
  
  # catch_vars <- catch_vars %>%
  #   mutate(Index = as.numeric(gsub(paste0(code, "_Catch_FC"), "", name))) %>%
  #   left_join(fleet_key %>% dplyr::select(Code, Index))
  
  catch <- purrr::map(catch_vars$name,ncdf4::ncvar_get,nc=this_nc) 
  
  catch_df <- as.data.frame((catch))
  
  # add box_id
  catch_df <- catch_df %>% mutate(box_id = 0:108)
  
  # reshape
  catch_df_long <- catch_df %>%
    pivot_longer(-box_id, names_to = "ts", values_to = "mt")
  
  # turn ts column to integer
  catch_df_long <- catch_df_long %>%
    mutate(ts = gsub("X","",ts)) %>%
    mutate(ts = as.numeric(ts)) %>%
    mutate(ts = ts - 1) # start numbering ts from 0
  
  catch_box_df <- catch_df_long
  
  # drop ts = 0
  catch_box_df <- catch_box_df %>%
    filter(ts > 0)
  
  # add species
  catch_box_df <- catch_box_df %>%
    mutate(Name = fg)
  
  catch_nc_tot_ls[[n]] <- catch_box_df
  
}

catch_nc_TOT <- bind_rows(catch_nc_tot_ls)

# aggregate fleets in FC frame
catch_nc_FC_nofleets <- catch_nc_FC %>%
  group_by(ts, box_id, Name) %>%
  summarise(mt = sum(mt, na.rm = T))

# compare
catch_comp_nc <- catch_nc_TOT %>%
  left_join(catch_nc_FC_nofleets, by = c('ts','box_id','Name')) %>%
  mutate(ratio = mt.x/mt.y)

summary(catch_comp_nc$ratio)

# the two nc files are identical in their reported catch by box

# So, conclusions for now are that:
# CatchPerFishery.txt and CATCH.nc contain the same catch by fleet FOR NON-MIGRATING SPECIES
# CATCH.nc and CATCHTOT.nc contain the same information per box

# Where does this leave us? is CATCH.nc safe to use, at least for non-migrating species?
# I guess we can't say that... Dig deeper into this with Isaac
# If the text file is the truth, let's compare the CATCHTOT.nc file to that, aka forget about the fleets
# This would allow us to make a matching key species-box-port, maybe, which is what we need in the end


# Compare CATCHTOT.nc to txt ----------------------------------------------
# This is a comparison of total catch, because CATCHTOT.nc does not hold fleet information

catch_comp_2 <- catch_txt %>%
  group_by(year, Code) %>%
  summarize(mt = sum(mt, na.rm = T)) %>%
  left_join(grps %>% select(Code, Name)) %>%
  filter(Name %in% unique(catch_nc_TOT$Name)) %>%
  left_join(catch_nc_TOT %>%
              rename(year = ts) %>%
              group_by(year, Name) %>%
              summarise(mt = sum(mt, na.rm=T)),
            by = c('year','Name')) %>%
  mutate(ratio = mt.x / mt.y)

summary(catch_comp_2$ratio)

# these are not identical, but fairly close for most. Exceptions are, as usual, migrating species
# So:
# 1. CatchPerFishery.txt is the truth, but it does not have info by box
# 2. Catch per fleet per species is comparable between CATCH.nc (FC variables) and CatchPerFishery.txt, but only for non-migrating species
# 3. Catch per box per species is identical between CATCH.nc and CATCHTOT.nc
# 4. Total catch per species is comparable between CATCHTOT.nc and CatchPerFishery.txt, but not for migrating species

# There is no txt file with information by box, therefore groundtruthing the spatial information in the NC files is impossible
# Apart from that, discrepancies concern migrating species

# For the purpose of creating a key boxes to ports, which NC file should we rely on?
# Should we just accept that reported catch by box is not reliable in any output?
# Question for Isaac. For now, make a key box-species-port from Adam's data, and use that to check outputs from CATCHTOT.nc (that is, remove the fleet layer)
 