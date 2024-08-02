# Alberto Rovellini
# 08/02/2024
# Compare flat text to nc files for catch output

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
  summarise(mt = sum(mt), na.rm = T) %>%
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
