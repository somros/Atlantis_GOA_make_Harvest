# Alberto Rovellini
# 04/04/2024
# This code will calculate recent F for each species for each fleet
# Here, fleets are the "metiers" from Adam's work
# To get F we will need catch (we have the catch reconstruction from Adam), and biomass estimates
# So far we have been calculating F only on the selected age classes
# It sure should not be on total biomass. Should it be on spawning biomass? It seemed like fishing intensity should be on selected age classes
# Because we end up calibrating mFC, this is possibly not as important, but check with Isaac
# Which biomass "data" should we use?
# We do not have reliable hindcast from Atlantis GOA that we can use as "internal" estimates
# This leaves us with stock assessment output that we used for model parametrization
# We need to leverage the same method we have used to parameterize the initial conditions to break biomass streams into age classes
# Because the life history parameters are time-invariant, we should take the decomposition we have already made (with the Canada approximation)
# And just carry it forward based on the trends from the stock assessment models
# We will be left with several species which are not assessed for which we need F
# ASSUMPTION: biomass in Canada follows the same trends as in Alaska. This is not true, but the other option is to perform a careful reconstruction with stock assessment biomass from AK and BC

# Ingredients:
# Catch streams from Adam
# Initial conditions for Base model (for age class biomass proportions), which reflect 1990
# Output from stock assessments
# age at selectivity

# Method:
# 1. Extract 1990 biomass by age class for all functional groups (we have functions that do that elsewhere)
# 2. Get estimates for following year assuming that proportions remain the same 9as do the life history params in Atlantis)
# 3. Filter out unselected age classes and aggregate biomass into total biomass of selected age classes
# 4. Do catch / biomass for each fleet to get harvest rate from each fleet
# 5. Convert to F and mFC
# 6. Take the last few years as an average
# 7. Compare to FOFL from assessments 

# a hack to skip data extraction from the NetCDF is to use the output of any (recent) model run and take biomass by age class at T=0. That is equivalent to the intial conditions

library(tidyverse)
library(tidync)
library(ncdf4)
library(readxl)
library(viridis)

# Read in data ------------------------------------------------------------

# read in RDS object from Adam
# read fleets and fleet keys
fleets <- readRDS("fleets/fleet_total_catch_atl.RDS")
fleet_key <- read.csv("fleets/fleet_Atlantis.csv")
grps <- read.csv("data/GOA_Groups.csv")

# read in initial conditions
init_file <- "data/GOA_cb_summer.nc"
init <- nc_open(init_file)
tidy_init <- tidync(init_file)
# skip the nc files and use the biomass output at t0
# biom_age <- read.table("data/outputGOA01517_testAgeBiomIndx.txt", sep = " ", header = T)
# ABORT - you cannot do this, because of... CANADA!!! :D :D :D

# read in stock assessment data streams
# Notes: 
# - these all have a different initial age
# - they only refer to Alaska. We assume that BC is an added quantity and that the dynamics are the same
# - some of the 4-yr cycle assessments streams stop in 2017. By now these have been updated, and we should use updated values for all. However, updating the streams will also require updating the life-history parameters...
biom_assessments <- read_excel("data/biomass_from_stock_assessment.xlsx", sheet = 1, range = "A1:AM54")

# read in selectivity information to subset age classes for F calculation
selex <- read.csv("data/age_at_selex_new.csv")


# # 1990 biomass by age class
# biom_age_1990 <- biom_age %>% 
#   filter(Time == 0) %>%
#   pivot_longer(-Time, names_to = "Code.Age", values_to = "mt") %>%
#   separate(Code.Age, into = c("Code", "Age"), sep = "\\.") %>%
#   mutate(Age = as.numeric(Age)) %>%
#   mutate(Time = 1990)

# Stock assessment output -------------------------------------------------
# UPDATE 04/05/2024
# Data streams from the most recent assessments are avaialble from https://apps-st.fisheries.noaa.gov/stocksmart?app=download-data
# Note, that the models used to produce these are different from the models we extracted LH info from to build the model
# Meaning, if we use these time series, they may not be representative of the "same" stocks
# However, if the purpose is to rescale 1990 Atlantis biomass based on relative changes in the assessment, we don't care about age classes. 
# Of course, the higher the lowest age in the assessment, the less precise this is. 
# It also always rests on the assumption that proportions among age classes are the same over time


# rework the assessment data frame
# long format, assign species to Codes, add up biomasses (where possible)
# express all biomasses as relative to 1990
biom_assessments_long <- biom_assessments %>%
  pivot_longer(-Year, names_to = "Stock", values_to = "mt")

# write a key for stocks to codes
key_ak <- data.frame('Stock'=c('Pacific cod',
                               'Alaska plaice',
                               'ATF',
                               'Butter sole',
                               'Dover sole', 
                               'Dusky rockfish',                   
                               'English sole',
                               'Flathead sole',
                               'Northern Rock Sole',                
                               'Pollock',
                               'Rex sole',
                               'Southern Rock Sole',                
                               'Starry flounder',
                               'Yellowfin sole',
                               'Northern rockfish',                 
                               'POP',
                               'Sablefish',  
                               'Halibut',   
                               'Giant Grenadier',
                               'Big skate',  
                               'Dogfish',   
                               'Longnose skate',
                               'Other skate',
                               'Thornyhead',
                               'Shortraker rockfish',
                               'Rougheye and blackspotted rockfish',
                               'Bigmouth sculpin',
                               'Great sculpin',
                               'Plain sculpin',
                               'Yellow Irish Lord',
                               'Harlequin rockfish',
                               'Redstripe rockfish',
                               'Sharpchin rockfish',
                               'Silvergray rockfish',
                               'Redbanded rockfish',
                               'Yelloweye rockfish',
                               'Yelloweye WGOA',
                               'Yelloweye EGOA'),
                     'Name'=c('Cod',
                              'Flatfish_shallow',
                              'Arrowtooth_flounder',
                              'Flatfish_shallow',
                              'Flatfish_deep',
                              'Rockfish_pelagic_shelf',
                              'Flatfish_shallow',
                              'Flathead_sole',
                              'Flatfish_shallow',
                              'Pollock',
                              'Rex_sole',
                              'Flatfish_shallow',
                              'Flatfish_shallow',
                              'Flatfish_shallow',
                              'Rockfish_slope',
                              'Pacific_ocean_perch',
                              'Sablefish',
                              'Halibut',
                              'Deep_demersal',
                              'Skate_big',
                              'Dogfish',
                              'Skate_longnose',
                              'Skate_other',
                              'Thornyhead',
                              'Rockfish_slope',
                              'Rockfish_slope',
                              'Sculpins',
                              'Sculpins',
                              'Sculpins',
                              'Sculpins',
                              'Rockfish_slope',
                              'Rockfish_slope',
                              'Rockfish_slope',
                              'Rockfish_slope',
                              'Rockfish_demersal_shelf',
                              'Rockfish_demersal_shelf',
                              'Rockfish_demersal_shelf',
                              'Rockfish_demersal_shelf'))

# join, group, add up
biom_assessments_long <- biom_assessments_long %>%
  left_join(key_ak, by = "Stock") %>%
  left_join(grps %>% select(Name, Code), by = "Name") %>%
  group_by(Year, Code, Name) %>%
  summarise(mt = sum(mt, na.rm = F)) %>%
  ungroup() %>%
  filter(Year >= 1990)

# subset to 1990
biom_assessment_1990 <- biom_assessments_long %>% 
  filter(Year == 1990) %>%
  rename(mt_1990 = mt) %>% 
  select(-Year)

# get relative biomass in the assessments
biom_assessment_relative <- biom_assessments_long %>%
  left_join(biom_assessment_1990) %>%
  mutate(relbiom = mt / mt_1990) %>%
  select(Year, Code, relbiom)

# Biomass in 1990 from Atlantis -------------------------------------------
# Do this one group at a time. It is worth doing it only for the groups tha are in the assessments for now
fg_to_do <- biom_assessment_relative %>%
  select(Code) %>%
  distinct() %>%
  left_join(grps %>% select(Code, Name)) %>%
  pull(Name)

# structure of the intial conditions file is very different from the structure of the output
# resN and structN are as fillvalues
# there is no depth distribution yet

biom_ak_ls <- list()

for(n in 1:length(fg_to_do)){
  
  print(paste("Doing", fg_to_do[n]))
  
  # args for the function below
  fg <- fg_to_do[n] # this needs to use the "Name" to pull from the NC file
  
  #Extract from the output .nc file the appropriate reserve N time series variables
  resN_vars <- hyper_vars(tidy_init) %>% # all variables in the .nc file active grid
    filter(grepl("_ResN",name)) %>% # filter for reserve N
    filter(grepl(fg,name)) # filter for specific functional group
  
  #Extract from the output .nc file the appropriate structural N time series variables
  strucN_vars <- hyper_vars(tidy_init) %>% # all variables in the .nc file active grid
    filter(grepl("_StructN",name)) %>% # filter for structural N
    filter(grepl(fg,name)) # filter for specific functional group
  
  # Get numbers by box
  abun_vars <- hyper_vars(tidy_init) %>% # all variables in the .nc file active grid
    filter(grepl("_Nums",name)) %>% # filter for abundance variables
    filter(grepl(fg,name)) # filter for specific functional group
  
  if(nrow(resN_vars)==0) {return("no data.")}
  else {
    # # Actually pull the data from the .nc
    # here we can collapse the depth layers but need to keep the boxes
    get_fillval <- function(ncfile, ncvar){
      this_fillval <- ncatt_get(ncfile, ncvar, "_FillValue")
      this_fillval <- this_fillval[[2]]
      return(this_fillval)
    }
    N <- data.frame("Name" = fg, 
                    "resN_var" = resN_vars$name, 
                    "strucN_var" = strucN_vars$name, 
                    "age" = 0:(length(resN_vars$name)-1))
    
    # extract reserve and structural nitrogen from the fillvalues
    N <- N %>% 
      rowwise() %>% 
      mutate(resN = get_fillval(ncfile = init, ncvar = resN_var),
             strucN = get_fillval(ncfile = init, ncvar = strucN_var)) %>%
      ungroup() %>%
      mutate(totN = resN + strucN,
             mt = totN * 20 * 5.7 / 1000000000)
    
    # now pull numbers
    # in init.nc they are all squished into one layer, they have not yet been distributed by vert parameters
    nums <-purrr::map(abun_vars$name,ncdf4::ncvar_get,nc=init) #numbers by age group,box,layer,time
    
    # collapse arrays
    nums_flat <- bind_rows(lapply(nums, function(x) as.data.frame(matrix(colSums(x, na.rm = T), nrow = 1))))
    nums_flat <- nums_flat %>%
      mutate(age = 1:nrow(.)) %>%
      mutate(age = age-1)
    
    nums_long <- nums_flat %>%
      pivot_longer(-age, names_to = "box_id", values_to = "nums") %>%
      mutate(box_id = gsub("V","",box_id),
             box_id = as.numeric(box_id)-1)
    
    # bring biomass in
    biom <- nums_long %>%
      left_join(N %>% select(age, mt), by = "age") %>%
      mutate(mt_box = mt * nums)
    
    # keep AK only
    biom_ak <- biom %>%
      filter(box_id < 92)
    
    # write out
    biom_ak <- biom_ak %>%
      mutate(Name = fg) %>%
      select(Name, age, box_id, mt_box)
    
    biom_ak_ls[[n]] <- biom_ak
    
  }
}

biom_age_1990 <- bind_rows(biom_ak_ls) # this is still by box
biom_age_1990 <- biom_age_1990 %>% left_join(grps %>% select(Code, Name)) # add species code

# Rescale Atlantis biomass based on stocks assessments --------------------

# this only makes sense for stocks that we have assessments for
biom_age_ts <- biom_assessment_relative %>%
  full_join(biom_age_1990) %>%
  mutate(mt_dyn = mt_box * relbiom)

# have a look
biom_age_ts %>%
  group_by(Year, Code, Name, age) %>%
  summarise(mt_dyn = sum(mt_dyn)) %>%
  ggplot(aes(x = Year, y = mt_dyn, color = factor(age)))+
  geom_line(linewidth = 1)+
  geom_point() +
  theme_bw()+
  scale_color_viridis_d()+
  facet_wrap(~Code, scales = "free")

# now bring in age at selectivity, and subset to the selected age classes only
biom_age_ts <- biom_age_ts %>%
  left_join(selex, by = "Code") %>%
  mutate(idx = age - age_class_selex) %>%
  filter(idx >= 0) %>%
  select(Year, Code, age, mt_dyn)

# now aggregate over ages
# this is the biomass of the selected age classes
biom_ts <- biom_age_ts %>%
  group_by(Year, Code) %>%
  summarise(mt_dyn = sum(mt_dyn, na.rm = F)) %>%
  filter(Year >= 2008)

# now join with the catch
# F is not box specific so first aggregate over space
# for now this only works with grouups for which we have biomass estimates
catch_biom <- fleets %>%
  filter(year >= 2008) %>%
  filter(spp %in% unique(biom_ts$Code)) %>%
  group_by(year, fleet, spp) %>%
  summarize(weight_mton = sum(weight_mton)) %>%
  left_join(biom_ts, by = c("spp"="Code", "year" = "Year"))

# get exploitation rate by fleet as catch/biomass
catch_biom <- catch_biom %>%
  mutate(mu = weight_mton/mt_dyn)

# get F / mFC:
# annual_mu = mean_catch / biomass_mt, # exploitation rate
# annual_F = (-1 * log(1-annual_mu)), # F
# mFC = 1-exp(-annual_F / 365)
catch_biom <- catch_biom %>%
  mutate(annual_F = (-1 * log(1-mu)),
         mFC = 1-exp(-annual_F / 365))

# how is F looking over time across fleets (last few years? For example for Cod?)
f_dyn <- catch_biom %>%
  group_by(year, spp) %>%
  summarize(f = sum(annual_F))

f_dyn %>%
  ggplot(aes(x = year, y = f))+
  geom_line()+
  geom_point()+
  theme_bw()+
  facet_wrap(~spp, scales = "free")

# get terminal F/mFC


# Relationships between F and catch ---------------------------------------

# We do not have biomass estimates for most functional groups
# In those cases, calculate the proportion of catch by fleet and apply that to the background F
# Before doing so, check that proportion of F and proportion of catch are equal for groundfish



# If this holds, I am starting to think that a far better approach would be:
# Pull F from Stock SMART
# Get proportion of catch by fleet from Adam's work
