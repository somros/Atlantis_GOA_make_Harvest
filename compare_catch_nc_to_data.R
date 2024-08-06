# Alberto Rovellini
# 6/18/2024
# this script extracts and views catch in space from netcdf files
# allow for comparisons between runs (e.g., Base no MPA, MPA fully open, in the future partially closed MPAs, etc)
# allows for comparisons between a run and Adam's catch streams (aka "the truth")
# gets residuals between compared runs
# TODO: add to the function an option of getting a time series instead of the terminal average
# In fact, that should just be the default, and then any averaging could be done afterwards

# Important caveat:
# Adam allocated catches to boundary boxes and island boxes, likely as a result of a spatial join between statistical areas and Atlantis boxes
# For island boxes, catch may come from any number of the boxes around an island
# For boundary boxes, catch may be a misallocation due to the chunky box shape or it could legitimately be catch occurring off the shelf
# At the moment (6/12/2024) it seems unlikely that we will be able to fully re-run his analysis, so we need to make assumptions for both cases
# Island boxes: split catch among neighboring boxes based on the existing proportions
# boundary boxes: assume that this is catch that occurred outside the Atlantis domain and drop it from the reconstruction

pacman::p_load(tidync, ncdf4, tidyverse, sf, rbgm, viridis, here, RColorBrewer)

# runs to compare and other settings
old_run <- 1517
new_run <- 1574
select <- dplyr::select
dir <- here()
data_dir <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/" # directory with the Atlantis runs

# Read data ---------------------------------------------------------------

# nc files
catch_nc_file_old <- paste0(data_dir, "/out_", old_run, "/outputGOA0", old_run, "_testCATCH.nc") # base
catch_nc_file_new <- paste0(data_dir, "/out_", new_run, "/outputGOA0", new_run, "_testCATCH.nc") # mpa

# Atlantis groups file
grps <- read.csv(here(dir, "data/GOA_Groups.csv"))
all_fg <- grps %>% filter(GroupType %in% c("FISH", "SHARK")) %>% pull(Name)

# geometry
goa_bgm <- read_bgm(here(dir, "data/GOA_WGS84_V4_final.bgm"))
goa_sf <- goa_bgm %>% box_sf()

# catch reconstruction from data
catch_dat <- readRDS(here(dir, "fleets/fleet_total_catch_atl.RDS"))
fleet_key <- read.csv(here(dir, "data/GOA_fisheries.csv"))

# load functions
source("catch_compare_functions.R")

# handle the catch reconstruction
catch_data <- prepare_catch_data(catch_dat = catch_dat, do_boundary = TRUE, do_islands = TRUE)

# which groups do you want to view in the plot comparisons?
to_plot <- c("Pollock","Cod","Arrowtooth_flounder","Flatfish_shallow","Flatfish_deep","Rex_sole","Flathead_sole","Pacific_ocean_perch","Sablefish","Halibut")
# which fleets do you want to see?
to_drop <- c("Canada", # none of these appears in the data
             "background",
             "historical",
             paste0("dummy", 1:9))
to_keep <- setdiff(fleet_key$Code, to_drop)

# make a folder for the plots
plotdir <- here(dir, "fleets", "catch_plots", new_run)
dir.create(plotdir, recursive = T)

# total catch comparison plots --------------------------------------------

scalars <- plot_total_catch(nc_old = catch_nc_file_old, 
                 nc_new = catch_nc_file_new, 
                 fleet_struc = F, 
                 relative = T, 
                 old_run = old_run, 
                 new_run = new_run, 
                 key = fleet_key, 
                 plotdir = plotdir,
                 write_scalars = T)

# optionally calibrate the mFC scalars based on the difference between the two runs

# which harvest.prm file do you want to work on
# harvest_old <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1573/GOA_harvest_fleets_mpa_CASE2_10YR.prm"
# harvest_new <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1573/GOA_harvest_fleets_mpa_CASE2_10YR_calibrated.prm"
# 
# calibrate_mfc_total(harvest_old, harvest_new, scalars) # this may lead to incorrect spatial patterns; and it may have no effect if a fishery fails in a box because its target isn't there

# Compare spatial patterns in the catch between run and data --------------

# This only makes sense either in relative terms or when the model operates under realistic F

# by species
plot_spatial_catch(nc_new = catch_nc_file_new, 
                   fleet_struc = F, 
                   relative = T, 
                   new_run = new_run, 
                   key = fleet_key, 
                   plotdir = plotdir,
                   write_scalars = T, 
                   catch_data = catch_data, 
                   by_species = T)
# by fleet
scalars <- plot_spatial_catch(nc_new = catch_nc_file_new, 
                              fleet_struc = T, 
                              relative = T, 
                              new_run = new_run, 
                              key = fleet_key, 
                              plotdir = plotdir,
                              write_scalars = T,
                              catch_data = catch_data, 
                              by_species = F)

# optionally calibrate the MPAYYY scalars based on the difference between the two runs
# which harvest.prm file do you want to work on?
# TODO: what is happening with BC here? Are you counting it in getting the scalars by box?
# harvest_old <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1562/GOA_harvest_fleets_mpa_v3.prm"
# harvest_new <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1562/GOA_harvest_fleets_mpa_v4.prm"
# 
# calibrate_mfc_total(harvest_old, harvest_new, scalars) # this may lead to total catches that are too high or too low


# Catch by port comparison plots ------------------------------------------
# view catch by port as model output compared to catch reconstruction (last 10 year average)
# This hinges on the matching key built in make_box_to_port_key.R
# read in matching key
box_port_key <- readRDS("data/box_to_port_key.RDS")

data_dir <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/" # directory with the Atlantis runs

# nc file
catch_nc_file_TOT <- paste0(data_dir, "/out_", new_run, "/outputGOA0", new_run, "_testTOTCATCH.nc") # full nc catch file

plot_port_catch(nc_tot = catch_nc_file_TOT, run = new_run, key = box_port_key)

# Catch composition comparison plots --------------------------------------

# Check that catch split by fleet is similar to the data
# by species
plot_catch_composition(nc_new = catch_nc_file_new, 
                       fleet_struc = T, 
                       relative = T, 
                       new_run = new_run, 
                       key = fleet_key, 
                       plotdir = plotdir,
                       write_scalars = T, 
                       catch_data = catch_data, 
                       by_species = T)

# by fleet
plot_catch_composition(nc_new = catch_nc_file_new, 
                       fleet_struc = T, 
                       relative = T, 
                       new_run = new_run, 
                       key = fleet_key, 
                       plotdir = plotdir,
                       write_scalars = T, 
                       catch_data = catch_data, 
                       by_species = F)
