
# Utility functions -------------------------------------------------------

#' renames fleet codes from the catch reconstruction to match the codes used by Atlantis
#'
#' @param original_string fleet name in the catch reconstruction 
#'
#' @description 
#' renames fleet codes from the catch reconstruction to match the codes used by Atlantis
#' @return name compliant with the fleet codes used by Atlantis
#' @export
#' 
#' 
#' 

# deal with codes for the fleets being different
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


# Manipulate catch data reconstruction ------------------------------------

#' perform some operations on the catch reconstruction from Adam
#'
#' @param catch catch reconstruction by fleet from Adam
#' @param do_boundary TRUE/FALSE do you want to turn the catch in boundary boxes to 0?
#' @param do_islands TRUE/FALSE do you want to redistribute the catch in island boxes to neighboring boxes based on catch proportions in the neighbors?
#'
#' @description 
#' Takes catch reconstruction and it modifies the catch values in boundary boxes and island boxes
#' These were non-zero in Adam's reconstruction, likely because of the the spatial join bertween the Atlantis geometry and the statistical areas 
#' We assume that:
#' Catch assigned from boundary boxes is effectively from areas off the shelf (1000m), so we discard it
#' Catch assigned to island boxes is to be redistributed to the neighbors of the island box based on catch proportions (by species, fleet, year) among the neighbors
#' @return catch reconstruction data frame
#' @export
#' 
#' 
#' 

prepare_catch_data <- function(catch_dat, do_boundary = FALSE, do_islands = FALSE){
  
  fleets <- catch_dat
  
  # zero-out BBs (TODO: thinks if this makes sense)
  if(do_boundary){
    
    bboxes <- goa_sf %>% filter(boundary == TRUE) %>% pull(box_id)
    fleets <- fleets %>%
      rowwise() %>%
      mutate(weight_mton = ifelse(box_id %in% bboxes,0,weight_mton)) %>%
      ungroup()
    
  }
  
  # handle islands: attribute catch from an island box to its neighbors based on the existing catch props
  if(do_islands){
    
    # 1. identify island boxes
    island_boxes <- goa_sf %>% filter(botz == 0) %>% pull(box_id) # 21, 44, 99 (Shumagin, Kodiak, Haida Gwaii)
    
    # 2. identify their neighbors
    neighbors <- st_touches(goa_sf) # identify THE INDEX of the touching geometries
    neighbors <- neighbors[island_boxes+1] # keep islands only
    
    # Initialize an empty data frame
    neighbors_key <- data.frame(island_box = integer(0), neighbor = integer(0))
    
    # Loop through the list and create the data frame
    for (i in seq_along(island_boxes)) {
      this_island_box <- island_boxes[i]
      these_neighbors <- neighbors[[i]] - 1 # turn from INDEX to box_id
      
      # Create a temporary data frame for the current focal geometry
      temp_df <- data.frame(island_box = this_island_box, neighbor = these_neighbors)
      
      # Append the temporary data frame to the result data frame
      neighbors_key <- rbind(neighbors_key, temp_df)
    }
    
    # nest this
    neighbors_key <- nest(neighbors_key, .by = neighbor, .key = "island_box")
    
    # 3. for each island box, year, fleet, species, port, get the catch value from the neighboring boxes
    # First subset the data to the boxes of interest - neighbors and their anchor islands
    # Then, for each island box and by accounting for grouping vars, look at the split of the catch between the neighbors
    fleets_islands <- fleets %>%
      left_join(neighbors_key, by = c("box_id" = "neighbor")) %>%
      unnest(cols = c(island_box)) 
    
    # 4. Turn said catch values into proportions
    fleets_islands_prop <- fleets_islands %>%
      group_by(year, island_box, fleet, spp, port) %>%
      mutate(total_weight = sum(weight_mton, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(prop = weight_mton / total_weight) %>%
      select(year, island_box, fleet, spp, port, box_id, prop) %>%
      rename(neighbors = box_id) %>%
      rename(box_id = island_box)
    
    # check that these all add up to 1
    # fleets_islands_prop %>%
    #   group_by(year, island_box, fleet, spp, port) %>%
    #   summarise(check = sum(prop)) %>%
    #   pull(check) %>%
    #   summary() # OK
    
    # 5. split catch from the island box into chunks based on calculated proportions
    fleets_islands_reallocated <- fleets_islands_prop %>%
      left_join(fleets, by = c("year", "box_id", "fleet", "spp", "port")) %>%
      mutate(weight_mton_2 = weight_mton * prop) %>% # this is the new weight to add to the neighbors
      drop_na()
    
    # now let's do some checks. The total catch in column weight_mton_2 summed by year needs to be the same as the raw catch
    # if it is not, we lost some of it on the way (for instance what if a metier/species.port combo fished in the island box but nowhere else? Then it would not appear here)
    # Expect differences
    # check1 <- fleets %>%
    #   filter(box_id %in% island_boxes) %>%
    #   group_by(year, fleet) %>%
    #   summarise(tot = sum(weight_mton))
    # 
    # check2 <- fleets_islands_reallocated %>%
    #   group_by(year, fleet) %>%
    #   summarise(tot = sum(weight_mton_2))
    # 
    # check3 <- left_join(check1, check2, by = c("year", "fleet")) %>%
    #   mutate(lost = tot.x - tot.y,
    #          lost_percent = (tot.x - tot.y) / tot.x * 100)
    
    # we lose some catch this way 
    # the worst cases are: 
    # 75% of CG_OTH_SMX_ROA in 1998
    # 65% of CG_OTH_SMX_ROA in 1997
    # 21% of CG_OTH_SSO_KI in 2019
    # 20% of CG_OTH_SMX_ROA in 1995
    # 17% of CG_OTH_SSO_KI in 2021
    # 4.4% of CG_OTH_SMX_ROA in 2006
    
    # 6. Add the new catch to the respective neighbors
    fleets_islands_reallocated <- fleets_islands_reallocated %>%
      select(year, neighbors, fleet, spp, port, weight_mton_2)
    
    fleets_2 <- fleets %>%
      left_join(fleets_islands_reallocated, by = c("year","box_id" = "neighbors", "fleet", "spp", "port")) %>%
      mutate(weight_mton_2 = replace_na(weight_mton_2, 0)) %>% # new weight to add is 0 if you are not the neighbor of an island
      mutate(weight_mton_3 = weight_mton + weight_mton_2) %>%
      select(-weight_mton, -weight_mton_2) %>%
      rename(weight_mton = weight_mton_3)
    
    # 7. set catch in all island boxes to 0
    fleets_2 <- fleets_2 %>%
      rowwise() %>%
      mutate(weight_mton = ifelse(box_id %in% island_boxes, 0, weight_mton)) %>%
      ungroup()
    
    # now check
    # check4 <- fleets %>%
    #   group_by(year, fleet) %>%
    #   summarise(tot  = sum(weight_mton))
    # 
    # check5 <- fleets_2 %>%
    #   group_by(year, fleet) %>%
    #   summarise(tot  = sum(weight_mton))
    # 
    # check6 <- left_join(check4,check5, by = c("year","fleet")) %>%
    #   mutate(lost = tot.x - tot.y,
    #          lost_percent = (tot.x - tot.y) / tot.x * 100)
    # small amount lost, revisit this in the future
    
    # replace the original fleet
    fleets <- fleets_2
    
  }
  
  return(fleets)
  
}


# Extract catch from Atlantis nc files ------------------------------------

# TODO: drop this or minimally rename them all such that this is v2 or something about age

#' Reads bio and catch nc files for a specific run
#'
#' @param catch_nc path to the catch NetCDF file
#' @param bio_nc path to the bio NetCDF file
#' @param age_struc TRUE/FALSE boolean - whether you want to see catch by age class 
#' @param relative TRUE/FALSE whether you want relative proportion of catch per box (T) or raw catch (F)
#' @param run number of calibration run
#'
#' @description 
#' 1. Reads bio and catch nc files for a specific run
#' 2. Calculates mean weight at age per species per boxm weighted by the number of individuals in the box / layer
#' 3. Extracts catch output per cell (in numbers) and converts it to tons
#' 4. Optionally, collapses across age classes for total catch
#' 5. Optionally, it turns it into a relative catch index to examine spatial distributions of catch#' 
#' @return data frame with catch information in space per species
#' @export
#' 
#' 
#' 

# build_catch_output <- function(catch_nc, bio_nc, age_struc, relative, run){
#   
#   # catch_nc = catch_nc_file_mpa
#   # bio_nc = bio_nc_file_mpa
#   # age_struc = F
#   # relative = T
#   # run = 1553
#   
#   # catch
#   this_tidync <- tidync(catch_nc)
#   this_nc <- ncdf4::nc_open(catch_nc)
#   #bio
#   this_tidync_bio <- tidync(bio_nc)
#   this_nc_bio <- ncdf4::nc_open(bio_nc)
#   
#   catch_nc_ls <- list()
#   
#   for(n in 1:length(all_fg)){
#     fg <- all_fg[n] # this needs to use the "Name" to pull from the NC file
#     
#     #Extract from the output .nc file the appropriate catch time series variables
#     catch_vars <- this_tidync %>%
#       activate("D1,D0") %>%
#       hyper_vars() %>% # all variables in the .nc file active grid
#       filter(grepl("_Catch",name)) %>% # filter for reserve N
#       filter(grepl(fg,name)) # filter for specific functional group
#     
#     #Extract from the output .nc file the appropriate reserve N time series variables
#     resN_vars <- hyper_vars(this_tidync_bio) %>% # all variables in the .nc file active grid
#       filter(grepl("_ResN",name)) %>% # filter for reserve N
#       filter(grepl(fg,name)) # filter for specific functional group
#     
#     #Extract from the output .nc file the appropriate structural N time series variables
#     strucN_vars <- hyper_vars(this_tidync_bio) %>% # all variables in the .nc file active grid
#       filter(grepl("_StructN",name)) %>% # filter for structural N
#       filter(grepl(fg,name)) # filter for specific functional group
#     
#     # Get numbers by box
#     # going to need these for weighting over depth layers (which we need to average over because catch is by box and not by depth)
#     abun_vars <- hyper_vars(this_tidync_bio) %>% # all variables in the .nc file active grid
#       filter(grepl("_Nums",name)) %>% # filter for abundance variables
#       filter(grepl(fg,name)) # filter for specific functional group
#     
#     if(nrow(resN_vars)==0) {return("no data.")}
#     else {
#       # # Actually pull the data from the .nc
#       # start from bio information
#       resN <- purrr::map(resN_vars$name,ncdf4::ncvar_get,nc=this_nc_bio) 
#       strucN <- purrr::map(strucN_vars$name,ncdf4::ncvar_get,nc=this_nc_bio)
#       nums <-purrr::map(abun_vars$name,ncdf4::ncvar_get,nc=this_nc_bio) #numbers by age group,box,layer,time
#       
#       # number of time steps
#       t_steps <- dim(resN[[1]])[3]
#       
#       # sum res and struct and then take averages over the water column, ending up with a value per box per time step
#       weighted_avg_list <- vector("list", length = length(resN)) # Initialize the list to store the weighted averages
#       
#       # Perform the operation (A+B) and then calculate the weighted average using C
#       weighted_mean_func <- function(x, w) {
#         sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
#       }
#       
#       weighted_waa <- list()
#       for(i in 1:length(resN)) {
#         # Calculate the sum of A and B
#         sum_ab <- (resN[[i]] + strucN[[i]]) * 20 * 5.7 / 1000000000 # from from mgN to tons
#         this_nums <- nums[[i]] # numbers at this age over the spatial structure
#         
#         # Initialize an array to store the results
#         this_weighted_waa <- array(dim = c(109, t_steps))
#         
#         # Calculate weighted mean across the first dimension
#         for(j in 1:109) {
#           for(k in 1:t_steps) {
#             # Extract the slice for the current i, j across all 7 values
#             sum_ab_slice <- sum_ab[, j, k]
#             nums_slice <- this_nums[, j, k]
#             
#             # Calculate the weighted mean for this slice
#             this_weighted_waa[j, k] <- weighted_mean_func(sum_ab_slice, nums_slice)
#           }
#         }
#         
#         # turn to df and pivot
#         this_weighted_waa <- this_weighted_waa %>%
#           as.data.frame() %>%
#           mutate(box_id = 0:108) %>%
#           pivot_longer(-box_id, values_to = "mt", names_to = "ts") %>%
#           mutate(ts = gsub("V","",ts)) %>%
#           mutate(ts = as.numeric(ts)) %>%
#           mutate(ts = ts - 1) %>%
#           mutate(age = i - 1) # count from 0
#         
#         weighted_waa[[i]] <- this_weighted_waa
#       }
#       
#       # turn it ito a df
#       weighted_waa_df <- bind_rows(weighted_waa)
#       
#       # drop t 0 and then take every 5th 
#       weighted_waa_df <- weighted_waa_df %>%
#         filter(ts > 0) %>%
#         filter(ts %in% seq(5,max(ts),5)) %>%
#         mutate(ts = ts / 5)
#       
#       # now let's extract the catch for this fg
#       # here we can collapse the depth layers but need to keep the boxes
#       catch <- purrr::map(catch_vars$name,ncdf4::ncvar_get,nc=this_nc) 
#       
#       # now turn to data frame
#       catch_box_ls <- list()
#       
#       # Loop over each matrix in the collapsed list to fill the data frame
#       for(i in 1:length(catch)) {
#         # Get the current matrix
#         mat <- catch[[i]]
#         
#         # Convert matrix to a data frame
#         mat_df <- as.data.frame((mat))
#         
#         # add box_id
#         mat_df <- mat_df %>% mutate(box_id = 0:108)
#         
#         # reshape
#         mat_df_long <- mat_df %>%
#           pivot_longer(-box_id, names_to = "ts", values_to = "nums")
#         
#         # turn ts column to integer
#         mat_df_long <- mat_df_long %>%
#           mutate(ts = gsub("V","",ts)) %>%
#           mutate(ts = as.numeric(ts)) %>%
#           mutate(ts = ts - 1) # start numbering ts from 0
#         
#         # add age
#         mat_df_long <- mat_df_long %>%
#           mutate(age = i-1) # number from 0 for consistency with age at selex and age mat
#         
#         catch_box_ls[[i]] <- mat_df_long
#         
#       }
#       
#       catch_box_df <- bind_rows(catch_box_ls) # this is in numbers... so we'd need to multiply this by WAA, pulling that in
#       
#       # drop ts = 0
#       catch_box_df <- catch_box_df %>%
#         filter(ts > 0)
#       
#       # now join and multiply to get mt per box
#       catch_box_df <- catch_box_df %>%
#         left_join(weighted_waa_df) %>%
#         mutate(mt_tot = nums * mt) %>%
#         mutate(Name = fg)
#       
#     }
#     
#     catch_nc_ls[[n]] <- catch_box_df
#   }
#   
#   catch_nc_df <- bind_rows(catch_nc_ls)
#   
#   # turn NaN to NA
#   catch_nc_df$mt[is.nan(catch_nc_df$mt)] <- NA
#   catch_nc_df$mt_tot[is.nan(catch_nc_df$mt_tot)] <- NA
#   
#   # get the equilibrium for the last 5 years
#   catch_nc_eq <- catch_nc_df %>%
#     filter(ts > (max(ts)-5)) %>%
#     group_by(box_id, Name, age) %>% # get average of last 5 years of the run
#     summarise(mt_tot = mean(mt_tot, na.rm = T)) %>%
#     ungroup()
#   
#   # make a spatial version
#   catch_nc_spatial <- goa_sf %>%
#     select(box_id) %>%
#     full_join(catch_nc_eq %>%
#                 select(box_id, Name, age, mt_tot),
#               by = "box_id")
#   
#   # if age_struc, keep ages, otherwise sum across them
#   if(!age_struc){
#     catch_nc_spatial <- catch_nc_spatial %>%
#       group_by(box_id, Name) %>%
#       summarize(mt_tot = sum(mt_tot, na.rm = T)) %>%
#       ungroup()
#     
#     # if relative catch is wanted, rescale and get proportions
#     if(relative){
#       catch_nc_spatial <- catch_nc_spatial %>%
#         group_by(Name) %>%
#         mutate(mt_goa = sum(mt_tot, na.rm = T)) %>%
#         ungroup() %>%
#         mutate(prop = mt_tot / mt_goa)
#     }
#     
#   } else {
#     
#     # if relative catch is wanted, rescale and get proportions
#     if(relative){
#       catch_nc_spatial <- catch_nc_spatial %>%
#         group_by(Name, age) %>%
#         mutate(mt_goa = sum(mt_tot, na.rm = T)) %>%
#         ungroup() %>%
#         mutate(prop = mt_tot / mt_goa)
#       
#     }
#     
#   }
#   
#   # add run information
#   catch_nc_spatial <- catch_nc_spatial %>%
#     mutate(run = run)
#   
#   return(catch_nc_spatial)
#   
# }

#' Reads bio and catch nc files for a specific run
#'
#' @param catch_nc path to the catch NetCDF file
#' @param bio_nc path to the bio NetCDF file
#' @param fleet_struc TRUE/FALSE whether you want to see catch by fleet (T) or total (F) 
#' @param relative TRUE/FALSE whether you want relative proportion of catch per box (T) or raw catch (F)
#' @param run number of calibration run
#' @param key the fisheries.csv file from Atlantis 
#'
#' @description 
#' 1. Reads bio and catch nc files for a specific run
#' 2. Calculates mean weight at age per species per boxm weighted by the number of individuals in the box / layer
#' 3. Extracts catch output per cell (in numbers) and converts it to tons
#' 4. Optionally, collapses across age classes for total catch
#' 5. Optionally, it turns it into a relative catch index to examine spatial distributions of catch#' 
#' @return data frame with catch information in space per species
#' @export
#' 
#' 
#' 

# second version using the FC variables instead of catch 
# compare the output of the two
build_catch_output_v2 <- function(catch_nc, bio_nc, fleet_struc, relative, run, key){
  
  # catch
  this_tidync <- tidync(catch_nc)
  this_nc <- ncdf4::nc_open(catch_nc)
  
  catch_nc_ls <- list()
  
  for(n in 1:length(all_fg)){
    fg <- all_fg[n] # this needs to use the "Name" to pull from the NC file
    code <- grps %>% filter(Name == fg) %>% pull(Code)
    
    #Extract from the output .nc file the appropriate catch time series variables
    catch_vars <- this_tidync %>%
      activate("D1,D0") %>%
      hyper_vars() %>% # all variables in the .nc file active grid
      filter(grepl("_Catch_FC.",name)) %>% # _Catch_FC
      filter(grepl(code,name)) # filter for specific functional group
    
    # variables are labeled numerically - I assume in the order that fleets appear in the harvest file / fisheries.csv
    # contrary to regular Atlantis labeling, it seems as though indices start from 1 - though it may be because I have them that way in the fisheries.csv file
    catch_vars <- catch_vars %>%
      mutate(Index = as.numeric(gsub(paste0(code, "_Catch_FC"), "", name))) %>%
      left_join(key %>% select(Code, Index))
    
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
  
  catch_nc_df <- bind_rows(catch_nc_ls)
  
  # get the equilibrium for the last 5 years
  # catch_nc_eq <- catch_nc_df %>%
  #   filter(ts > (max(ts)-5)) %>%
  #   group_by(box_id, Name, fleet) %>% # get average of last 5 years of the run
  #   summarise(mt_tot = mean(mt, na.rm = T)) %>%
  #   ungroup()
  
  # make a spatial version
  # catch_nc_spatial <- goa_sf %>%
  #   select(box_id) %>%
  #   full_join(catch_nc_df %>%
  #               select(ts, box_id, Name, fleet, mt),
  #             by = "box_id")
  
  catch_nc_df <- catch_nc_df %>%
    select(ts, box_id, Name, fleet, mt)
  
  # if fleet_struc, keep fleets, otherwise sum across them
  if(!fleet_struc){
    catch_nc_df <- catch_nc_df %>%
      group_by(ts, box_id, Name) %>%
      summarize(mt = sum(mt, na.rm = T)) %>%
      ungroup()
    
    # if relative catch is wanted, rescale and get proportions
    if(relative){
      catch_nc_df <- catch_nc_df %>%
        group_by(ts, Name) %>%
        mutate(mt_goa = sum(mt, na.rm = T)) %>%
        ungroup() %>%
        mutate(prop = mt / mt_goa)
    }
    
  } else {
    
    # if relative catch is wanted, rescale and get proportions
    if(relative){
      catch_nc_df <- catch_nc_df %>%
        group_by(ts, Name, fleet) %>%
        mutate(mt_goa = sum(mt, na.rm = T)) %>%
        ungroup() %>%
        mutate(prop = mt / mt_goa)
    }
  }
  
  # add run information
  catch_nc_df <- catch_nc_df %>%
    mutate(run = run)
  
  return(catch_nc_df)
  
}


# Plotting functions ------------------------------------------------------

#' Plot catch from NC files comparing two Atlantis runs. 
#'
#' @param nc_old path to the netcdf file of the old run
#' @param nc_new path to the netcdf file of the new run
#' @param fleet_struc TRUE/FALSE whether you want to see results by fleet (T) or total (F)
#' @param relative TRUE/FALSE whether you want relative proportion of catch per box (T) or raw catch (F)
#' @param old_run number of previous run
#' @param new_run number of current run
#' @param key the fisheries.csv file from Atlantis 
#' @param plotdir path to the directory where the output plots will be stored
#' @param write_scalars writes out scalars to calibrate mFC based on differences between the two runs
#' 
#' @description 
#' Extracts catch output from Atlantis catch.nc files from two runs to compare. It performs comparisons of total catch and catch by box.
#' Optionally, it calculates the ration between the old and new run and gets
#' scalars that, if multiplied by the mFC values in the new run, will yield the catches in the old run. In this context, "old" is the target you strive for
#' while new is where you are now.
#' @return Optional data frame of scalars by species for mFC calibration
#' @export
#' 
#' 
#' 

plot_total_catch <- function(nc_old, nc_new, fleet_struc = F, relative = F, old_run, new_run, key, plotdir, write_scalars = F){
  
  catch_nc_old <- build_catch_output_v2(catch_nc = nc_old, 
                                        fleet_struc = fleet_struc,
                                        relative = relative,
                                        run = old_run,
                                        key = key)
  
  catch_nc_new <- build_catch_output_v2(catch_nc = nc_new, 
                                        fleet_struc = fleet_struc,
                                        relative = relative,
                                        run = new_run,
                                        key = key)
  
  
  # get residuals (it gets hazy for proportions - what's a big residual and how do you translate that to catch?)
  catch_diff <- catch_nc_old %>%
    left_join(catch_nc_new, by = c("ts", "box_id", "Name")) %>%
    select(ts, Name, mt_goa.x, mt_goa.y) %>%
    distinct() %>%
    rename(old = mt_goa.x, new = mt_goa.y) 
  
  # for plotting
  catch_diff_long <- catch_diff %>%
    pivot_longer(-c(ts,Name), names_to = "run", values_to = "mt")
  
  # view
  catch_diff_long %>%
    filter(Name %in% to_plot) %>%
    ggplot()+
    geom_line(aes(x = ts, y = mt, color = run), linewidth = 1.5)+
    facet_wrap(~Name, scales = "free")
  
  ggsave(paste0(plotdir, "/", new_run, "_vs_", old_run, ".png"), width = 10, height = 8)
  
  # One key question here is: How short a run can we get away with?
  # The plot above indicates that, while there can still be some divergence at the end of the run,
  # for most groups it is fairly apparent already at the start of the run that the fishery can't hit the non-spatial quota
  # So in the interest of ironing out the big kinks fast, let's do some bried (1-2 years) runs
  
  # for ease of visualization, add a plot for the first time step only (you will need it when you have 1-yr long runs)
  catch_diff_long %>%
    filter(Name %in% all_fg) %>%
    filter(ts == 1) %>%
    ggplot()+
    geom_bar(aes(x = run, y = mt), stat = "identity")+
    facet_wrap(~Name, scales = "free")
  
  ggsave(paste0(plotdir, "/", new_run, "_vs_", old_run, "_ts1.png"), width = 10, height = 8)
  
  # compare catch in space
  # get residuals (it gets hazy for proportions - what's a big residual and how do you translate that to catch?)
  catch_diff2 <- catch_nc_old %>%
    left_join(catch_nc_new, by = c("ts", "box_id", "Name")) %>%
    mutate(residual_mt = mt.x - mt.y,
           residual_prop = prop.x - prop.y)
  
  # add space
  catch_diff2 <- goa_sf %>%
    select(box_id) %>%
    left_join(catch_diff2, by = "box_id")
  
  # view
  
  catch_diff2 %>%
    filter(ts == 15, Name %in% to_plot) %>% # 15 just as test
    ggplot()+
    geom_sf(aes(fill = residual_mt))+
    scale_fill_viridis()+
    theme_bw()+
    facet_wrap(~Name, ncol = 3)
  ggsave(paste0(plotdir, "/", new_run, "_vs_", old_run, "_spatial.png"), width = 10, height = 6.5)
  
  # other view
  catch_diff2 %>%
    filter(ts == 15, Name %in% to_plot) %>%
    ggplot(aes(x = mt.x, y = mt.y, color = box_id))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")+  # Add 1:1 line
    labs(x = "Catch in non-spatial run", y = "Catch in spatial run")+
    facet_wrap(~Name, scales = "free")
  ggsave(paste0(plotdir, "/", new_run, "_vs_", old_run, "_by_box.png"), width = 10, height = 9)
  
  # Optionally get scaling factors to calibrate mFC
  if(write_scalars){
    
    scalars <- catch_diff %>%
      filter(ts == 1) %>%
      mutate(scalar = old / new) %>%
      select(Name, scalar) %>%
      left_join(grps %>% select(Name, Code))
    
    return(scalars)
    
  }
  
}

# TODO: could ALL these functions become one with args for which plot / slicing you want?

#' Compares catch from the nc file of one Atlantis run with data in space. 
#'
#' @param nc_new path to the netcdf file of the new run
#' @param fleet_struc TRUE/FALSE whether you want to see results by fleet (T) or total (F)
#' @param relative TRUE/FALSE whether you want relative proportion of catch per box (T) or raw catch (F)
#' @param old_run number of previous run
#' @param new_run number of current run
#' @param key the fisheries.csv file from Atlantis 
#' @param plotdir path to the directory where the output plots will be stored
#' @param write_scalars writes out scalars to calibrate mFC based on differences between the two runs
#' @param catch_data data frame with the catch reconstruction
#' @param by_species TRUE/FALSE whether you want to make plots by species (T) or by fleet (F)
#' 
#' 
#' @description 
#' Extracts catch output from Atlantis catch.nc files from one model run and compares it with data. 
#' It performs comparisons of catches in space, either by species or by fleet.
#' It is set up to operate on relative catch distributions rather than raw tonnage.
#' Optionally, it calculates the ration between the old and new run and gets scalars that, if multiplied by the MPAYYY values in the new run, 
#' will yield the catches in the old run. In this context, "old" is the target you strive for while new is where you are now.
#' It will produce maps.
#' @return Optional data frame of scalars by species for mFC calibration
#' @export
#' 
#' 
#'

plot_spatial_catch <- function(nc_new, fleet_struc = F, relative = F, new_run, key, plotdir, write_scalars = F, catch_data = catch_data, by_species = T){
  
  if(by_species){
    
    catch_nc_new <- build_catch_output_v2(catch_nc = nc_new, 
                                          fleet_struc = fleet_struc,
                                          relative = relative,
                                          run = new_run,
                                          key = key)
    
    # average of end of the run
    catch_nc_end <- catch_nc_new %>%
      filter(ts > (max(ts)-5)) %>%
      group_by(box_id, Name) %>%
      summarize(mt = mean(mt)) %>%
      ungroup() %>%
      group_by(Name) %>%
      mutate(tot_mt = sum(mt)) %>%
      ungroup() %>%
      mutate(prop = mt / tot_mt)
    
    # start wih no fleets here - space only
    catch_dat <- catch_data %>%
      ungroup() %>%
      left_join(grps %>% select(Code,Name), by = c("spp"="Code")) %>%
      select(year, box_id, Name, fleet, weight_mton) %>%
      mutate(weight_mton = replace_na(weight_mton, 0)) %>%
      group_by(year, box_id, Name) %>%
      summarise(mt = sum(weight_mton)) %>% # sum up across fleets
      ungroup() %>%
      filter(year > (max(year)-5)) %>%
      group_by(box_id, Name) %>%
      summarize(mt = mean(mt)) %>%
      ungroup()
    
    # add missing combinations
    dummy_df <- expand.grid("box_id" = 0:108, "Name" = unique(catch_dat$Name))
    catch_dat <- catch_dat %>%
      full_join(dummy_df) %>%
      mutate(mt = replace_na(mt, 0))
    
    # get relative catch
    catch_dat <- catch_dat %>%
      group_by(Name) %>%
      mutate(tot_mt = sum(mt, na.rm = T)) %>%
      ungroup() %>%
      mutate(prop = mt / tot_mt)
    
    # check differences
    catch_diff <- catch_dat %>%
      left_join(catch_nc_end, by = c("box_id", "Name"))# %>%
    # mutate(residual = mt.x - mt.y,
    #        ratio = mt.x / mt.y) # observed - predicted
    
    # view
    catch_diff <- catch_diff %>%
      select(box_id, Name, prop.x, prop.y) %>%
      rename(data = prop.x, model = prop.y) %>%
      pivot_longer(-c(box_id, Name), names_to = "Type", values_to = "Prop")
    
    # add space
    catch_diff <- goa_sf %>%
      select(box_id, boundary) %>%
      left_join(catch_diff, by = "box_id")
    
    # make a directory
    dir.create(paste0(plotdir, "/spatial/species/"), recursive = T)
    
    for(i in 1:length(all_fg)){
      this_fg <- all_fg[i]
      
      if(this_fg %in% unique(catch_diff$Name)){
        print(paste("doing", this_fg))
        
        catch_diff %>%
          filter(Name == this_fg) %>%
          ggplot()+
          geom_sf(aes(fill = Prop))+
          scale_fill_viridis()+
          theme_bw()+
          labs(title = this_fg)+
          facet_grid(Name~Type)
        ggsave(paste0(plotdir, "/spatial/species/data_vs_", new_run, this_fg, ".png"), width = 8, height = 4)
        
      }
      
    }
    
  } else { # do this by fleet
    
    if(!fleet_struc)stop("You want to plot fleet catch in space but you are not keeping the fleet structure")
    
    catch_nc_new <- build_catch_output_v2(catch_nc = nc_new, 
                                          fleet_struc = fleet_struc,
                                          relative = relative,
                                          run = new_run,
                                          key = key)
    
    # first - look at total catch across the model and species
    # Is the split across fleets the same in the data and the model?
    
    tot_atlantis <- catch_nc_new %>%
      filter(fleet != "background", fleet != "Canada") %>%
      filter(ts > (max(ts)-5)) %>%
      group_by(box_id, fleet, Name) %>%
      summarize(mt = mean(mt)) %>% # mean over last 5 years
      ungroup() %>%
      group_by(fleet) %>%
      summarize(mt = sum(mt)) %>% # sum over species and boxes
      ungroup() %>%
      mutate(tot = sum(mt)) %>% # total catch
      mutate(prop = mt / tot) %>%
      mutate(Type = "model") %>%
      select(Type, fleet, prop)
    
    bg <- catch_nc_new %>%
      filter(fleet == "background")
    
    tot_data <- catch_data %>%
      ungroup() %>%
      left_join(grps %>% select(Code,Name), by = c("spp"="Code")) %>%
      select(year, box_id, Name, fleet, weight_mton) %>%
      mutate(weight_mton = replace_na(weight_mton, 0)) %>%
      group_by(year, fleet) %>%
      summarise(mt = sum(weight_mton)) %>% # sum up across boxes and species
      ungroup() %>%
      filter(year > (max(year)-5)) %>% # keep last 5 years of data
      group_by(fleet) %>%
      summarize(mt = mean(mt)) %>% # average over last 5 years
      ungroup() %>%
      mutate(tot = sum(mt), # total catch
             prop = mt / tot) %>%
      rowwise() %>%
      mutate(fleet = rewrite_codes(fleet))%>%
      mutate(Type = "data") %>%
      select(Type, fleet, prop)
    
    catch_diff <- rbind(tot_atlantis, tot_data) %>%
      left_join(fleet_key %>% select(Code, Name), by = c("fleet"="Code"))
    
    # view
    # this one belongs with the total catch function above
    colourCount <- length(to_keep) 
    getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
    
    catch_diff %>%
      filter(fleet %in% to_keep) %>%
      ggplot(aes(x = Type, y = prop, fill = Name))+
      geom_bar(stat = "identity", position = "stack")+
      scale_fill_manual(values = getPalette(colourCount))+
      labs(title = "Fleet composition of the catch (data vs model")
    ggsave(paste0(plotdir, "/data_vs_", new_run, "_total_catch_byfleet.png"), width = 8, height = 4)
    
    # Now produce spatial plots but fleet-by-fleet
    # model output
    catch_nc_end <- catch_nc_new %>%
      filter(ts > (max(ts)-5)) %>% # keep last 5 years of run
      group_by(box_id, Name, fleet) %>%
      summarize(mt = mean(mt)) %>% # means across last 5 years
      group_by(box_id, fleet) %>%
      summarise(mt = sum(mt)) %>% # sum across species
      group_by(fleet) %>%
      mutate(tot_mt = sum(mt)) %>% # get total annual catch from a fleet
      ungroup() %>% 
      mutate(Prop = mt / tot_mt) %>% # get proportions by box
      mutate(Type = "model") %>%
      select(box_id, fleet, Prop, Type)
    
    # data
    catch_dat <- catch_data %>%
      ungroup() %>%
      left_join(grps %>% select(Code,Name), by = c("spp"="Code")) %>%
      select(year, box_id, Name, fleet, weight_mton) %>%
      mutate(weight_mton = replace_na(weight_mton, 0)) %>%
      filter(year >= (max(year)-5)) %>%
      group_by(box_id, Name, fleet) %>%
      summarize(mt = mean(weight_mton)) %>% # averages over time
      group_by(box_id, fleet) %>%
      summarise(mt = sum(mt)) %>% # sum across species
      ungroup()
    
    # add missing combinations
    dummy_df <- expand.grid("box_id" = 0:108, "fleet" = unique(catch_dat$fleet))
    catch_dat <- catch_dat %>%
      full_join(dummy_df) %>%
      mutate(mt = replace_na(mt, 0))
    
    # get relative catch
    catch_dat <- catch_dat %>%
      group_by(fleet) %>%
      mutate(tot_mt = sum(mt)) %>% # get total annual catch from a fleet
      ungroup() %>%
      mutate(Prop = mt / tot_mt) %>% # get proportions by box
      mutate(Type = "data") %>%
      select(box_id, fleet, Prop, Type) %>%
      rowwise() %>%
      mutate(fleet = rewrite_codes(fleet)) %>%
      ungroup()
    
    # bind data and model output
    catch_diff <- rbind(catch_nc_end, catch_dat) %>%
      left_join(fleet_key %>% select(Code, Name), by = c("fleet" = "Code"))
    
    # add space
    catch_diff_sf <- goa_sf %>%
      select(box_id) %>%
      left_join(catch_diff, by = "box_id")
    
    # make a directory
    dir.create(paste0(plotdir, "/spatial/fleet/"), recursive = T)
    
    # make plots
    for(i in 1:length(to_keep)){
      this_fleet <- to_keep[i]
      this_fleet_name <- catch_diff_sf %>%
        filter(fleet == this_fleet) %>%
        pull(Name) %>%
        unique()
      
      print(paste("doing", this_fleet))
      
      catch_diff_sf %>%
        filter(fleet == this_fleet) %>%
        ggplot()+
        geom_sf(aes(fill = Prop))+
        scale_fill_viridis()+
        theme_bw()+
        facet_grid(~Type)+
        labs(title = this_fleet_name)
      ggsave(paste0(plotdir, "/spatial/fleet/data_vs_", new_run, this_fleet, ".png"), width = 8, height = 4)
      
    }
    
    if(write_scalars){
      # based on the bit above, we can identify differences in proportional catch between data and model
      # In theory, you can rescale MPAYYY based on these discrepancies
      # I am  not convinced that relative quantities are going to work here, you could rescale biomasses but this seems odd
      # let's try it 
      
      scalars <- catch_diff %>%
        select(box_id, fleet, Type, Prop) %>%
        filter(fleet %in% to_keep) %>%
        pivot_wider(names_from = "Type", values_from = "Prop") %>%
        mutate(scalar = data / model) %>%
        rowwise() %>%
        mutate(scalar = ifelse(is.nan(scalar), 1, scalar)) %>% # turn NaN's to 1's, meaning no scaling
        ungroup()
      
      return(scalars)
    }
    
  }
  
  
}

# TODO: the args of this are identical to the one above - These should all be the same function

#' Compares catch compositions by species and fleet from the nc file of one Atlantis run and from the data. 
#'
#' @param nc_new path to the netcdf file of the new run
#' @param fleet_struc TRUE/FALSE whether you want to see results by fleet (T) or total (F)
#' @param relative TRUE/FALSE whether you want relative proportion of catch per box (T) or raw catch (F)
#' @param old_run number of previous run
#' @param new_run number of current run
#' @param key the fisheries.csv file from Atlantis 
#' @param plotdir path to the directory where the output plots will be stored
#' @param write_scalars writes out scalars to calibrate mFC based on differences between the two runs
#' @param catch_data data frame with the catch reconstruction
#' @param by_species TRUE/FALSE whether you want to make plots by species (T) or by fleet (F)
#' 
#' 
#' @description 
#' Extracts catch output from Atlantis catch.nc files from one model run and compares it with data. 
#' It performs comparisons of catch composition either by species or by fleet.
#' It is set up to operate on relative catch distributions rather than raw tonnage.
#' It will produce bar charts.
#' @return NA
#' @export
#' 
#' 
#'

plot_catch_composition <- function(nc_new, fleet_struc = T, relative = T, new_run, key, plotdir, write_scalars = F, catch_data = catch_data, by_species = T){
  
  catch_nc_new <- build_catch_output_v2(catch_nc = nc_new, 
                                        fleet_struc = fleet_struc,
                                        relative = relative,
                                        run = new_run,
                                        key = key)
  
  if(by_species){ # if you want to slice this by species
    
    # average of end of the run
    catch_nc_end <- catch_nc_new %>%
      filter(ts > (max(ts)-5)) %>%
      filter(box_id < 92) %>% # keep AK only
      group_by(fleet, Name) %>%
      summarize(mt = mean(mt)) %>%
      ungroup() %>%
      group_by(Name) %>% # group by species
      mutate(tot_mt = sum(mt)) %>%
      ungroup() %>%
      mutate(prop = mt / tot_mt)
    
    # no space - fleets only
    catch_dat <- catch_data %>%
      ungroup() %>%
      left_join(grps %>% select(Code,Name), by = c("spp"="Code")) %>%
      select(year, box_id, Name, fleet, weight_mton) %>%
      mutate(weight_mton = replace_na(weight_mton, 0)) %>%
      group_by(year, fleet, Name) %>%
      summarise(mt = sum(weight_mton)) %>% # sum up across boxes
      ungroup() %>%
      filter(year > (max(year)-5)) %>%
      group_by(fleet, Name) %>%
      summarize(mt = mean(mt)) %>% # mean of last 5 years
      ungroup()
    
    # add missing combinations
    dummy_df <- expand.grid("fleet" = unique(catch_dat$fleet), "Name" = unique(catch_dat$Name))
    catch_dat <- catch_dat %>%
      full_join(dummy_df) %>%
      mutate(mt = replace_na(mt, 0))
    
    # get relative catch
    catch_dat <- catch_dat %>%
      group_by(Name) %>%
      mutate(tot_mt = sum(mt, na.rm = T)) %>%
      ungroup() %>%
      mutate(prop = mt / tot_mt)
    
    # deal with code names
    catch_dat <- catch_dat %>% 
      rowwise() %>%
      mutate(fleet = rewrite_codes(as.character(fleet))) %>%
      ungroup()
    
    # stack
    catch_dat <- catch_dat %>% mutate(Type = "data")
    catch_nc_end <- catch_nc_end %>% mutate(Type = "model")
    catch_diff <- rbind(catch_dat, catch_nc_end)
    catch_diff <- catch_diff %>%
      left_join(fleet_key %>% select(Code, Name), by = c("fleet" = "Code"))
    
    # make a directory
    dir.create(paste0(plotdir, "/compositions/species/"), recursive = T)
    
    # make a plot per species
    for(i in 1:length(all_fg)){
      this_fg <- all_fg[i]
      
      if(this_fg %in% unique(catch_dat$Name)){
        
        print(paste("Doing", this_fg))
        
        tmp <- catch_diff %>%
          filter(Name.x == this_fg) %>%
          filter(fleet %in% to_keep) %>%
          filter(prop > 0)
        
        colourCount <- length(unique(tmp$Name.y)) 
        getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
        
        tmp %>%
          ggplot(aes(x = Type, y = prop, fill = Name.y))+
          geom_bar(stat = "identity", position = "stack")+
          scale_fill_manual(values = getPalette(colourCount))+
          theme_bw()+
          facet_wrap(~Name.x)
        ggsave(paste0(plotdir, "/compositions/species/data_vs_", new_run, this_fg, ".png"), width = 10, height = 8)
        
      }
      
    }
    
  } else { # if you want to slice this by fleet
    
    if(!fleet_struc)stop("You want to plot fleet catch in space but you are not keeping the fleet structure")
    
    # average of end of the run
    catch_nc_end <- catch_nc_new %>%
      filter(ts > (max(ts)-5)) %>%
      filter(box_id < 92) %>% # keep AK only
      group_by(fleet, Name) %>%
      summarize(mt = mean(mt)) %>%
      ungroup() %>%
      group_by(fleet) %>% # group by fleets
      mutate(tot_mt = sum(mt)) %>%
      ungroup() %>%
      mutate(prop = mt / tot_mt)
    
    # no space - fleets only
    catch_dat <- catch_data %>%
      ungroup() %>%
      left_join(grps %>% select(Code,Name), by = c("spp"="Code")) %>%
      select(year, box_id, Name, fleet, weight_mton) %>%
      mutate(weight_mton = replace_na(weight_mton, 0)) %>%
      group_by(year, fleet, Name) %>%
      summarise(mt = sum(weight_mton)) %>% # sum up across boxes
      ungroup() %>%
      filter(year > (max(year)-5)) %>%
      group_by(fleet, Name) %>%
      summarize(mt = mean(mt)) %>%
      ungroup()
    
    # add missing combinations
    dummy_df <- expand.grid("fleet" = unique(catch_dat$fleet), "Name" = unique(catch_dat$Name))
    catch_dat <- catch_dat %>%
      full_join(dummy_df) %>%
      mutate(mt = replace_na(mt, 0))
    
    # get relative catch
    catch_dat <- catch_dat %>%
      group_by(fleet) %>%
      mutate(tot_mt = sum(mt, na.rm = T)) %>%
      ungroup() %>%
      mutate(prop = mt / tot_mt)
    
    # deal with codes for the fleets being different
    catch_dat <- catch_dat %>% 
      rowwise() %>%
      mutate(fleet = rewrite_codes(as.character(fleet))) %>%
      ungroup()
    
    # stack
    catch_dat <- catch_dat %>% mutate(Type = "data")
    catch_nc_end <- catch_nc_end %>% mutate(Type = "model")
    catch_diff <- rbind(catch_dat, catch_nc_end)
    catch_diff <- catch_diff %>%
      left_join(fleet_key %>% select(Code, Name), by = c("fleet" = "Code"))
    
    # make a directory
    dir.create(paste0(plotdir, "/compositions/fleet/"), recursive = T)
    
    # make a plot per species
    for(i in 1:length(to_keep)){
      this_fleet <- to_keep[i]
      
      print(paste("Doing", this_fleet))
      
      catch_diff %>%
        filter(fleet == this_fleet) %>%
        filter(Name.x %in% all_fg) %>% # do vertebrates only
        filter(prop > 0) %>%
        ggplot(aes(x = Type, y = prop, fill = Name.x))+
        geom_bar(stat = "identity", position = "stack")+
        theme_bw()+
        facet_wrap(~Name.y)
      ggsave(paste0(plotdir, "/compositions/fleet/data_vs_", new_run, this_fleet, ".png"), width = 10, height = 8)
      
      
    }
    
  }
  
}



# Calibration functions ---------------------------------------------------


calibrate_mfc_total <- function(harvest_old, harvest_new, scalars){
  
  harvest_tab <- readLines(harvest_old)
  
  for(i in 1:nrow(scalars)){
    
    this_code <- scalars[i,]$Code
    this_scalar <- scalars[i,]$scalar
    
    # old mFC
    mfc_old <- harvest_tab[grep(paste0("mFC_", this_code, " 33"), harvest_tab)+1]
    mfc_old_vec <- as.numeric(unlist(strsplit(mfc_old, " ")))
    
    # new vec
    mfc_new_vec <- mfc_old_vec * this_scalar
    mfc_new <- paste(as.character(mfc_new_vec), collapse = " ")
    
    # replace relevant line
    harvest_tab[grep(paste0("mFC_", this_code, " 33"), harvest_tab)+1] <- mfc_new
  }
  
  # write out
  writeLines(harvest_tab, con = harvest_new)
  
}


calibrate_MPAYYY <- function(harvest_old, harvest_new, scalars){
  
  harvest_tab <- readLines(harvest_old)
  
  for(i in 1:length(to_keep)){
    
    this_code <- to_keep[i]
    this_scalar_vec <- scalars %>%
      filter(fleet == this_code) %>%
      arrange(box_id) %>%
      pull(scalar)
    
    # old mFC
    MPAYYY_old <- harvest_tab[grep(paste0("MPA", this_code, " 109"), harvest_tab)+1]
    MPAYYY_old_vec <- as.numeric(unlist(strsplit(MPAYYY_old, " ")))
    
    # fix fleets that are fishing in Canada and should not
    if(MPAYYY_old_vec[93] != 0) MPAYYY_old_vec[93] <- 0
    
    # new vec
    MPAYYY_new_vec <- MPAYYY_old_vec * this_scalar_vec
    MPAYYY_new <- paste(as.character(MPAYYY_new_vec), collapse = " ")
    
    # replace relevant line
    harvest_tab[grep(paste0("MPA", this_code, " 109"), harvest_tab)+1] <- MPAYYY_new
  }
  
  # write out
  writeLines(harvest_tab, con = harvest_new)
  
}

