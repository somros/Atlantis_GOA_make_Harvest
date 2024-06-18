
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
