
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

build_catch_output <- function(catch_nc, bio_nc, age_struc, relative, run){
  
  # catch_nc = catch_nc_file_mpa
  # bio_nc = bio_nc_file_mpa
  # age_struc = F
  # relative = T
  # run = 1553
  
  # catch
  this_tidync <- tidync(catch_nc)
  this_nc <- ncdf4::nc_open(catch_nc)
  #bio
  this_tidync_bio <- tidync(bio_nc)
  this_nc_bio <- ncdf4::nc_open(bio_nc)
  
  catch_nc_ls <- list()
  
  for(n in 1:length(all_fg)){
    fg <- all_fg[n] # this needs to use the "Name" to pull from the NC file
    
    #Extract from the output .nc file the appropriate catch time series variables
    catch_vars <- this_tidync %>%
      activate("D1,D0") %>%
      hyper_vars() %>% # all variables in the .nc file active grid
      filter(grepl("_Catch",name)) %>% # filter for reserve N
      filter(grepl(fg,name)) # filter for specific functional group
    
    #Extract from the output .nc file the appropriate reserve N time series variables
    resN_vars <- hyper_vars(this_tidync_bio) %>% # all variables in the .nc file active grid
      filter(grepl("_ResN",name)) %>% # filter for reserve N
      filter(grepl(fg,name)) # filter for specific functional group
    
    #Extract from the output .nc file the appropriate structural N time series variables
    strucN_vars <- hyper_vars(this_tidync_bio) %>% # all variables in the .nc file active grid
      filter(grepl("_StructN",name)) %>% # filter for structural N
      filter(grepl(fg,name)) # filter for specific functional group
    
    # Get numbers by box
    # going to need these for weighting over depth layers (which we need to average over because catch is by box and not by depth)
    abun_vars <- hyper_vars(this_tidync_bio) %>% # all variables in the .nc file active grid
      filter(grepl("_Nums",name)) %>% # filter for abundance variables
      filter(grepl(fg,name)) # filter for specific functional group
    
    if(nrow(resN_vars)==0) {return("no data.")}
    else {
      # # Actually pull the data from the .nc
      # start from bio information
      resN <- purrr::map(resN_vars$name,ncdf4::ncvar_get,nc=this_nc_bio) 
      strucN <- purrr::map(strucN_vars$name,ncdf4::ncvar_get,nc=this_nc_bio)
      nums <-purrr::map(abun_vars$name,ncdf4::ncvar_get,nc=this_nc_bio) #numbers by age group,box,layer,time
      
      # number of time steps
      t_steps <- dim(resN[[1]])[3]
      
      # sum res and struct and then take averages over the water column, ending up with a value per box per time step
      weighted_avg_list <- vector("list", length = length(resN)) # Initialize the list to store the weighted averages
      
      # Perform the operation (A+B) and then calculate the weighted average using C
      weighted_mean_func <- function(x, w) {
        sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
      }
      
      weighted_waa <- list()
      for(i in 1:length(resN)) {
        # Calculate the sum of A and B
        sum_ab <- (resN[[i]] + strucN[[i]]) * 20 * 5.7 / 1000000000 # from from mgN to tons
        this_nums <- nums[[i]] # numbers at this age over the spatial structure
        
        # Initialize an array to store the results
        this_weighted_waa <- array(dim = c(109, t_steps))
        
        # Calculate weighted mean across the first dimension
        for(j in 1:109) {
          for(k in 1:t_steps) {
            # Extract the slice for the current i, j across all 7 values
            sum_ab_slice <- sum_ab[, j, k]
            nums_slice <- this_nums[, j, k]
            
            # Calculate the weighted mean for this slice
            this_weighted_waa[j, k] <- weighted_mean_func(sum_ab_slice, nums_slice)
          }
        }
        
        # turn to df and pivot
        this_weighted_waa <- this_weighted_waa %>%
          as.data.frame() %>%
          mutate(box_id = 0:108) %>%
          pivot_longer(-box_id, values_to = "mt", names_to = "ts") %>%
          mutate(ts = gsub("V","",ts)) %>%
          mutate(ts = as.numeric(ts)) %>%
          mutate(ts = ts - 1) %>%
          mutate(age = i - 1) # count from 0
        
        weighted_waa[[i]] <- this_weighted_waa
      }
      
      # turn it ito a df
      weighted_waa_df <- bind_rows(weighted_waa)
      
      # drop t 0 and then take every 5th 
      weighted_waa_df <- weighted_waa_df %>%
        filter(ts > 0) %>%
        filter(ts %in% seq(5,max(ts),5)) %>%
        mutate(ts = ts / 5)
      
      # now let's extract the catch for this fg
      # here we can collapse the depth layers but need to keep the boxes
      catch <- purrr::map(catch_vars$name,ncdf4::ncvar_get,nc=this_nc) 
      
      # now turn to data frame
      catch_box_ls <- list()
      
      # Loop over each matrix in the collapsed list to fill the data frame
      for(i in 1:length(catch)) {
        # Get the current matrix
        mat <- catch[[i]]
        
        # Convert matrix to a data frame
        mat_df <- as.data.frame((mat))
        
        # add box_id
        mat_df <- mat_df %>% mutate(box_id = 0:108)
        
        # reshape
        mat_df_long <- mat_df %>%
          pivot_longer(-box_id, names_to = "ts", values_to = "nums")
        
        # turn ts column to integer
        mat_df_long <- mat_df_long %>%
          mutate(ts = gsub("V","",ts)) %>%
          mutate(ts = as.numeric(ts)) %>%
          mutate(ts = ts - 1) # start numbering ts from 0
        
        # add age
        mat_df_long <- mat_df_long %>%
          mutate(age = i-1) # number from 0 for consistency with age at selex and age mat
        
        catch_box_ls[[i]] <- mat_df_long
        
      }
      
      catch_box_df <- bind_rows(catch_box_ls) # this is in numbers... so we'd need to multiply this by WAA, pulling that in
      
      # drop ts = 0
      catch_box_df <- catch_box_df %>%
        filter(ts > 0)
      
      # now join and multiply to get mt per box
      catch_box_df <- catch_box_df %>%
        left_join(weighted_waa_df) %>%
        mutate(mt_tot = nums * mt) %>%
        mutate(Name = fg)
      
    }
    
    catch_nc_ls[[n]] <- catch_box_df
  }
  
  catch_nc_df <- bind_rows(catch_nc_ls)
  
  # turn NaN to NA
  catch_nc_df$mt[is.nan(catch_nc_df$mt)] <- NA
  catch_nc_df$mt_tot[is.nan(catch_nc_df$mt_tot)] <- NA
  
  # get the equilibrium for the last 5 years
  catch_nc_eq <- catch_nc_df %>%
    filter(ts > (max(ts)-5)) %>%
    group_by(box_id, Name, age) %>% # get average of last 5 years of the run
    summarise(mt_tot = mean(mt_tot, na.rm = T)) %>%
    ungroup()
  
  # make a spatial version
  catch_nc_spatial <- goa_sf %>%
    select(box_id) %>%
    full_join(catch_nc_eq %>%
                select(box_id, Name, age, mt_tot),
              by = "box_id")
  
  # if age_struc, keep ages, otherwise sum across them
  if(!age_struc){
    catch_nc_spatial <- catch_nc_spatial %>%
      group_by(box_id, Name) %>%
      summarize(mt_tot = sum(mt_tot, na.rm = T)) %>%
      ungroup()
    
    # if relative catch is wanted, rescale and get proportions
    if(relative){
      catch_nc_spatial <- catch_nc_spatial %>%
        group_by(Name) %>%
        mutate(mt_goa = sum(mt_tot, na.rm = T)) %>%
        ungroup() %>%
        mutate(prop = mt_tot / mt_goa)
    }
    
  } else {
    
    # if relative catch is wanted, rescale and get proportions
    if(relative){
      catch_nc_spatial <- catch_nc_spatial %>%
        group_by(Name, age) %>%
        mutate(mt_goa = sum(mt_tot, na.rm = T)) %>%
        ungroup() %>%
        mutate(prop = mt_tot / mt_goa)
      
    }
    
  }
  
  # add run information
  catch_nc_spatial <- catch_nc_spatial %>%
    mutate(run = run)
  
  return(catch_nc_spatial)
  
}

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
