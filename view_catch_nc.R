# this script extracts and views catch in space
# allow for comparisons between runs (e.g., Base no MPA, MPA fully open, in the future partially closed MPAs, etc)

library(tidync)
library(ncdf4)
library(tidyverse)
library(sf)
library(rbgm)
library(viridis)

select <- dplyr::select

# nc files
# base
catch_nc_file_base <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1517/outputGOA01517_testCATCH.nc"
bio_nc_file_base <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1517/outputGOA01517_test.nc"
# mpa
catch_nc_file_mpa <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1555/outputGOA01555_testCATCH.nc"
bio_nc_file_mpa <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1555/outputGOA01555_test.nc"
# groups
grps <- read.csv("data/GOA_Groups.csv")
all_fg <- grps %>% filter(GroupType %in% c("FISH", "SHARK")) %>% pull(Name)

# geometry
goa_bgm <- read_bgm("data/GOA_WGS84_V4_final.bgm")
goa_sf <- goa_bgm %>% box_sf()

# writing a function that:
# reads bio and catch nc files
# calculates mean weight at age per species per boxm weighted by the number of individuals in the box / layer
# extracts catch output per cell (in numbers) and converts it to tons
# optionally, collapses across age classes for total catch
# optionally, it turns it into a relative catch index to examine spatial distributions of catch

build_catch_output <- function(catch_nc, bio_nc, age_struc, relative, run){
  
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
    resN_vars <- hyper_vars(out_bio) %>% # all variables in the .nc file active grid
      filter(grepl("_ResN",name)) %>% # filter for reserve N
      filter(grepl(fg,name)) # filter for specific functional group
    
    #Extract from the output .nc file the appropriate structural N time series variables
    strucN_vars <- hyper_vars(out_bio) %>% # all variables in the .nc file active grid
      filter(grepl("_StructN",name)) %>% # filter for structural N
      filter(grepl(fg,name)) # filter for specific functional group
    
    # Get numbers by box
    # going to need these for weighting over depth layers (which we need to average over because catch is by box and not by depth)
    abun_vars <- hyper_vars(out_bio) %>% # all variables in the .nc file active grid
      filter(grepl("_Nums",name)) %>% # filter for abundance variables
      filter(grepl(fg,name)) # filter for specific functional group
    
    
    if(nrow(resN_vars)==0) {return("no data.")}
    else {
      # # Actually pull the data from the .nc
      # start from bio information
      resN <- purrr::map(resN_vars$name,ncdf4::ncvar_get,nc=this_nc_bio) 
      strucN <- purrr::map(strucN_vars$name,ncdf4::ncvar_get,nc=this_nc_bio)
      nums <-purrr::map(abun_vars$name,ncdf4::ncvar_get,nc=this_nc_bio) #numbers by age group,box,layer,time
      
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
        this_weighted_waa <- array(dim = c(109, 251))
        
        # Calculate weighted mean across the first dimension
        for(j in 1:109) {
          for(k in 1:251) {
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
        filter(ts %in% seq(5,250,5)) %>%
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
  
  catch_nc <- bind_rows(catch_nc_ls)
  
  # turn NaN to NA
  catch_nc$mt[is.nan(catch_nc$mt)] <- NA
  catch_nc$mt_tot[is.nan(catch_nc$mt_tot)] <- NA
  
  # get the equilibrium for the last 5 years
  catch_nc_eq <- catch_nc %>%
    filter(ts > 44) %>%
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

catch_nc_base <- build_catch_output(catch_nc = catch_nc_file_base, 
                                    bio_nc = bio_nc_file_base, 
                                    age_struc = F,
                                    relative = T,
                                    run = 1517)

catch_nc_mpa <- build_catch_output(catch_nc = catch_nc_file_mpa, 
                                   bio_nc = bio_nc_file_mpa,
                                   age_struc = F,
                                   relative = T,
                                   run = 1555)

# get residuals
catch_diff <- catch_nc_base %>%
  left_join(catch_nc_mpa %>% st_set_geometry(NULL), by = c("box_id", "Name")) %>%
  mutate(residual_mt = mt_tot.x - mt_tot.y,
         residual_prop = prop.x - prop.y)
  
# view
catch_diff %>%
  filter(Name != "Pacific_hake") %>%
  ggplot()+
  geom_sf(aes(fill = residual_prop))+
  scale_fill_viridis()+
  theme_bw()+
  facet_wrap(~Name)

# When the MPA cells are fully open, the spatial distribution seems to be similar
# However, the total catch is much lower, so why is that the case?

# catch seems much lower in the file with MPAs
# This can also be seen in the CSV files
# catch_flat_base <- read.table("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1517/outputGOA01517_testCatch.txt", sep = " ", header = T)
# catch_flat_mpa <- read.csv("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1555/outputGOA01555_testCatch.txt", sep = " ", header = T)
# 
# catch_flat_base_long <- catch_flat_base %>%
#   pivot_longer(-Time, names_to = "Code", values_to = "mt") %>%
#   mutate(run = "base")
# 
# catch_flat_mpa_long <- catch_flat_mpa %>%
#   pivot_longer(-Time, names_to = "Code", values_to = "mt") %>%
#   mutate(run = "mpa")
# 
# catch_flat <- rbind(catch_flat_base_long, catch_flat_mpa_long)
# 
# # view
# catch_flat %>%
#   filter(Code %in% (grps %>% filter(Name %in% all_fg) %>% pull(Code)))  %>%
#   ggplot(aes(x = Time, y = mt, color = run))+
#   geom_line()+
#   theme_bw()+
#   facet_wrap(~Code, scales = "free")
