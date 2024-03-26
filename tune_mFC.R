# Calibrate mFC values
library(tidyverse)
library(data.table)

# REad data ---------------------------------------------------------------
# groups
fg <- read.csv('data/GOA_Groups.csv', header = T)
fg_codes <- fg %>% pull(Code)

# fleets
fleets <- read.csv('data/GOA_fisheries.csv')
fleet_codes <- fleets %>% pull(Code)
fleet_names <- fleets %>% pull(Name)
n_fleets <- length(fleet_codes)

# pull the catch by fishery file from the most recent output folder
run <- 1516

tune_mFC <- function(run){
  outdir <- paste0('C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_', run)
  outfiles <- list.files(outdir, full.names = T)
  
  biomfile <- outfiles[grepl('_testBiomIndx', outfiles)]
  cbf_file <- outfiles[grepl('_testCatchPerFishery', outfiles)]
  harvestPrm_file <- outfiles[grepl('GOA_harvest_background.prm', outfiles)]
  
  # initial biomass from Atlantis output
  # this is tied to initial conditions, still tie it to the correct file
  biom <- read.table(biomfile, sep = ' ', header = T) %>%
    filter(Time == 0) %>%
    pivot_longer(-Time, names_to = 'Code', values_to = 'Biomass') %>%
    filter(Code %in% fg$Code) %>%
    select(-Time)
  
  catchout <- read.table(cbf_file, sep = ' ', header = T)
  
  # only keep time step 1
  catchout <- catchout %>%
    filter(Time == 365) %>%
    distinct() %>%
    select(-Time)
  
  # fishery data
  catchdata_tmp <- read.csv('data/all_fisheries_goa_25FMSY_OY.csv')
  catchdata_tmp <- catchdata_tmp %>% select(-c(year,box,month,tot_catch_mt))
  
  # set this in the same format as catchout
  catchdata_tmp <- catchdata_tmp %>%
    filter(atlantis_fg %in% names(catchout[,-1])) %>% # keep groups that got fished %>%
    left_join(fleets %>% select(Code, Name), by = c('gear_name' = 'Name')) %>%
    select(-gear_name) %>%
    arrange(factor(atlantis_fg, levels=names(catchout[,-1])))
  
  catchdata <- catchdata_tmp %>%
    pivot_wider(names_from = atlantis_fg, values_from = catch_tons) %>%
    rename(Fishery = Code)
  
  # get ratios, and transpose 
  expToObs <- t(catchdata[,-1]/catchout[,-1]) %>% data.frame()
  colnames(expToObs) <- catchdata %>% pull(Fishery)
  expToObs <- expToObs %>% mutate(Code = row.names(.))
  expToObs <- expToObs %>% 
    mutate_all(~ifelse(is.nan(.), 0, .))
  
  # split this into groups to get a range of how bad it was at first
  tier3 <- c('POL','COD','ATF','SBF','POP','REX','FFS','FHS','RFS','FFD','HAL')
  tier_45 <- c('DOG','SKL','SKB','SKO','RFP','RFD','THO')
  migrating <- c('HAK','SPI','SCH','SCM','SCO','SPI','SSO')
  topverts <- c('WHT','WHB','DOL','KWT','KWR','SSL','PIN','BSF','BSI','BDF','BDI','SHD','WHH','WHG','SHP')
  forage <- c('SAN','CAP','FOS','EUL','HER')
  other <- c(setdiff(expToObs$Code, c(tier3, tier_45, migrating, topverts, forage)))

  expToObs %>% filter(Code %in% tier3) %>% pull(background) %>% range()
  expToObs %>% filter(Code %in% tier_45) %>% pull(background) %>% range()
  expToObs %>% filter(Code %in% migrating) %>% pull(background) %>% range()
  expToObs %>% filter(Code %in% topverts) %>% pull(background) %>% range()
  expToObs %>% filter(Code %in% forage) %>% pull(background) %>% range()
  expToObs %>% filter(Code %in% other) %>% pull(background) %>% range()
  
  # now for the tricky part, alter the harvest.prm
  harvestPrm <- readLines(harvestPrm_file)
  
  mFC_list <- list()
  for(i in 1:length(fg_codes)) {
    
    grp <- fg_codes[i]
    idx <- grep(paste0('mFC_', grp, ' 33'), harvestPrm)
    
    parname <- harvestPrm[idx]
    parvals <- harvestPrm[idx + 1] %>% 
      strsplit(' ') %>% 
      unlist() %>% 
      as.numeric() %>% 
      matrix(nrow=1) %>%
      data.frame() %>%
      select_if(~ !any(is.na(.))) %>%
      set_names(fleet_codes)
    
    this_mFC <- cbind(parname, parvals)
    mFC_list[[i]] <- this_mFC
  }
  
  mFC_frame_old <- mFC_frame_new <- rbindlist(mFC_list)
  
  # function to modify mFC based on expToObs scalars
  for(i in 1:nrow(mFC_frame_old)){
    this_fg <- gsub(' 33.*', '', gsub('mFC_', '', mFC_frame_old[i,1]))
    
    if(length(grep(this_fg, expToObs$Code)) == 1){
      mFC_frame_new[i,] <- cbind(mFC_frame_old[i,1],
                                 matrix(as.numeric(mFC_frame_old[i,-1]) * as.numeric(expToObs[grep(this_fg, expToObs$Code),-ncol(expToObs)]), nrow=1))
    }
  }
  
  # turn mFC_frame_new back to section of prm file and write out
  newfile <- paste0('mFC_tuning/mFC_from_',run, '.prm')
  file.create(newfile)
  
  for(i in 1:nrow(mFC_frame_new)){
    
    cat(mFC_frame_new[i,]$parname, file=newfile, append=TRUE,'\n')
    cat(as.numeric(mFC_frame_new[i,-1]), file=newfile, append=TRUE, '\n')
  }
  
}

tune_mFC(run)

