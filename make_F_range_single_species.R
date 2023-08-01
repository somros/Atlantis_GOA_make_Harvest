# Alberto Rovellini
# 7/31/2023
# This code takes the harvest.prm file from a simulation and manipulates values of mFC
# the intent is to create a number of harvest.prm files with varying levels of mFC
# This has several applications, such as:
# 1. Do single-species tests of levels of fishing, with mFC for all other species being fixed to a baseline
# 2. Set up harvest.prm files for simulations where all fishing varies as a multiplier of a reference F (such as status quo ot a reference point)
# 3. Manpulate the mFC matrix by row (species) or column (fleet) (not applicable for now as we still have one F for all fleets)

library(tidyverse)

# read groups.csv file
grps <- read.csv('data/GOA_Groups.csv')
# filter by those that are impacted
grps_all <- grps %>% pull(Code)
grps_imp <- grps %>% filter(IsImpacted == 1) %>% pull(Code)
grps_fish <- grps %>% filter(GroupType %in% c('FISH','SHARK')) %>% pull(Code)

# how many fleets do you have?
n_fleets <- 33

# read harvest.prm
harvest_prm_file <- '../../Parametrization/output_files/data/out_1328/GOA_harvest_background.prm'
harvest_prm <- readLines(harvest_prm_file)

# construct data frame with nrow = n species and ncol = n fleets, each entry is mFC
mfc_list <- list()

for(i in 1:length(grps_all)){
 
  this_code <- grps_all[i]
  this_handle <- paste0('mFC_', this_code, ' ', n_fleets)
  this_mFC <- data.frame(matrix(as.numeric(unlist(strsplit(harvest_prm[grep(this_handle, this_prm)+1], ' '))), nrow = 1))
  colnames(this_mFC) <- 0:(n_fleets-1)
  mfc_list[[i]] <- cbind('name' = this_handle, this_mFC)
}

mfc_frame <- bind_rows(mfc_list)


# Single-species F --------------------------------------------------------
# make mFC vectors based on following range of F:
# F = 0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.25, 1.5, 2

# Formula for mFC is mfc = 1-exp(-F / 365)

f_vals <- c(0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.00, 1.25, 1.50, 2.00)
mfc_vals <- 1-exp(-f_vals / 365)

# create a file per value of F per fish species, leaving the rest as is. This is 408 files (and as many runs)
# If we want to do it also for the all impacted groups (mammals, birds, and invertebrates), it will be 720 runs
dir.create('single_species_F')

for (i in 1:length(grps_fish)){
  this_code <- grps_fish[i]
  print(paste('Doing', this_code, sep = ' '))
  
  this_handle <- paste0('mFC_', this_code, ' ', n_fleets)
  newdir <- paste0('single_species_F/', this_code)
  dir.create(newdir)
  
  for(j in 1:length(mfc_vals)){
    
    this_mfc_val <- mfc_vals[j]
    print(paste('Doing', f_vals[j], sep = ' '))
    
    # create a file label
    file_label <- gsub('\\.', '', as.character(f_vals[j]))
    
    # create new file
    newfile <- paste0(newdir, '/', 'mFC_', file_label, '.prm')
    file.create(newfile)
    
    # keep the front and the back of the prm file as is, and only manipulate the 
    part1 <- paste0(harvest_prm[1:grep(this_handle, harvest_prm)], collapse = '\n')
    part2 <- paste(c(as.character(mfc_vals[j]), 
               unlist(strsplit(harvest_prm[grep(this_handle, harvest_prm)+1],' '))[-1]), collapse = ' ')
    part3 <- paste0(harvest_prm[(grep(this_handle, harvest_prm)+2):length(harvest_prm)], collapse = '\n')
    
    # concatenate
    cat(part1, file=newfile, append=TRUE,'\n')
    cat(part2, file=newfile, append=TRUE, '\n') 
    cat(part3, file=newfile, append=TRUE, '\n') 
  }
}




