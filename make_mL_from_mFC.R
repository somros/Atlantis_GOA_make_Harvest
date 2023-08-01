# Alberto Rovellini
# 3/3/2023
# Script to change mFC values to mL values
# The need for this arose to avoid applying too heavy of an mFC during model spinup and before applying catch.ts to the GOA
# Higher mFC helps with age structures, but relying on high mFC sets up for having issues when F is lifted and catch.ts is applied
# Transfer some (or all) of that mFC to mL on adults
# Importantly, both mFC and mL to aid with age structure equilibrium is "cheating". The underlying issue is a lack of predation mortality on larger individuals (gape size refuge, underestimation in the diets, etc.), especially but not only for high trophic levels
# 

# Read in mFC vector
# this is the fully-selected 1/4 FMSY mFC vector created in the atHarvest code, pre-tuning
this_mFC <- readLines('C:/Users/Alberto Rovellini/Documents/GOA/F/at_harvestR_GOA/mFC_tuning/mFC_to_mL.txt')

this_mL <- readLines('C:/Users/Alberto Rovellini/Documents/GOA/F/at_harvestR_GOA/mFC_tuning/mL_template.txt')

# groups that are being affected by mFC
groups_mFC <- this_mFC[grep('mFC', this_mFC)]
groups_mFC <- gsub(' 33.*', '', gsub('mFC_', '', groups_mFC))

# and which verts groups do we have in mL
groups_mL <- this_mL[grep('_mL', this_mL)]
groups_mL <- gsub('_mL.*', '', groups_mL)

# for each group in groups_mL, go dig out the value of mFC, modify it, and and write it out 

multiplier <- 0.5
do_all <- TRUE
do_adults <- TRUE
do_juveniles <- TRUE
new_mL_file <- paste('C:/Users/Alberto Rovellini/Documents/GOA/F/at_harvestR_GOA/mFC_tuning/','mL_vec_BY_', multiplier, '_all_', do_all, '_AJ.txt', sep = '')
file.create(new_mL_file)

# optional list of groups

to_do <- c("COD", "FHS", "REX", "FFS", "FFD", "SKL", "SKB", "SKO", "SBF", 
           "RFS", "RFP", "THO", "DFS", "DFD", "SCU", "CAP", "SAN", "FOS",
           "SCH", "SCM", "SCO", "SSO", "SPI", "HAK")

for( i in 1:length(groups_mL)){
  
  print(paste('Doing ', groups_mL[i]))
  
  this_group <- groups_mL[i]
  
  if(do_all){
    
    this_group <- groups_mL[i]
    
    this_mFC_value <- as.numeric(unlist(strsplit(this_mFC[grep(this_group, this_mFC)+1], ' '))[1])
    
    new_mL <- juv_mL <- adult_mL <- this_mFC_value * multiplier
    
    if(!do_juveniles){
      juv_mL <- 0
    }
    
    if(!do_adults){
      adult_mL <- 0
    }
    
    XXX_mL_line <- this_mL[grep(this_group, this_mL)]
    new_mL_line <- paste(juv_mL, adult_mL, sep = ' ')
    
    cat(XXX_mL_line, file=new_mL_file, append=TRUE,'\n\n')
    cat(new_mL_line, file=new_mL_file, append=TRUE, '\n\n\n')
    
  } else {
    
    this_mFC_value <- as.numeric(unlist(strsplit(this_mFC[grep(this_group, this_mFC)+1], ' '))[1])
    
    new_mL <- juv_mL <- adult_mL <- this_mFC_value * multiplier
    
    if(!do_juveniles){
      juv_mL <- 0
    }
    
    if(!do_adults){
      adult_mL <- 0
    }
    
    XXX_mL_line <- this_mL[grep(this_group, this_mL)]
    
    if(this_group %in% to_do){
      new_mL_line <- paste(juv_mL, adult_mL, sep = ' ')
    } else {
      new_mL_line <- paste(0, 0, sep = ' ')
    }
    
    cat(XXX_mL_line, file=new_mL_file, append=TRUE,'\n\n')
    cat(new_mL_line, file=new_mL_file, append=TRUE, '\n\n\n')
    
  }
  
}

# What's the plan for invertebrates? It should be marginal but expect differences if we remove F
