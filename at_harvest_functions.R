# this function makes the mFC vectors 
# there is one vector per species, with as many entries as fleets
# it needs time series of catch, that it averages from
# in the simple case of one background F at 25% MSY, we can get catch based on B1990 and FMSY 
# Previous conversations hinted that FMSY in 1990 for many GOA stocks is probably similar to present days, so just use values from Martin
# Make time series of catch based on this. Homogeneous by box (it all gets lumped anyway)
# Make it for one year only as this is a dummy and supposed to be constant
# When we will prepare actual mFc values we will need to map catch.ts to this somehow

make_mfc <- function(all.goa.fisheries, fleet.list, group.list, goa.biomass){
  
  all.gears <- all.goa.fisheries %>% distinct(gear_name)
  
  # This makes a df of mean catch per group per fleet
  # it requires time series of fished biomass by fleet to work
  # Take a time period
  # get average catch over said time period per fleet per species per year (tot model domain)
  # then average this catch over the selected years
  tot.fishing.average <- all.goa.fisheries %>% 
    #filter(year %in% c(2009:2019)) %>% # may need this for real catch
    filter(catch_tons!=0) %>% # but this drops the fisheries that did not catch anything, only to add them in later. It also drops everything that is not a vertebrate
    group_by(year,atlantis_fg,gear_name) %>% #sum across for each box and month
    summarise(tot_catch_mt=sum(catch_tons)) %>% 
    ungroup %>% 
    group_by(atlantis_fg,gear_name) %>% #sum across boxes for each year
    summarise(mean_catch=mean(tot_catch_mt)) #get average across years
  
  # this sums across gears
  # so equivalent to avg tot catch per species per year (i.e. out catch.ts files)
  tot.fishing.initial <- tot.fishing.average %>% 
    group_by(atlantis_fg) %>% 
    summarise(catch_tons=sum(mean_catch))
  
  # write_csv(tot.fishing.initial, "data/initial_mfc_catch.csv")
  # 
  # write_csv(tot.fishing.average, "data/mfc_average_catch.csv")
  
  # by species by fleet, calculated based on average catch over a period
  # The calculation to go from annual mu to daily F seems to be correct
  tot.fishing.mort <- tot.fishing.average %>%
    left_join(goa.biomass, by="atlantis_fg") %>% 
    mutate(annual_mu = mean_catch / biomass_mt, # exploitation rate
           annual_F = (-1 * log(1-annual_mu)), # F
           mFC = 1-exp(-annual_F / 365)) %>% # daily probability of being caught (see manual part II)
    ungroup()
  
  #update the fleets so that it includes all in the TS even if they are 0 in this point
  mfc.data <- fleet.list %>% 
    dplyr::select(gear_name) %>% 
    bind_rows(tot.fishing.mort) %>% 
    group_by(gear_name,atlantis_fg) %>% # not sure what this line and the one below do
    summarise(mFC = sum(mFC)) %>% 
    ungroup %>% 
    complete(atlantis_fg,gear_name, fill=list(mFC=0)) %>% 
    left_join(fleet.list, by="gear_name")
  
  fleet.order <- fleet.list %>% # these are the fleet codes
    pull(Code)
  
  write_csv(mfc.data, "data/goa_fishing_mortality.csv")
  
  species.list <- group.list %>% 
    pull(Code) # all groups, even detritus etc
  
  file.name <- "data/mfc_vector_SELEX.prm"
  file.create(file.name)
  
  # mfc.data %>% 
  #   distinct(gear_name) %>% 
  #   arrange(gear_name) %>% 
  #   write_csv("data/mfc_gears.csv")
  
  mfc.list <- list()
  
  for(eachentry in 1:length(species.list)) {
    
    eachspecies <- species.list[eachentry]
    print(eachspecies)
    
    this.data <- mfc.data %>% 
      filter(atlantis_fg==eachspecies) %>% 
      arrange(match(Code, fleet.order)) %>% 
      pull(mFC)
    
    print(eachspecies)
    
    if(length(this.data)==0){
      
      this.data <- rep(0,33) # because there are 33 fleets here - could make this dynamic or just change here
    }
    
    cat(paste("mFC_",eachspecies," 33",sep=""),file=file.name,append=TRUE, '\n')
    
    res.tibble <- this.data %>% 
      t%>% 
      as.data.frame() %>% 
      mutate(mfc=paste("mFC_",eachspecies," 59",sep="")) %>% 
      dplyr::select(mfc,everything())
    
    cat(paste0(this.data,collapse=" "), file=file.name,append=TRUE, '\n')
    
    mfc.list[[eachentry]] <- res.tibble
    
  }
  
  mfc.frame <- mfc.list %>% 
    bind_rows()
  
  return(mfc.frame)
}

# this function creates the harvest.prm file
# it will probably miss a couple of parameters
make_at_harvest <- function(grp.file, fsh.file, temp, bgm.file, cum.depths, harvest.file.name, run.type, this.mfc) {
  
  n.lay <- length(cum.depths) - 1 # number of depth layers
  df.grp <- read_csv(file = grp.file) # groups.csv
  df.fsh <- read_csv(file = fsh.file) # fleets
  tpe    <- read_csv(file = temp) #template - if we are missing parameters we need to add them here # also contains values - VERY IMPORTANT
  bgm  <- readLines(bgm.file) # bgm file
  n.box    <- as.numeric(gsub('nbox', '', grep("nbox", bgm, value = TRUE)))
  
  n.prms <- nrow(tpe) # how many params
  n.grp <- nrow(df.grp) # how many groups
  n.fsh <- nrow(df.fsh) # how many fisheries
  a.grp <- df.grp[which(df.grp$NumCohorts>1), ] # age structured groups
  hab   <- c(df.grp$Code[which(df.grp$IsCover == 1)], 'reef', 'flat', 'soft', 'canyon') # habitats
  n.ports <- 1 # this would be 0 for us for the dummy file and then however many ports / borroughs we use
  
  file.create(harvest.file.name)
  
  print(paste("Starting loop over rows",n.prms))
  
  for(p in 1:n.prms){
    
    print(p)
    
    
    if(!is.na(tpe$Legend[p])) {
      
      cat('\n', file=harvest.file.name,append=TRUE, '\n')
      cat(tpe$Legend[p], file=harvest.file.name,append=TRUE, '\n')
      
      
    } 
    
    if(tpe$Template[p]==1){
      
      if(is.na(tpe$TextPost[p])) tpe$TextPost[p] <- ""
      if(is.na(tpe$TextPre[p])) tpe$TextPre[p] <- ""
      if(!is.na(tpe$Legend[p])) tpe$Description[p] <- ""
      
      print("In template")
      
      # if(tpe$Flags[p] == '[TTL]'){
      #   ## Titles
      #   nch <- (nchar(tpe$Description[p])/2)
      #   cat('##', rep('~', nch), file=harvest.file.name,append=TRUE, '#\n')
      #   cat(paste0(tpe$Description[p]), file=harvest.file.name,append=TRUE, '\n')
      #   cat('##', rep('~', nch), file=harvest.file.name,append=TRUE, '#\n')
      # }
      
      if(tpe$Flags[p] == '[FLAGG]') {
        
        cat(paste0(tpe$TextPre[p], '\t', n.grp,'\t'), file=harvest.file.name,append=TRUE,'\n')
        cat(rep(tpe$Value[p], n.grp), file=harvest.file.name,append=TRUE, '\n')
        
      }
      
      if(tpe$Flags[p] == '[FLAG]'){
        ## Flags
        if(tpe$Vector[p]){
          ## Vectors
          if(tpe$Category[p] == 'Box'){
            for( f in 1 : n.fsh){
              
              print("Box")
              
              cat(paste0(tpe$TextPre[p], df.fsh$Code[f], tpe$TextPost[p],'\t' , n.box),file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], n.box), file=harvest.file.name,append=TRUE,'\n')
              
            }
          } else if(tpe$Category[p] == 'Fishery'){
            for( f in 1 : n.fsh){
              
              print("Fishery")
              
              if(f == 1) cat('## ', paste(df.grp$Code,'|'), file=harvest.file.name,append=TRUE,'\n')
              cat(paste0(tpe$TextPre[p], df.fsh$Code[f], tpe$TextPost[p],'\t' , n.grp, '\t## ',
                         tpe$Description[p],' ', df.fsh$Name[f]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], n.grp), file=harvest.file.name,append=TRUE,'\n')
              
            }
          } else if(tpe$Category[p] == 'GroupBox'){
            
            print("GroupBox")
            
            for(gr in 1 : n.grp){
              
              # if(gr == 1) cat('## ', paste(df.fsh$Code,'|'), file=harvest.file.name,append=TRUE,'\n')
              
              cat(paste0(tpe$TextPre[p], df.grp$Code[gr], tpe$TextPost[p],'\t', n.box), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], n.box), file=harvest.file.name,append=TRUE, '\n')
            }
            
            
          } else if(tpe$Category[p] == 'Group'){
            
            print("Group")
            
            for(gr in 1 : n.grp){
              
              # if(gr == 1) cat('## ', paste(df.fsh$Code,'|'), file=harvest.file.name,append=TRUE,'\n')
              
              cat(paste0(tpe$TextPre[p], df.grp$Code[gr], tpe$TextPost[p],'\t', n.fsh,'\t##', tpe$Description[p], df.grp$`Long Name`[gr]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], n.fsh), file=harvest.file.name,append=TRUE, '\n')
            }
            
            
          } else if(tpe$Category[p] == 'Cohort'){
            
            print("Cohort")
            
            for(agr in 1:nrow(a.grp)){
              cat(paste0(tpe$TextPre[p], a.grp$Code[agr], tpe$TextPost[p],'\t' , a.grp$NumCohorts[agr],'\t## ', tpe$Description[p], a.grp$Name[agr]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], a.grp$NumCohorts[agr]), file=harvest.file.name,append=TRUE, '\n')
            }
          } else if(tpe$Category[p] == 'SPSP'){
            for(gr in 1:n.grp){
              if(gr == 1) cat('## ', paste(df.grp$Code,'|'), file=harvest.file.name,append=TRUE,'\n')
              cat(paste0(tpe$TextPre[p], df.grp$Code[gr], tpe$TextPost[p],'\t' , n.grp,'\t## ', tpe$Description[p], ' ', df.grp$`Long Name`[gr]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], n.grp), file=harvest.file.name,append=TRUE, '\n')
            }
          } else if(tpe$Category[p] == 'Habitat'){
            for( f in 1 : n.fsh){
              if(f == 1) cat('## ', paste(hab,'|'), file=harvest.file.name,append=TRUE,'\n')
              cat(paste0(tpe$TextPre[p], df.fsh$Code[f], tpe$TextPost[p],'\t' , length(hab), '\t## ', tpe$Description[p], ' ', df.fsh$Name[f]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], length(hab)), file=harvest.file.name,append=TRUE, '\n')
            }
          } else if(tpe$Category[p] == 'FishFish'){
            for( f in 1 : n.fsh){
              if(f == 1) cat('## ', paste(df.fsh$Code,'|'), file=harvest.file.name,append=TRUE,'\n')
              cat(paste0(tpe$TextPre[p], df.fsh$Code[f], tpe$TextPost[p], '\t', n.fsh,'\t## ', tpe$Description[p], ' ', df.fsh$Name[f]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], n.fsh), file=harvest.file.name,append=TRUE, '\n')
            }
          }
          else if(tpe$Category[p] == 'Season'){
            for( f in 1 : n.fsh){
              cat(paste0(tpe$TextPre[p], df.fsh$Code[f], tpe$TextPost[p], '\t', 4,'\t## ', tpe$Description[p], ' ', df.fsh$Name[f]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], 4), file=harvest.file.name,append=TRUE, '\n')
            }
          } else if(tpe$Category[p] == 'FishSeason'){
            for( f in 1 : n.fsh){
              cat('## ', tpe$Description[p], ' ', df.fsh$Name[f], file=harvest.file.name,append=TRUE, '\n')
              for(s in 1 : 4){
                cat(paste0(tpe$TextPre[p], df.fsh$Code[f], tpe$TextPost[p], s, '\t', n.box), file=harvest.file.name,append=TRUE,'\n')
                cat(rep(tpe$Value[p], n.box), file=harvest.file.name,append=TRUE, '\n')
              }
              cat('\n')
            }
          } else if(tpe$Category[p] == 'FishVertical'){
            for( f in 1 : n.fsh){
              cat(paste0(tpe$TextPre[p], df.fsh$Code[f], tpe$TextPost[p], '\t', n.lay,'\t## ', tpe$Description[p], ' ', df.fsh$Name[f]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], n.lay), file=harvest.file.name,append=TRUE, '\n')
            }
          } else if(tpe$Category[p] == 'Port'){
            for( f in 1 : n.fsh){
              cat(paste0(tpe$TextPre[p], df.fsh$Code[f], tpe$TextPost[p], '\t', n.ports,'\t## ', tpe$Description[p], ' ', df.fsh$Name[f]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], n.ports), file=harvest.file.name,append=TRUE, '\n')
            }
          } else if(tpe$Category[p] == 'Value'){
            cat(paste0(tpe$TextPre[p], tpe$TextPost[p], '\t', tpe$Value[p], '\t## ', tpe$Description[p]), file=harvest.file.name,append=TRUE, '\n')
            cat(rep(1, tpe$Value[p]), file=harvest.file.name,append=TRUE, '\n')
          } else if(tpe$Category[p] == 'ValueGroup'){
            for(gr in 1 : n.grp){
              cat(paste0(tpe$TextPre[p], df.grp$Code[gr], tpe$TextPost[p], '\t', tpe$Value[p],'\t## ', tpe$Description[p], ' ', df.grp$`Long Name`[gr]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(0, tpe$Value[p]), file=harvest.file.name,append=TRUE, '\n')
            }
          } else if(tpe$Category[p] == 'HabitatGroup'){
            for(gr in 1 : n.grp){
              cat(paste0(tpe$TextPre[p], df.grp$Code[gr], tpe$TextPost[p], '\t', length(hab),'\t## ', tpe$Description[p], ' ', df.grp$`Long Name`[gr]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], length(hab)), file=harvest.file.name,append=TRUE, '\n')
            }
          } else if(tpe$Category[p] == 'ValueFish'){
            for( f in 1 : n.fsh){
              cat(paste0(tpe$TextPre[p], df.fsh$Code[f], tpe$TextPost[p], '\t', tpe$Value[p],'\t## ', tpe$Description[p], ' ', df.fsh$Name[f]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(0, tpe$Value[p]), file=harvest.file.name,append=TRUE, '\n')
            }
          } else if(tpe$Category[p] == 'ValueFishTRUE'){
            for( f in 1 : n.fsh){
              cat(paste0(tpe$TextPre[p], df.fsh$Code[f], tpe$TextPost[p], '\t', tpe$Value[p],'\t## ', tpe$Description[p], ' ', df.fsh$Name[f]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(1, tpe$Value[p]), file=harvest.file.name,append=TRUE, '\n')
            }
          }else if(tpe$Category[p] == 'ValuePort'){
            for( f in 1 : n.ports){
              cat(paste0(tpe$TextPre[p], f, tpe$TextPost[p], '\t', tpe$Value[p],'\t## ', tpe$Description[p]), file=harvest.file.name,append=TRUE, '\n')
              cat(rep(0, tpe$Value[p]), file=harvest.file.name,append=TRUE, '\n')
            }
          } else if(tpe$Category[p] == 'HabitatFish'){
            for( f in 1 : n.fsh){
              cat(paste0(tpe$TextPre[p], df.fsh$Code[f], tpe$TextPost[p], '\t', length(hab),'\t## ', tpe$Description[p], ' ', df.fsh$Name[f]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], length(hab)), file=harvest.file.name,append=TRUE, '\n')
            }
          } else if(tpe$Category[p] == 'GroupSeason'){
            for( gr in 1 : n.grp){
              if(gr == 1) cat('## Adults\n')
              cat(paste0(tpe$TextPre[p], df.grp$Code[gr], tpe$TextPost[p], '\t', 4,'\t## ', tpe$Description[p], ' ', df.grp$`Long Name`[gr]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], 4), file=harvest.file.name,append=TRUE, '\n')
            }
            for( gr in 1 : n.grp){
              if(gr == 1) cat('## Juveniles\n')
              cat(paste0(tpe$TextPre[p], 'j', df.grp$Code[gr], tpe$TextPost[p], '\t', 4,'\t## ', tpe$Description[p], ' ', df.grp$`Long Name`[gr]), file=harvest.file.name,append=TRUE,'\n')
              cat(rep(tpe$Value[p], 4), file=harvest.file.name,append=TRUE, '\n')
            }
          }
        } else if(!tpe$Vector[p]){
          ## Scalar
          if(tpe$Category[p] == 'Fishery'){
            for(f in 1 : n.fsh){
              cat(paste0(tpe$TextPre[p], df.fsh$Code[f], tpe$TextPost[p], '\t', tpe$Value[p], '\t## ', tpe$Description[p], ' ', df.fsh$Name[f]), file=harvest.file.name,append=TRUE, '\n')
            }
          } else if(tpe$Category[p] == 'Group') {
            for(gr in 1 : n.grp){
              cat(paste0(tpe$TextPre[p], df.grp$Code[gr], tpe$TextPost[p], '\t', tpe$Value[p], '\t## ', tpe$Description[p], ' ', df.grp$`Long Name`[gr]), file=harvest.file.name,append=TRUE, '\n')
            }
          }
        }
      }
      if (tpe$Flags[p] == '[GEN]'){
        if(tpe$Vector[p]){
          if(tpe$Category[p] == 'Port'){
            cat(paste0(tpe$TextPre[p], tpe$TextPost[p], '\t', n.ports, '\t## ', tpe$Description[p]), file=harvest.file.name,append=TRUE, '\n')
            cat(rep(tpe$Value[p], n.ports), file=harvest.file.name,append=TRUE, '\n')
          } else if(tpe$Category[p] == 'Group'){
            cat(paste0(tpe$TextPre[p], tpe$TextPost[p], '\t', n.grp, '\t## ', tpe$Description[p]), file=harvest.file.name,append=TRUE, '\n')
            cat(rep(tpe$Value[p], file=harvest.file.name,append=TRUE, n.grp), '\n')
          }
        } else if(!tpe$Vector[p]){
          cat(paste0(tpe$TextPre[p],  '\t', tpe$Value[p], '\t ## ', tpe$Description[p]), file=harvest.file.name,append=TRUE, '\n')
        }
      }
      
      
    } else if(tpe$Template[p]==0){
      
      print("External vector")
      
      if(!is.na(tpe$TextPre[p])) {
        
        if(tpe$TextPre[p] == "flag"){
          
          print("flag")
          
          flag.text <- df.fsh %>% 
            mutate(this_value = paste("flag",Code,"day ",day,sep="")) %>% 
            pull(this_value)
          
          cat("\n",file=harvest.file.name, append=TRUE)
          cat(paste0(flag.text,collapse='\n'), file=harvest.file.name,append=TRUE, '\n')
        }
        
        if(tpe$TextPre[p] == "k_cover"){
          
          print("kcover")
          
          flag.text <- rep(1,n.box) #number of boxes
          
          for(eachfleet in df.fsh$Code) {
            
            cat(paste(tpe$TextPre[p],eachfleet," ",n.box,sep=""), file=harvest.file.name, append=TRUE, '\n')
            cat(paste0(flag.text), file=harvest.file.name,append=TRUE, '\n')
          }
          
        }
        
        if(tpe$TextPre[p] == "flagfish"){
          
          print("flagfish")
          
          flag.text <- df.grp %>% 
            mutate(fish_flag = paste("flagfish", Code," ",IsFished," ## [0] Inactive or [1] Active Fishery", sep="")) %>% 
            pull(fish_flag)
          
          cat("\n",file=harvest.file.name, append=TRUE)
          cat(paste0(flag.text,collapse='\n'), file=harvest.file.name,append=TRUE, '\n')
        }
        
        if(tpe$TextPre[p] == 'target_'){
          
          print("target")
          
          target.data <- read_delim("data/target_vector.prm", delim = "//n", col_names = FALSE)
          
          
          for(eachrow in 1:nrow(target.data)){
            
            this.data <- target.data[eachrow,]
            
            cat(paste0(this.data,collapse='\n'), file=harvest.file.name,append=TRUE, '\n')
          }
          
        }
        
        if(tpe$TextPre[p] == 'CatchTS_agedistrib'){
          
          print("CatchTS_agedistrib")
          
          target.data <- read_delim("data/agedistrib_vector.prm", delim = "//n", col_names = FALSE)
          
          for(eachrow in 1:nrow(target.data)){
            
            this.data <- target.data[eachrow,]
            
            cat(paste0(this.data,collapse='\n'), file=harvest.file.name,append=TRUE, '\n')
          }
        }
        
        
        if(tpe$TextPre[p] == 'mFC_'){
          
          print("mFC")
          
          target.data <- read_delim(this.mfc, delim = "//n", col_names = FALSE)
          
          
          for(eachrow in 1:nrow(target.data)){
            
            this.data <- target.data[eachrow,]
            
            # if(run.type == "historical"){
            #   
            #   if(!grepl("mFC",this.data)){
            #     
            #     this.data <- rep(0,times=59)
            #     
            #   }
            # } 
            #comment out if using mFC, currently set to 0 for historical run
            
            cat(paste0(this.data), file=harvest.file.name,append=TRUE, '\n')
          }
          
          flagf.data <- mfc.tibble %>% 
            mutate(mfc = gsub("mFC_","flagF_",mfc)) %>% 
            mutate_if(is.numeric, funs(ifelse(.>0, 1, 0)))
          
          
          for(eachrow in 1:nrow(flagf.data)){
            
            this.data <- flagf.data[eachrow,] %>% dplyr::select(-mfc) 
            this.name <- flagf.data[eachrow,] %>% dplyr::select(mfc) %>% pull(mfc)
            
            cat(this.name, file=harvest.file.name,append=TRUE, '\n')    
            cat(paste0(this.data), file=harvest.file.name,append=TRUE, '\n')
          }
        }
        if(tpe$TextPre[p] == 'q_'){
          
          print("q_")
          
          target.data <- read_delim("data/q_vector.prm", delim = "//n", col_names = FALSE)
          
          for(eachrow in 1:nrow(target.data)){
            
            this.data <- target.data[eachrow,]
            
            cat(paste0(this.data), file=harvest.file.name,append=TRUE, '\n')
          }
          
        }
        
        if(tpe$TextPre[p] == 'MPA'){
          
          print("MPA")
          
          target.data <- read_delim("data/mpa_vector.prm", delim = "//n", col_names = FALSE)
          
          for(eachrow in 1:nrow(target.data)){
            
            this.data <- target.data[eachrow,]
            
            cat(paste0(this.data,collapse='\n'), file=harvest.file.name,append=TRUE, '\n')
          }
        }
      }
      
      if(!is.na(tpe$TextPost[p])) {
        
        if(tpe$TextPost[p] == "_flagdempelfishery"){
          
          print("_flagdempelfishery")
          
          flag.text <- df.fsh %>% 
            mutate(this_value = paste(Code,"_flagdempelfishery ",flagdempelfishery," #",Name, sep="")) %>% 
            pull(this_value)
          
          cat(paste0(flag.text,collapse='\n'), file=harvest.file.name,append=TRUE, '\n')
        }
        
        
        if(tpe$TextPost[p] == '_mFC_startage'){
          
          print("_mFC_startage")
          
          target.data <- read_delim("data/startage_vector.prm", delim = "//n", col_names = FALSE)
          
          for(eachrow in 1:nrow(target.data)){
            
            this.data <- target.data[eachrow,]
            
            cat(paste0(this.data,collapse='\n'), file=harvest.file.name,append=TRUE, '\n')
          }
          
        }
        
      }
      
    } 
    
  }
}