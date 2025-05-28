# Alberto Rovellini
# 05/02/2025
# streamline HCR plotting functions
# Purpose: 
# Visualize catch, F, biomass, and HCR plots if the HCR is on
# Produce plots and metrics by species

# add conditions for possible caps, weight structures, climate types, etc etc
# difficulty here is that there are so many different things we are looking at
# For exmaple, different caps by different HCRs etc etc
# Seems difficult to profuce one sinlge function

library(tidyverse)
library(ggh4x)
library(patchwork)

rm(list = ls())

grps <- read.csv("data/GOA_Groups.csv")
codes <- grps %>% pull(Code)

# species and fleets
# pull these from a reference OY run - comparing them makes sense only if all species are consistent
ref_run <- 2060
oy_dir <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", ref_run)
harvest_prm <- list.files(oy_dir)[grep("GOA_harvest_.*.prm", list.files(oy_dir))]
harvest <- readLines(paste(oy_dir, harvest_prm, sep = "/"))
cap_vec <- as.numeric(unlist(strsplit(harvest[grep("FlagSystCapSP", harvest)+1], split = " ")))
oy_species <- codes[which(cap_vec>0)]

oy_fleets <- "background" # manually set this, it will just be bg for the foreseeable future

# get spp that are managed with the HCRs
hcr_spp <- c()
for(sp in codes){
  # Look for "tierXX\t" pattern to match the tab-delimited format
  pattern <- paste0("tier", sp, "\t")
  matches <- harvest[grep(pattern, harvest, fixed = TRUE)]
  if(length(matches) > 0) {
    tier <- as.numeric(gsub("[^0-9]", "", matches[1]))
    if(tier > 0){hcr_spp <- c(hcr_spp, sp)}
  }
}

# Run properties ----------------------------------------------------------
# Still debating on the best way to do this - typing in here or reading in an excel sheet. Pros and cons to both
# model runs
run <- c(2060, 2061, 2064, 2065, 2066, 2067)

# caps
cap <- c(400, 200, 600, 400, 200, 600) * 1000

# weight scheme
wgts <- c(rep("equal", 3), rep("binary", 3))

# combine all run factors in a label key
set_key <- function(run, cap = NULL, wgts = NULL, env = NULL, other = NULL) { # gamma = NULL, 
  # Get the length of run vector
  n <- length(run)
  
  # Set defaults for all arguments except run
  if (is.null(cap)) cap <- rep(NA, n)
  if (is.null(wgts)) wgts <- rep(NA, n)
  #if (is.null(gamma)) gamma <- rep(NA, n) # this is only useful for HCR shape 5
  if (is.null(env)) env <- rep(NA, n)
  if (is.null(other)) other <- rep(NA, n) # this could be rec devs or or off, etc
  
  key_config <- data.frame(run, cap, wgts, env, other) # gamma, 
  return(key_config)
  
}

key_config <- set_key(run, cap, wgts)



# function to extract fishery info for stocks that are part of the OY
pull_fishery_info <- function(this_run){
  
  wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  biom_file <- paste0("outputGOA0", this_run, "_testAgeBiomIndx.txt")
  catch_file <- paste0("outputGOA0", this_run, "_testCatch.txt")
  harvest_prm <- list.files(wd)[grep("GOA_harvest_.*.prm", list.files(wd))]
  ecocap_file <- paste0(wd, "/outputGOA0", this_run, "_test_EcosystemCapResult.txt")
  
  biom <- read.csv(paste(wd, biom_file, sep = "/"), sep = " ", header = T)
  catch <- read.csv(paste(wd, catch_file, sep = "/"), sep = " ", header = T)
  harvest <- readLines(paste(wd, harvest_prm, sep = "/"))
  ecocap_report <- read.delim(ecocap_file, sep = " ")
  
  res_df <- data.frame()
  for(sp in oy_species){
    
    sp_idx <- grep(sp, codes)
    
    # get start age for selex
    startage <- as.numeric(unlist(strsplit(harvest[grep(paste0(sp, "_mFC_startage"), harvest)+1], split = " ")))[1]
    
    # reference points: pull them from PRM if the species is managed with an HCR
    estbo <- NA
    if(sp %in% hcr_spp){
      # get estbo
      estbo_vec <- as.numeric(unlist(strsplit(harvest[grep("estBo\t", harvest)+1], split = " ")))
      estbo <- estbo_vec[sp_idx]
      
      # get Fref
      fref_vec <- as.numeric(unlist(strsplit(harvest[grep("Fref\t", harvest)+1], split = " ")))
      fref <- -log(1 - fref_vec[sp_idx]) # turn to real F as there were entered as mu
      
    } else { # fref for non-HCR species is from mFC
      
      mfc <- as.numeric(unlist(strsplit(harvest[grep(paste0("mFC_",sp," "), harvest)+1], split = " ")))[1]
      fref <- -(365) * log(1-mfc)
      
    }
    
    # biomass and catch
    # selected biomass
    biom_selex <- biom %>%
      pivot_longer(-Time, names_to = "Code.Age", values_to = "mt") %>%
      separate(`Code.Age`, into = c("Code", "Age"), sep = "\\.") %>%
      filter(Code == sp) %>%
      filter(Age >= startage) %>%
      group_by(Time) %>%
      summarise(biom_mt_selex = sum(mt))
    
    # total biomass for B/B0
    biom_tot <- biom %>%
      pivot_longer(-Time, names_to = "Code.Age", values_to = "mt") %>%
      separate(`Code.Age`, into = c("Code", "Age"), sep = "\\.") %>%
      filter(Code == sp) %>%
      group_by(Time) %>%
      summarise(biom_mt_tot = sum(mt))
    
    # catch
    catch_sp <- catch %>% select(Time, all_of(sp)) %>% rename(catch_mt = sp)
    
    # put it all together
    tmp <- left_join(biom_selex, catch_sp, by = "Time") %>% 
      left_join(biom_tot, by = "Time") %>%
      filter(Time > 0, Time < max(Time)) %>% # drop first and last time step
      mutate(mu = catch_mt / biom_mt_selex, # selected exploitation rate
             f = -log(1-mu), # selected F, as emerging from Atlantis
             biom_frac = biom_mt_tot / estbo, # biomass fraction, computed with total biomass because that is what Atlantis uses by default
             fref = fref,
             Code = sp,
             run = this_run)
    
    res_df <- rbind(res_df, tmp)
    
  }
  
  # pull the oy rescaling factor
  # remove col in excess
  ecocap_report <- ecocap_report %>%
    select(Time, SpeciesName, FisheryName, SpBased_Frescale, PostSystCap_Frescale)
  
  # filter and compute rescaling, then thin out
  ecocap_report <- ecocap_report %>%
    filter(SpeciesName %in% oy_species,
           FisheryName %in% oy_fleets) %>%
    mutate(Time = Time,
           oy_rescale = PostSystCap_Frescale / SpBased_Frescale) %>%
    select(Time, SpeciesName, oy_rescale) %>%
    rename(Code = SpeciesName)
  
  # tie to other info
  res_df <- res_df %>%
    left_join(ecocap_report, by = c("Time", "Code"))
  
  return(res_df)
}

catch_df <- bind_rows(lapply(run, pull_fishery_info))

# bind to key, and to group names
catch_df <- catch_df %>%
  left_join(key_config, by = "run") %>%
  left_join(grps %>% select(Code, Name), by = "Code")

# Plots -------------------------------------------------------------------

plot_fishery <- function(catch_df){
  
  # make a plot directory
  current_base_run <- unique(catch_df$run)[1]
  
  plotdir <- paste0("plots/oy/", current_base_run, "_", Sys.Date())
  if(!dir.exists(plotdir)){
    dir.create(plotdir)
  } else {
    print("This directory exists")
  }
  
  # total catch against cap
  p1 <- catch_df %>%
    group_by(Time,run,cap,wgts,env,other) %>%
    summarise(catch_tot = sum(catch_mt)) %>%
    ggplot(aes(x = Time/365, y = catch_tot, color = factor(cap), shape = wgts))+
    geom_point()+
    scale_shape_manual(values = c(1,2))+
    geom_hline(aes(yintercept = cap), linetype = "dashed")+
    theme_bw()+
    scale_y_continuous(limits = c(0,NA))+
    labs(x = "Year", 
         y = "Total catch (mt)",
         color = "Cap",
         shape = "Weight scheme",
         title = "Total catch")

  ggsave(paste0(plotdir, "/", "oy_tot.png"), p1, 
         width = 8, height = 5,  
         units = "in", dpi = 300)
  
  # by species, biom fraction (need B0), f fraction, catch, biomass (selected), exploitation rate, ...?
  # maybe it would be best to produce this plot one species at a time
  # Loop through each code and create a plot
  
  for(i in seq_along(oy_species)) {
    current_code <- oy_species[i]
    current_name <- grps %>% filter(Code == current_code) %>% pull(Name)
    
    p_tmp <- catch_df %>%
      filter(Code == current_code) %>%
      mutate(f_frac = f / fref) %>%
      select(-fref) %>%
      pivot_longer(-c(Time, Code, Name, run, cap, wgts, env, other)) %>%
      ggplot(aes(x = Time/365, y = value, color = factor(cap), linetype = wgts)) +
      geom_line() +
      facet_wrap(~ name, scales = "free_y") +
      labs(title = current_name,
           x = "Year", 
           y = "",
           color = "Cap",
           linetype = "Weight scheme") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 7)
      ) 
    
    ggsave(paste0(plotdir, "/by_spp/", current_code, ".png"), p_tmp, 
           width = 8, height = 5,  
           units = "in", dpi = 300)
    
  }
  
  # for HCR species only: HCR plots
  p2 <- catch_df %>%
    filter(Code %in% hcr_spp, Time >= 10*365) %>%
    ggplot(aes(x = biom_frac, y = f, color = Time/365))+
    geom_point(aes(shape = factor(wgts)), size = 2)+
    scale_shape_manual(values = c(1,2))+
    scale_color_viridis_c(option = "inferno")+
    geom_vline(xintercept = 0.4, linetype = "dashed", color = "blue", linewidth = 1)+
    geom_vline(xintercept = 0.02, linetype = "dashed", color = "red", linewidth = 1)+
    geom_vline(xintercept = 0.2, linetype = "dotted", color = "orange", linewidth = 1)+
    geom_hline(aes(yintercept = fref), linetype = "dashed", linewidth = 1, color = "steelblue")+
    scale_y_continuous(limits = c(0,NA))+
    theme_bw()+
    labs(x = "B/B0", 
         y = "F", 
         color = "Year",
         shape = "Weight scheme")+
    facet_grid(Name~cap)
  
  ggsave(paste0(plotdir, "/", "hcr.png"), p2, 
         width = 8, height = 5, 
         units = "in", dpi = 300)
  
  # ecocap rescaling
  p3 <- catch_df %>%
    ggplot(aes(x = Time / 365, y = oy_rescale, color = wgts))+
    geom_line()+
    theme_bw()+
    facet_grid(Code~cap)
  
  ggsave(paste0(plotdir, "/", "rescale_factor.png"), p3, 
         width = 8, height = 24, 
         units = "in", dpi = 300)
  
}

plot_fishery(catch_df)
























# Plot total biomass ------------------------------------------------------
# pull total catch of the oy species
plot_biom <- function(this_run){
  
  #this_run <- 1906      
  
  wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  biom_file <- paste0("outputGOA0", this_run, "_testAgeBiomIndx.txt")
  harvest_prm <- list.files(wd)[grep("GOA_harvest_.*.prm", list.files(wd))]
  
  biom <- read.csv(paste(wd, biom_file, sep = "/"), sep = " ", header = T)
  #harvest <- readLines(paste(wd, harvest_prm, sep = "/"))
  #startage <- as.numeric(unlist(strsplit(harvest[grep(paste0(this_sp, "_mFC_startage"), harvest)+1], split = " ")))[1]
  
  #estbo <- estbo_key %>% filter(run == this_run) %>% pull(estbo)
  
  # isolate POL
  #biom <- biom %>% select(Time, POL) %>% rename(biom_mt = POL)
  # selected biomass
  # biom_selex <- biom %>%
  #   pivot_longer(-Time, names_to = "Code.Age", values_to = "mt") %>%
  #   separate(`Code.Age`, into = c("Code", "Age"), sep = "\\.") %>%
  #   filter(Code == this_sp) %>%
  #   filter(Age >= startage) %>%
  #   group_by(Time) %>%
  #   summarise(biom_mt_selex = sum(mt))
  
  # total biomass for B/B0
  biom_tot <- biom %>%
    pivot_longer(-Time, names_to = "Code.Age", values_to = "mt") %>%
    separate(`Code.Age`, into = c("Code", "Age"), sep = "\\.") %>%
    filter(Code %in% oy_species) %>%
    #filter(Code == this_sp) %>%
    group_by(Time) %>%
    mutate(biom_mt_tot = sum(mt)) %>%
    ungroup() %>%
    mutate(run = this_run)#,
  #biom_frac = biom_mt_tot / estbo)
  
  # tmp <- left_join(biom_selex, biom_tot, by = "Time") %>%
  #   mutate(run = this_run,
  #          biom_frac = biom_mt_tot / estbo)
  return(biom_tot)
  
}

biom_df <- bind_rows(lapply(all_runs, plot_biom)) %>%
  left_join(label_key, by = "run")

# order follows the 
biom_df <- biom_df %>% left_join(key_wgts)
biom_df$wgts <- factor(biom_df$wgts, levels = c("Equal","Binary","Value-based"))

# order species
biom_df$Code <- factor(biom_df$Code, levels = oy_species)

# plot all
p1 <- biom_df %>%
  filter(Time >= 50*365, Time <= 100*365) %>%
  #filter(run %in% c(1877,1878)) %>%
  ggplot(aes(x = Time/365, y = biom_mt_tot, color = factor(cap, levels = (unique(catch_df$cap)))))+
  geom_point(size = 2)+#, position = position_dodge(width = 2))+
  scale_color_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7"))+
  geom_hline(yintercept = this_cap, color = "black", linetype = "dashed", linewidth = 0.5)+
  theme_bw()+
  #scale_y_continuous(limits = c(0,400000), breaks = seq(0,4000000,50000))+
  labs(x = "Year", 
       y = "Catch (mt)", 
       title = paste(paste(oy_species, collapse = ", "), "in OY"),
       color = "Cap (mt)")+
  facet_wrap(~wgts, nrow = 1)
p1 
#ggsave("../../FHL_Workshop/talk/total_catch.png", p1, width = 9, height = 4)

# Terminal catch plot -----------------------------------------------------
# instead of plotting time series, let's see if this looks better
terminal_biom <- biom_df %>%
  group_by(Time,run,cap,wgts,Code)%>%
  summarise(mt = sum(mt)) %>%
  group_by(run,cap,wgts,Code) %>%
  slice_max(Time, n=5) %>%
  summarise(mt = mean(mt)) %>%
  ungroup()

# make a reference frame for the catch when the cap is not reached (400k)
ref_frame <- terminal_biom %>%
  filter(cap == 400000) %>%
  rename(mt_nocap = mt) %>%
  select(-run,-cap)

# join
terminal_biom <- terminal_biom %>%
  left_join(ref_frame) %>%
  mutate(rel_biom = mt / mt_nocap)

# plot
p4 <- terminal_biom %>%
  filter(cap<400000) %>%
  mutate(cap = factor(cap, levels = (unique(catch_df$cap)))) %>%
  ggplot(aes(x = Code, y = rel_biom, fill = wgts))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
  #scale_fill_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7"))+
  theme_bw()+
  labs(x = "Stock", 
       y = "Ratio to total biomass without cap", 
       title = "Terminal biomass",
       fill = paste("Weights"))+
  facet_grid(wgts~cap)
p4
ggsave("../../FHL_Workshop/talk/relbiom.png", p4, width = 9, height = 4.5)
