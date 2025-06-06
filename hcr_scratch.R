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
ref_run <- 2097
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
#run <- c(2060, 2061, 2064, 2065, 2066, 2067, 2072, 2073, 2074, 2075)
#run <- c(2072, 2073, 2074, 2075, 2080, 2081, 2082, 2083)
run <- c(2097:2112)

# caps
#cap <- c(400, 200, 600, 400, 200, 600, NA, 400, 200, 600) * 1000
#cap <- rep(c(NA, 400, 200, 600), 2) * 1000
cap <- rep(c(200,400,600,NA), 4) * 1000

# weight scheme
#wgts <- c(rep("equal", 3), rep("binary", 4), rep("ramp", 3))
#wgts <- c(rep("ramp", 8))
wgts <- c(rep("equal", 8), rep("binary", 8))

# climate forcings
# env <- c(rep(1999, 4), rep("2014", 4))
env <- rep(c(rep("ssp126", 4), rep("ssp585", 4)), 2)

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

key_config <- set_key(run, cap, wgts, env)

# function to extract fishery info for stocks that are part of the OY
pull_fishery_info <- function(this_run){
  
  print(this_run)
  
  # File paths
  wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  biom_file <- paste0("outputGOA0", this_run, "_testAgeBiomIndx.txt")
  catch_file <- paste0("outputGOA0", this_run, "_testCatch.txt")
  harvest_prm <- list.files(wd)[grep("GOA_harvest_.*.prm", list.files(wd))]
  ecocap_file <- paste0(wd, "/outputGOA0", this_run, "_test_EcosystemCapResult.txt")
  
  # Read files once
  biom <- read.csv(paste(wd, biom_file, sep = "/"), sep = " ", header = T)
  catch <- read.csv(paste(wd, catch_file, sep = "/"), sep = " ", header = T)
  harvest <- readLines(paste(wd, harvest_prm, sep = "/"))
  ecocap_report <- read.delim(ecocap_file, sep = " ")
  
  # Process biomass data once, outside the loop
  # Convert to long format once and use faster string splitting
  biom_long <- biom %>%
    pivot_longer(-Time, names_to = "Code.Age", values_to = "mt")
  
  # Use faster string splitting method
  # Split Code.Age using base R (much faster than separate())
  code_age_split <- strsplit(biom_long$Code.Age, "\\.", fixed = FALSE)
  biom_long$Code <- sapply(code_age_split, `[`, 1)
  biom_long$Age <- as.numeric(sapply(code_age_split, `[`, 2))
  
  # Pre-filter to only relevant species and times
  biom_filtered <- biom_long %>%
    filter(Code %in% oy_species,
           Time > 0, Time < max(Time)) %>%
    select(-Code.Age)  # Remove original column
  
  # Pre-process harvest parameters outside loop
  harvest_params <- list()
  for(sp in oy_species) {
    startage_line <- harvest[grep(paste0(sp, "_mFC_startage"), harvest) + 1]
    startage <- as.numeric(strsplit(startage_line, split = " ")[[1]])[1]
    
    if(sp %in% hcr_spp) {
      sp_idx <- grep(sp, codes)
      estbo_vec <- as.numeric(strsplit(harvest[grep("estBo\t", harvest) + 1], split = " ")[[1]])
      estbo <- estbo_vec[sp_idx]
      fref_vec <- as.numeric(strsplit(harvest[grep("Fref\t", harvest) + 1], split = " ")[[1]])
      fref <- -log(1 - fref_vec[sp_idx])
    } else {
      estbo <- NA
      mfc_line <- harvest[grep(paste0("mFC_", sp, " "), harvest) + 1]
      mfc <- as.numeric(strsplit(mfc_line, split = " ")[[1]])[1]
      fref <- -(365) * log(1 - mfc)
    }
    
    harvest_params[[sp]] <- list(startage = startage, estbo = estbo, fref = fref)
  }
  
  # Process all species at once using vectorized operations
  # Create biomass summaries for all species
  biom_selex_all <- biom_filtered %>%
    left_join(
      data.frame(Code = names(harvest_params),
                 startage = sapply(harvest_params, `[[`, "startage")),
      by = "Code"
    ) %>%
    filter(Age >= startage) %>%
    group_by(Time, Code) %>%
    summarise(biom_mt_selex = sum(mt), .groups = 'drop')
  
  biom_tot_all <- biom_filtered %>%
    group_by(Time, Code) %>%
    summarise(biom_mt_tot = sum(mt), .groups = 'drop')
  
  # Process catch data for all species
  catch_long <- catch %>%
    select(Time, all_of(oy_species)) %>%
    pivot_longer(-Time, names_to = "Code", values_to = "catch_mt")
  
  # join operation
  result_list <- list()
  for(sp in oy_species) {
    sp_params <- harvest_params[[sp]]
    
    sp_data <- biom_selex_all %>%
      filter(Code == sp) %>%
      left_join(filter(catch_long, Code == sp), by = c("Time", "Code")) %>%
      left_join(filter(biom_tot_all, Code == sp), by = c("Time", "Code")) %>%
      mutate(
        mu = catch_mt / biom_mt_selex,
        f = -log(1 - mu),
        biom_frac = biom_mt_tot / sp_params$estbo,
        fref = sp_params$fref,
        run = this_run
      )
    
    result_list[[sp]] <- sp_data
  }
  
  # OPTIMIZATION 7: Single rbind instead of iterative rbind
  res_df <- do.call(rbind, result_list)
  
  # Process ecocap report
  ecocap_report <- ecocap_report %>%
    select(Time, SpeciesName, FisheryName, SpBased_Frescale, PostSystCap_Frescale) %>%
    filter(SpeciesName %in% oy_species,
           FisheryName %in% oy_fleets) %>%
    mutate(oy_rescale = PostSystCap_Frescale / SpBased_Frescale) %>%
    select(Time, SpeciesName, oy_rescale) %>%
    rename(Code = SpeciesName)
  
  # Final join
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
    ggplot(aes(x = Time/365, y = catch_tot, color = factor(cap), shape = factor(env)))+
    geom_point()+
    scale_shape_manual(values = c(1:length(unique(catch_df$env))))+
    geom_hline(aes(yintercept = cap), linetype = "dashed")+
    theme_bw()+
    scale_y_continuous(limits = c(0,NA))+
    labs(x = "Year", 
         y = "Total catch (mt)",
         color = "Cap",
         shape = "Environment",
         title = "Total catch")+
    facet_grid(~wgts)

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
      filter(!is.na(catch_mt)) %>%
      mutate(f_frac = f / fref) %>%
      select(-fref, -biom_mt_tot, -biom_mt_selex, -mu) %>%
      pivot_longer(-c(Time, Code, Name, run, cap, wgts, env, other)) %>%
      ggplot(aes(x = Time/365, y = value, color = factor(cap), linetype = factor(env))) +
      geom_line() +
      scale_y_continuous(limits = c(0,NA)) +
      facet_grid2(wgts ~ name, scales = "free_y", independent = "y") +
      labs(title = current_name,
           x = "Year", 
           y = "",
           color = "Cap",
           linetype = "Environment") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 7)
      ) 
    
    ggsave(paste0(plotdir, "/by_spp/", current_code, ".png"), p_tmp, 
           width = 12, height = 5,  
           units = "in", dpi = 300)
    
  }
  
  # for HCR species only: HCR plots
  
  for(i in seq_along(hcr_spp)){
    
    current_code <- hcr_spp[i]
    current_name <- grps %>% filter(Code == current_code) %>% pull(Name)
    
    # TODO: sort out faceting - you want the interesting feature to be in the same panel (e.g., facet by weight scheme or by climate scenario?)
    p2 <- catch_df %>%
      filter(Code == current_code, Time >= 10*365, !is.na(catch_mt)) %>%
      ggplot(aes(x = biom_frac, y = f/fref, color = Time/365))+
      geom_point(aes(shape = factor(wgts)), size = 2)+
      scale_shape_manual(values = c(1:length(unique(catch_df$wgts))))+
      scale_color_viridis_c(option = "inferno")+
      geom_vline(xintercept = 0.4, linetype = "dashed", color = "blue", linewidth = 1)+
      geom_vline(xintercept = 0.02, linetype = "dashed", color = "red", linewidth = 1)+
      geom_vline(xintercept = 0.2, linetype = "dotted", color = "orange", linewidth = 1)+
      geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1, color = "steelblue")+
      scale_y_continuous(limits = c(0,NA))+
      theme_bw()+
      labs(x = "B/B0", 
           y = "F/Ftarg", 
           color = "Year",
           shape = "Weight scheme")+
      facet_grid(factor(env)~cap)
    
    ggsave(paste0(plotdir, "/", current_code, "_hcr.png"), p2, 
           width = 12, height = 5, 
           units = "in", dpi = 300)
    
  }
  
  # ecocap rescaling
  p3 <- catch_df %>%
    filter(!is.na(catch_mt)) %>%
    ggplot(aes(x = Time / 365, y = oy_rescale, color = env, linetype = wgts))+
    geom_line()+
    theme_bw()+
    labs(x = "Year", y = "OY-based catch rescaling")+
    facet_grid(Code~cap)
  
  ggsave(paste0(plotdir, "/", "rescale_factor.png"), p3, 
         width = 8, height = 24, 
         units = "in", dpi = 300)
  
}

plot_fishery(catch_df_fast)
























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
