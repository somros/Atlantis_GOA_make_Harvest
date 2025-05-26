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

grps <- read.csv("data/GOA_Groups.csv")
codes <- grps %>% pull(Code)

# pick runs
runs <- c(2060, 2061, 2064, 2065, 2066, 2067)

spp <- c("POL","COD","POP") # these are the species to look at for the catch, biomass, and HCR plots
oy_spp <- "POL" # these are the species in the OY cap

# combine all possible factors in a label key

# key 
# this probaby needs manual setup

set_key <- function(runs, cap = NULL, wgts = NULL, gamma = NULL, climate = NULL, other = NULL) {
  # Get the length of run vector
  n <- length(runs)
  
  # Set defaults for all arguments except run
  if (is.null(cap)) cap <- rep(NA, n)
  if (is.null(wgts)) wgts <- rep(NA, n)
  if (is.null(gamma)) gamma <- rep(NA, n)
  if (is.null(climate)) climate <- rep(NA, n)
  if (is.null(other)) other <- rep(NA, n)
  
  key_config <- data.frame(runs, cap, wgts, gamma, climate, other)
  return(key_config)
  
}

key_config <- set_key(runs)




# function to extract HCR-related info
hcr_spp <- c("POL","COD","POP")

hcr_tests <- function(this_run){
  
  wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  biom_file <- paste0("outputGOA0", this_run, "_testAgeBiomIndx.txt")
  catch_file <- paste0("outputGOA0", this_run, "_testCatch.txt")
  harvest_prm <- list.files(wd)[grep("GOA_harvest_.*.prm", list.files(wd))]
  
  biom <- read.csv(paste(wd, biom_file, sep = "/"), sep = " ", header = T)
  catch <- read.csv(paste(wd, catch_file, sep = "/"), sep = " ", header = T)
  harvest <- readLines(paste(wd, harvest_prm, sep = "/"))
  
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
  
  for(sp in hcr_spp){
    
    sp_idx <- grep(sp, codes)
    
    # get start age for selex
    startage <- as.numeric(unlist(strsplit(harvest[grep(paste0(sp, "_mFC_startage"), harvest)+1], split = " ")))[1]
    
    # get estbo
    estbo_vec <- as.numeric(unlist(strsplit(harvest[grep("estBo\t", harvest)+1], split = " ")))
    estbo <- estbo_vec[sp_idx]
    
    # get Fref
    fref_vec <- as.numeric(unlist(strsplit(harvest[grep("Fref\t", harvest)+1], split = " ")))
    fref <- -log(1 - fref_vec[sp_idx]) # turn to real F as there were entered as mu
    
    #biom <- biom %>% select(Time, POL) %>% rename(biom_mt = POL)
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
    
    catch <- catch %>% select(Time, all_of(sp)) %>% rename(catch_mt = sp)
    
    tmp <- left_join(biom_selex, catch, by = "Time") %>% 
      left_join(biom_tot, by = "Time") %>%
      #filter(Time > 365, Time < max(Time)) %>% # drop first 2 and last time step
      mutate(mu = catch_mt / biom_mt_selex,
             run = this_run,
             biom_frac = biom_mt_tot / estbo)
    
  }
  
  
  
  return(tmp)
}







label_key <- data.frame("run" = all_runs,
                        "cap" = c(this_cap))

# pull total catch of the oy species
plot_catch <- function(this_run){
  
  wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  catch_file <- paste0("outputGOA0", this_run, "_testCatch.txt")
  
  catch <- read.csv(paste(wd, catch_file, sep = "/"), sep = " ", header = T)
  
  catch <- catch %>% select(Time, all_of(oy_species))
  
  # pivot
  catch_long <- catch %>%
    pivot_longer(-Time, names_to = "Species", values_to = "mt")
  
  # add up
  catch_sum <- catch_long %>%
    group_by(Time) %>%
    mutate(tot_mt = sum(mt)) %>% 
    ungroup() %>%
    mutate(run = this_run)
  
  return(catch_sum)
}

catch_df <- bind_rows(lapply(all_runs, plot_catch)) %>%
  left_join(label_key, by = "run")


# aff climate key if relevant


# order species
catch_df$Species <- factor(catch_df$Species, levels = oy_species)

# plot all
p1 <- catch_df %>%
  filter(Time >= 50*365, Time <= 100*365) %>%
  #filter(run %in% c(1877,1878)) %>%
  ggplot(aes(x = Time/365, y = tot_mt, color = factor(cap, levels = (unique(catch_df$cap)))))+
  geom_point(size = 2)+#, position = position_dodge(width = 2))+
  scale_color_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7"))+
  #geom_hline(yintercept = this_cap, color = "black", linetype = "dashed", linewidth = 0.5)+
  theme_bw()+
  scale_y_continuous(limits = c(0,400000), breaks = seq(0,4000000,50000))+
  labs(x = "Year", 
       y = "Catch (mt)", 
       title = paste(paste(oy_species, collapse = ", "), "in OY"),
       color = "Cap (mt)")+
  facet_wrap(~wgts, nrow = 1)
p1 
ggsave("../../FHL_Workshop/talk/total_catch.png", p1, width = 9, height = 4)

# plot single species
p2 <- catch_df %>%
  filter(Time >= 50*365, Time <= 100*365) %>%
  #filter(wgts != "10,5,4,1,0.5") %>%
  ggplot()+
  #geom_point(aes(x = Time/365, y = mt, color = factor(cap, levels = (unique(catch_df$cap)))), size = 2)+#, position = position_dodge(width = 2))+
  geom_line(aes(x = Time/365, 
                y = mt, 
                color = factor(cap, levels = (unique(catch_df$cap))), 
                linetype = clim),
            linewidth = 1.2)+
  scale_linetype_manual(values = c(1,2,3))+
  scale_color_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7"))+
  scale_y_continuous(limits = c(0,NA))+
  theme_bw()+
  #geom_hline(yintercept = this_cap, color = "black", linetype = "dashed", linewidth = 1)+
  labs(x = "Year", 
       y = "Catch (mt)", 
       title = paste(paste(oy_species, collapse = ", "), "in OY"),
       color = "Cap (mt)")+
  facet_wrap(~Species, nrow = 2, scales = "free")
p2
ggsave("../../FHL_Workshop/talk/spp_catch_3.png", p2, width = 12, height = 8)

# Terminal catch plot -----------------------------------------------------
# instead of plotting time series, let's see if this looks better
terminal_catch <- catch_df %>%
  group_by(run,cap,wgts,Species)%>%
  slice_max(Time, n=5) %>%
  summarise(mt = mean(mt)) %>%
  ungroup()

# make a reference frame for the catch when the cap is not reached (400k)
ref_frame <- terminal_catch %>%
  filter(cap == 400000) %>%
  rename(mt_nocap = mt) %>%
  select(-run,-cap)

# join
terminal_catch <- terminal_catch %>%
  left_join(ref_frame) %>%
  mutate(rel_catch = mt / mt_nocap)

# plot
p3 <- terminal_catch %>%
  filter(cap<400000) %>%
  mutate(cap = factor(cap, levels = (unique(catch_df$cap)))) %>%
  ggplot(aes(x = Species, y = rel_catch, fill = wgts))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
  #scale_fill_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7"))+
  theme_bw()+
  labs(x = "Stock",
       y = "Ratio to catch without cap",
       title = "Terminal catch",
       fill = paste("Weights"))+
  facet_grid(wgts~cap)
p3
ggsave("../../FHL_Workshop/talk/relcatch.png", p3, width = 9, height = 4.5)

# p4 <- terminal_catch %>%
#   filter(cap < 400000) %>%
#   mutate(cap = factor(cap, levels = (unique(catch_df$cap)))) %>%
#   ggplot(aes(x = Species, y = rel_catch, fill = wgts, alpha = clim)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   scale_alpha_manual(values = c(1,0.5))+
#   geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
#   theme_bw() +
#   labs(x = "Stock", 
#        y = "Ratio to catch without cap", 
#        title = "Terminal catch",
#        fill = paste("Weights"),
#        alpha = "Climate") +
#   facet_grid(wgts ~ cap)
# p4
# ggsave("../../FHL_Workshop/talk/relcatch_climate.png", p4, width = 9, height = 6)

# HCR and reference point plots -------------------------------------------
estbo_key <- data.frame("run"=all_runs,
                        "estbo" = rep(1600000,3))


hcr_tests <- function(this_run, this_sp){
  
  wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
  biom_file <- paste0("outputGOA0", this_run, "_testAgeBiomIndx.txt")
  catch_file <- paste0("outputGOA0", this_run, "_testCatch.txt")
  harvest_prm <- list.files(wd)[grep("GOA_harvest_.*.prm", list.files(wd))]
  
  biom <- read.csv(paste(wd, biom_file, sep = "/"), sep = " ", header = T)
  catch <- read.csv(paste(wd, catch_file, sep = "/"), sep = " ", header = T)
  harvest <- readLines(paste(wd, harvest_prm, sep = "/"))
  startage <- as.numeric(unlist(strsplit(harvest[grep(paste0(this_sp, "_mFC_startage"), harvest)+1], split = " ")))[1]
  
  estbo <- estbo_key %>% filter(run == this_run) %>% pull(estbo)
  
  # isolate POL
  #biom <- biom %>% select(Time, POL) %>% rename(biom_mt = POL)
  # selected biomass
  biom_selex <- biom %>%
    pivot_longer(-Time, names_to = "Code.Age", values_to = "mt") %>%
    separate(`Code.Age`, into = c("Code", "Age"), sep = "\\.") %>%
    filter(Code == this_sp) %>%
    filter(Age >= startage) %>%
    group_by(Time) %>%
    summarise(biom_mt_selex = sum(mt))
  
  # total biomass for B/B0
  biom_tot <- biom %>%
    pivot_longer(-Time, names_to = "Code.Age", values_to = "mt") %>%
    separate(`Code.Age`, into = c("Code", "Age"), sep = "\\.") %>%
    filter(Code == this_sp) %>%
    group_by(Time) %>%
    summarise(biom_mt_tot = sum(mt))
  
  catch <- catch %>% select(Time, all_of(this_sp)) %>% rename(catch_mt = this_sp)
  
  tmp <- left_join(biom_selex, catch, by = "Time") %>% 
    left_join(biom_tot, by = "Time") %>%
    #filter(Time > 365, Time < max(Time)) %>% # drop first 2 and last time step
    mutate(mu = catch_mt / biom_mt_selex,
           run = this_run,
           biom_frac = biom_mt_tot / estbo)
  
  return(tmp)
}

# all_runs <- c(1580:1583)
tt <- bind_rows(lapply(all_runs, hcr_tests, "POL"))

# for plots
tt <- tt %>% 
  left_join(label_key) %>% 
  left_join(key_wgts) %>%
  left_join(key_gamma) %>%
  left_join(key_conf) %>%
  mutate(f = -(1/1)*log(1-mu))

# plot catch and biomass
tt %>%
  select(-mu) %>%
  #filter(Time >= 25*365, Time <= 100*365) %>%
  pivot_longer(-c(Time,run,cap,wgts), values_to = "mt", names_to = "var") %>%
  filter(var %in% c("catch_mt","biom_mt_selex","biom_frac","f")) %>%
  ggplot(aes(x=Time/365, y=mt))+#, color = factor(cap, levels = (unique(catch_df$cap)))))+
  geom_line(aes(), linewidth = 1.5)+#, position = position_dodge(width = 1))+
  scale_color_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7"))+
  theme_bw()+
  labs(x= "Year", y = "", title = "POL", color = "Cap (mt)")+
  facet_wrap(~var, scales = "free")+
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))

# kk <- tt %>%
#   select(-mu) %>%
#   filter(Time >= 30*365, Time <= 80*365) %>%
#   filter(run %in% c(1877,1878))
# plot
# tt %>% 
#   filter(Time > 365*10) %>%
#   #filter(run %in% c(1705, 1724)) %>%
#   ggplot(aes(x = biom_frac, y = mu, color = Time))+
#   geom_vline(xintercept = 0.4, linetype = "dashed", color = "blue", linewidth = 1)+
#   geom_vline(xintercept = 0.02, linetype = "dashed", color = "red", linewidth = 1)+
#   geom_vline(xintercept = 0.2, linetype = "dashed", color = "orange", linewidth = 1)+
#   geom_point()+
#   scale_x_continuous(limits = c(0,1.5))+
#   scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1))+
#   facet_wrap(~run, scales = "free")

hcr_df <- tt
hcr_df$cap <- factor(hcr_df$cap, levels = c("4e+05","3e+05","2e+05","1e+05"))
hcr_df$wgts <- factor(hcr_df$wgts, levels = c("Equal","Binary","Value-based"))
hcr_df$conf <- factor(hcr_df$conf, levels = c("Base","RecDev","RecDev+noCrash"))

p <- hcr_df %>%
  #mutate(cap = factor(cap, levels = c("No cap","100,000"))) %>%
  filter(Time >= 20*365, Time <= 100*365) %>%
  #ggplot()+
  #geom_point(aes(x = biom_frac, y = f))+
  ggplot(aes(x = biom_frac, y = f, color = Time/365))+
  geom_point(aes(shape = factor(gamma)), size = 2)+
  scale_shape_manual(values = c(1,2,3))+
  scale_color_viridis_c(option = "inferno")+
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "blue", linewidth = 1)+
  geom_vline(xintercept = 0.02, linetype = "dashed", color = "red", linewidth = 1)+
  geom_vline(xintercept = 0.2, linetype = "dotted", color = "orange", linewidth = 1)+
  geom_hline(yintercept = 0.3, linetype = "dashed", linewidth = 1, color = "steelblue")+
  # geom_hline(yintercept = 0.51, linetype = "dashed", linewidth = 1, color = "goldenrod")+
  # geom_hline(yintercept = 0.12, linetype = "dashed", linewidth = 1, color = "coral")+
  # scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1))+
  scale_y_continuous(limits = c(0,NA))+
  theme_bw()+
  labs(x = "B/B0", 
       y = "F", 
       title = "POL (realized) HCR",
       color = "Year",
       shape = "Gamma")+
  facet_wrap(~conf)
#facet_wrap(~(factor(cap, levels = (unique(catch_df$cap)))), ncol = 1)
p
ggsave("../../FHL_Workshop/talk/hcr_cod_all.png",p,width = 9,height = 5)

# some calcs
tt %>% filter(run == 1586, Time > 365*10) %>% pull(biom_mt) %>% min()
# F= - (1/t)ln(1-μ)
mu <- 0.34
f = -(1/1)*log(1-mu)
f

1-exp(-0.3)


tt1 <- tt %>%
  filter(run == 1596)

# catch by fleet
# what's happening to catch by fleet compositions with the HCR? DO they scale? How do they relate to mFC?

this_run <- 1659
wd <- paste0("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_", this_run)
catch_fleet_file <- paste0("outputGOA0", this_run, "_testCatchPerFishery.txt")
harvest_prm <- "GOA_harvest_fleets_v2_10YR_HCR.prm"

catch <- read.csv(paste(wd, catch_fleet_file, sep = "/"), sep = " ", header = T)
harvest <- readLines(paste(wd, harvest_prm, sep = "/"))

# isolate POL
catch %>% 
  select(Time, Fishery, POL) %>% rename(catch_mt = POL) %>%
  group_by(Time) %>%
  mutate(tot_mt = sum(catch_mt)) %>%
  ungroup() %>%
  mutate(prop = catch_mt /tot_mt) %>%
  ggplot(aes(x = Time, y = prop, fill = Fishery))+
  geom_bar(stat = "identity", position = "stack")


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
