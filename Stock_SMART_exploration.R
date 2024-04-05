# Alberto Rovellini
# 04/05/2024
# This code approaches the calculation of F by fleet in a different way
# It pulls F from https://apps-st.fisheries.noaa.gov/stocksmart?app=homepage
# it then apportions these F's across fleets by the same proportion of catch/tota catch in Adams' catch reconstruction

# One difficulty is understanding the nuances of different stocks (e.g. Rock sole vs Northern Rock sole; composite groups)
# Also this won't help us with stocks that are not here (forage, salmon, etc)

library(tidyverse)
library(readxl)

# Read in data ------------------------------------------------------------

# read in time series from Stock SMART
biom_ss <- read_excel("data/Assessment_TimeSeries_Data_SMART.xlsx", sheet  = 1, range = "A1:CL72", col_names = F)

# reshape this so that we can work with it
biom_ss <- t(biom_ss)
colnames(biom_ss) <- c(biom_ss[1,1:8], 1960:2023)
biom_ss <- biom_ss[-c(1,2),]
biom_ss <- biom_ss %>%
  as.data.frame() %>%
  select(-c("Stock ID", "Assessment ID", "Assessment Month" )) 

rownames(biom_ss) <- 1:nrow(biom_ss)

# pivot
biom_ss_long <- biom_ss %>%
  pivot_longer(-c("Stock Name", "Assessment Year", "Parameter", "Description", "Unit"),
               names_to = "Year", values_to = "Value")
  
# filter to Fmort
f_ss <- biom_ss_long %>%
  filter(Parameter == "Fmort")

# let's explore some properties of this data
f_ss %>% pull(`Stock Name`) %>% unique() 
f_ss %>% pull(Unit) %>% unique() 
# Some assessments are dated / we do not have the latest information / iteration
# There are potential duplicates: Skate, Rock sole, Rex sole. How to avoid double-counting?
# spatial aggregation is problematic: sablefish is all of Alaska
# F is available for only some stocks, not even all the T3 stocks
# Fmort comes in a bunch of different kinds, whose units are: "Fully Selected F", "Exploitation Rate", "Apical F", "Rate", "Fully-Selected F", "Fully-selected F", "Absolute F"
# Spelling aside, there are exploitation rates and proper "F", which we can convert to-from, but I don't know what "Apical F" means

# so, all in all, it's not as easy a data set to use as I thought
  