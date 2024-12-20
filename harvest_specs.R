# Alberto Rovellini
# 12/19/2024
# Script to explore TAC/ABC ratio in the EBS and GOA
# Based on harvest specification data from AKFIN
library(tidyverse)
library(readxl)

# the excel files for the harvest specifications are different enough that I'll just handle them separately

# read in GOA
goa_specs <- read_xlsx("data/GOA_harvest specs_1986-2024.xlsx", 
                       sheet = 1, 
                       range = "A3:DO131",
                       col_names = F) 

# set names for columns
yrs <- 2024:1986
tags <- paste(c("OFL","ABC","TAC"),rep(yrs,each=3),sep="_")
colnames(goa_specs) <- c("Stock","Area",tags)

# need to clean up all the asterisk and comma garbage first
goa_specs <- goa_specs %>%
  # Convert columns to numeric, handling special cases
  mutate(across(!c(Stock,Area), ~{
    x <- ifelse(tolower(.) %in% c("n/a", "n/a"), NA, .)
    # First remove asterisks
    x <- gsub("\\*", "", x)
    # Then remove commas
    x <- gsub(",", "", x)
    # Convert to numeric
    as.numeric(x)
  }))

# fill stock names
goa_specs <- goa_specs %>%
  fill(Stock, .direction = "down")

# melt
goa_specs_long <- goa_specs %>%
  pivot_longer(-c(Stock,Area), names_to = "Var_Year", values_to = "mt") %>%
  separate(col = Var_Year,
           into = c("Var", "Year"),
           sep = "_") %>%
  mutate(Year = as.numeric(Year))

# keep Total GOA only, and keep ABC and TAC only
goa_specs_tot <- goa_specs_long %>% 
  filter(Area == "Total", Var %in% c("ABC","TAC")) %>%
  drop_na() %>%
  select(-Area)

# pivot wider
goa_ratio <- goa_specs_tot %>%
  pivot_wider(id_cols = c(Stock,Year), names_from = Var, values_from = mt) %>%
  mutate(ratio = TAC/ABC)

# view
goa_ratio %>%
  ggplot(aes(x = Year, y = ratio))+
  geom_line()+
  facet_wrap(~Stock)


# Bering Sea --------------------------------------------------------------
# Some stocks have specs for BSAI, some BS and AI
# read in GOA
bsai_specs <- read_xlsx("data/BSAI_harvest specs_1986-2024.xlsx", 
                       sheet = 1, 
                       range = "A3:GO57",
                       col_names = F) 

# set names for columns
yrs <- 2024:1986
tags <- paste(c("OFL","ABC","TAC", "iTAC", "CDQ"),rep(yrs,each=5),sep="_")
colnames(bsai_specs) <- c("Stock","Area",tags)

# need to clean up all the asterisk and comma garbage first
bsai_specs <- bsai_specs %>%
  # Convert columns to numeric, handling special cases
  mutate(across(!c(Stock,Area), ~{
    x <- ifelse(tolower(.) %in% c("n/a", "n/a"), NA, .)
    # First remove asterisks
    x <- gsub("\\*", "", x)
    # Then remove commas
    x <- gsub(",", "", x)
    # Convert to numeric
    as.numeric(x)
  }))

# fill stock names
bsai_specs <- bsai_specs %>%
  fill(Stock, .direction = "down")

# rearrange a bit
bsai_1 <- bsai_specs %>%
  filter(Stock %in% c("Pollock","Pacific cod","Sablefish")) %>%
  filter(Area %in% c("BS","AI")) %>%
  mutate(Stock_Area = paste(Stock, Area, sep = "_"))

bsai_2 <- bsai_specs %>%
  filter(!Stock %in% c("Pollock","Pacific cod","Sablefish")) %>%
  filter(Area == "BSAI") %>%
  mutate(Stock_Area = paste(Stock, Area, sep = "_"))

# rejoin
bsai_specs <- rbind(bsai_1,bsai_2)

# melt
bsai_specs_long <- bsai_specs %>%
  pivot_longer(-c(Stock,Area,Stock_Area), names_to = "Var_Year", values_to = "mt") %>%
  separate(col = Var_Year,
           into = c("Var", "Year"),
           sep = "_") %>%
  mutate(Year = as.numeric(Year))

# keep ABC and TAC only
bsai_specs_tot <- bsai_specs_long %>% 
  filter(Var %in% c("ABC","TAC")) %>%
  drop_na() 

# pivot wider
bsai_ratio <- bsai_specs_tot %>%
  pivot_wider(id_cols = c(Stock,Area,Stock_Area,Year), names_from = Var, values_from = mt) %>%
  mutate(ratio = TAC/ABC)

# view
bsai_ratio %>%
  ggplot(aes(x = Year, y = ratio, color = Area))+
  geom_line()+
  facet_wrap(~Stock)

# for key stocks (t3) you'll want to view GOA and BSAI next to each other
