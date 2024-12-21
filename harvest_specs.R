# Alberto Rovellini
# 12/19/2024
# Script to explore TAC/ABC ratio in the EBS and GOA
# Based on harvest specification data from AKFIN
library(tidyverse)
library(readxl)

# the excel files for the harvest specifications are different enough that I'll just handle them separately


# Gulf of Alaska ----------------------------------------------------------
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
  mutate(Area = "GOA")

# pivot wider
goa_ratio <- goa_specs_tot %>%
  pivot_wider(id_cols = c(Stock,Area,Year), names_from = Var, values_from = mt) %>%
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
  filter(Area %in% c("BS","AI")) 

bsai_2 <- bsai_specs %>%
  filter(!Stock %in% c("Pollock","Pacific cod","Sablefish")) %>%
  filter(Area == "BSAI")

# rejoin
bsai_specs <- rbind(bsai_1,bsai_2)

# melt
bsai_specs_long <- bsai_specs %>%
  pivot_longer(-c(Stock,Area), names_to = "Var_Year", values_to = "mt") %>%
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
  pivot_wider(id_cols = c(Stock,Area,Year), names_from = Var, values_from = mt) %>%
  mutate(ratio = TAC/ABC)

# view
bsai_ratio %>%
  ggplot(aes(x = Year, y = ratio, color = Area))+
  geom_line()+
  facet_wrap(~Stock)

# where is POL TAC?
bsai_specs_long %>%
  filter(Stock == "Pollock", Var == "TAC") %>%
  ggplot(aes(x = Year, y = mt, color = Area))+
  geom_line(size = 1.2, alpha = 0.7)

# Historically, POL alone has made up between 50% and 75% of the BS 2MT cap

# Comparison --------------------------------------------------------------
# for key stocks (t3) you'll want to view GOA and BSAI next to each other
# find overlapping stock names between the two FMPs
sp <- c(
  intersect(unique(goa_ratio$Stock), unique( bsai_ratio$Stock)),
  "Rock Sole",
  "Yellowfin Sole",
  "Shallow-water Flatfish",
  "Deep-water Flatfish",
  "Rex Sole"
)

# filter and merge
ak_ratio <- rbind(goa_ratio,bsai_ratio) %>%
  filter(Stock %in% sp)

# view
ak_ratio %>%
  ggplot(aes(x = Year, y = ratio, color = Area))+
  geom_line(size = 1, alpha = 0.7)+
  labs(x = "", y = "TAC/ABC")+
  facet_wrap(~Stock)

# zoom into the last ~10 years
ak_ratio %>%
  filter(Year > 2010) %>%
  ggplot(aes(x = Year, y = ratio, color = Area))+
  geom_line(size = 1, alpha = 0.7)+
  labs(x = "", y = "TAC/ABC")+
  facet_wrap(~Stock)

# Some notes (last 10 years):

# BSAI (where ATTACH works, but does it apply to AI? I don't think so):
# - Pollock in BS: fluctuated between 50% and 100% since 2010
# - Pollock in AI: ~50% fairly stable over the last 10 years. In some years Pollock alone is way over the cap, so ratio has to go down
# - Cod in BS: near 100%
# - Cod in AI: ~70%
# - Sablefish has been high, near 100%, but declining in recent years, likely because of how much of it is there
# - ATF ratio is always low (<25%)
# - POP ratio is fairly high (>80%)
# - Rockfish ratio changes over time, pretty high now
# - Yellowfin sole (big stock) is >60%
# - Rock sole < 40% but getting higher

# GOA: 
# - Pollock: Pollock is always 100%
# - Cod is stable ~75%
# - Sablefish has been high but declining
# - POP, Northern rockfish, shortraker, Rex sole, Deep water flatfish are 100%
# - ATF is increasing and is now >75%
# - Shallow water flatfish (includes rock soles) at 80%
# - FLathead increasing

# Variable across stocks but in general TAC/ABC ratio is higher in the GOA
# There is no ecosystem cap reallocation in the GOA
