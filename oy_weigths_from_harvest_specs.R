# Alberto Rovellini
# 12/19/2024
# Script to explore TAC/ABC ratio in the EBS and GOA
# Based on harvest specification data from AKFIN
library(tidyverse)
library(readxl)
library(viridis)

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
  filter(Area %in% c("Total", "Total (GW)", "GW"), Var %in% c("ABC","TAC","OFL")) %>%
  drop_na() %>%
  mutate(Area = "GOA") %>%
  group_by(Stock,Area,Var,Year)%>%
  summarise(mt = sum(mt))

# plot for AMSS
tac_sp <- unique(goa_specs_tot$Stock)
length(unique(goa_specs_tot$Stock))
colors <- c(viridis(11)[2:10], inferno(11)[2:10], cividis(10)[2:9])
goa_specs_tot$Var <- factor(goa_specs_tot$Var, levels = c("OFL","ABC","TAC"))
p <- goa_specs_tot %>%
  filter(Year >1991) %>%
  filter(Var != "OFL") %>%
  ggplot(aes(x = Year, y = mt, fill = Stock))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = colors)+
  geom_hline(yintercept = 800000, linetype = "dashed", color = "red")+
  theme_bw()+
  labs(x = "", y = "Catch (mt)", fill = "")+
  theme(legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.spacing.y = unit(0, 'cm'))+
  guides(fill = guide_legend(nrow = 7))+
  facet_grid(~Var)
p

# add catch for AMSS
catch_data <- read.csv("data/Groundfish Total Catch.csv", fileEncoding = 'UTF-8-BOM')

# there are a lot of non TAC species reported here, as well as by-catch species, species that are in the FMP but are not groundfish, etc.
# For the purpose of comparing to ABC/TAC plots, we will filter only the species that have a TAC in the harvest specifications
# Map species in the catch to species in the harvest specification data set
tac_key <- read.csv("data/tac_catch_key.csv", header = T)

# process the catch data so that it can be mapped to the harvest specification data
catch_data_short <- catch_data %>%
  select(Year, Species.Group.Name, Catch..mt.) %>%
  left_join(tac_key, by = c("Species.Group.Name" = "Catch_sp")) %>%
  group_by(Year, TAC_sp) %>%
  summarise(mt = sum(Catch..mt., na.rm = T)) %>%
  mutate(Area = "GOA", Var = "Catch") %>%
  select(TAC_sp, Area, Var, Year, mt) %>%
  rename(Stock = TAC_sp)

goa_specs_tot_catch <- goa_specs_tot %>% rbind(catch_data_short) %>% drop_na()

# order factors
goa_specs_tot_catch$Var <- factor(goa_specs_tot_catch$Var, levels = c("OFL", "ABC", "TAC", "Catch"))

p <- goa_specs_tot_catch %>%
  filter(Year >1991) %>%
  filter(Var %in% c("TAC","Catch")) %>%
  ggplot(aes(x = Year, y = mt, fill = Stock))+
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = colors)+
  geom_hline(yintercept = 800000, linetype = "dashed", color = "red")+
  theme_bw()+
  labs(x = "", y = "Catch (mt)", fill = "")+
  theme(legend.position="bottom",
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.spacing.y = unit(0, 'cm'))+
  guides(fill = guide_legend(nrow = 7))+
  facet_grid(~Var)
p

# pivot wider
goa_ratio <- goa_specs_tot %>%
  pivot_wider(id_cols = c(Stock,Area,Year), names_from = Var, values_from = mt) %>%
  mutate(ratio = TAC/ABC)

# view
goa_ratio %>%
  ggplot(aes(x = Year, y = ratio))+
  geom_line()+
  facet_wrap(~Stock)


# Plot specs for GOA ------------------------------------------------------
goa_specs_tot %>%
  ggplot(aes(x = Year, y = mt, color = Var))+
  geom_line(size = 1.5, alpha = 0.7)+
  #scale_color_manual(values = c("#2E86AB", "#F24236"))+
  facet_wrap(~Stock, scales = 'free')

# Compute Catch/TAC for OY weights ----------------------------------------
# First need to map species in this data to the Atlantis groups
# Should this be TAC/ABC? Or Catch/ABC?

grp <- read.csv("data/GOA_Groups.csv")

key <- data.frame("Stock" = unique(goa_specs_tot_catch$Stock),
                  "Code" = c("ATF","DFS","SKB","FFD","RFP","FHS","SKL","RFS","OCT","FFS","RFD","SKO",NA,"POP","COD","RFP","POL","REX","RFS","SBF","SCU","FFS",NA,"RFS","RFS","SQD","THO"))

catch_to_tac <- goa_specs_tot_catch %>%
  filter(!(Stock %in% c("Octopus","Squids","Sharks","Other Species"))) %>%
  left_join(key, by = "Stock") %>%
  group_by(Code,Area,Var,Year) %>%
  summarise(mt = sum(mt)) %>%
  ungroup() %>%
  left_join(grp %>% select(Code,Name), by = "Code") %>%
  filter(Year > 2010, Year <= 2020) %>%
  group_by(Code,Name,Area,Var) %>%
  summarize(mean_mt = mean(mt, na.rm = T)) %>%
  ungroup() %>%
  filter(Var %in% c("TAC","Catch")) %>%
  pivot_wider(names_from = Var, values_from = mean_mt) %>%
  mutate(w = Catch / TAC) %>%
  arrange(-w)

# Case 1: weights are a simple progression from least to most valuable, or a linear transformation thereof
wgts1 <- catch_to_tac %>%
  arrange(w) %>%
  select(Code,Name,w) %>%
  mutate(w = row_number())

# Case 2: organize species into target and by-catch (aka high vs low value, or binary choice)
key_targ <- data.frame("Code" = catch_to_tac$Code,
                       "is_target" = c(1,1,1,1,0,1,1,0,1,1,1,0,0,0,0,0,0,0)) # assuming flatfish is not targeted, which is incorrect

wgts2 <- catch_to_tac %>%
  left_join(key_targ, by = "Code") %>%
  rowwise() %>%
  mutate(w = ifelse(is_target > 0, 10,1)) %>%
  ungroup() %>%
  select(Code,Name,w) %>%
  arrange(-w)
