library(tidyverse)
# add catch for AMSS
catch_data_raw <- read.csv("data/Groundfish Total Catch.csv", fileEncoding = 'UTF-8-BOM')

# add decade
catch_data <- catch_data_raw %>%
  dplyr::select(Year,Species.Group.Name,Catch..mt.) %>%
  mutate(Decade = floor(Year / 10) * 10) 

# filter out mollusks
catch_data <- catch_data %>%
  filter(!Species.Group.Name %in% c("Octopus","Squid"))

# testing region
this_sp <- "Atka Mackerel"
yrs <- unique(catch_data$Year)
catch_data %>%
  filter(Species.Group.Name == this_sp) %>%
  group_by(Year) %>%
  summarise(mt = sum(Catch..mt.)) %>%
  ungroup() %>%
  ggplot(aes(x = Year, y = mt))+
  geom_line()+
  scale_x_continuous(limits = c(min(yrs),max(yrs)))+
  scale_y_continuous(limits = c(0,NA))

# bring in Atlantis grps
key <- read.csv("data/Catch_Atlantis_key.csv", header = T)
catch_data <- catch_data %>%
  left_join(key, by = c("Species.Group.Name" = "Catch_sp"))

# this step gets rid of Halibut and Salmon bycatch
catch_data <- catch_data %>% drop_na() 

# clean the catch data so that we do not have repetitions (e.g. "GOA Big Skate" and "Big Skate")
catch_data <- catch_data %>%
  group_by(Decade, Year, Catch_sp_clean, LongName, is_focal) %>%
  summarise(Catch = sum(Catch..mt.)) %>%
  ungroup()

# check that there are 11 folcal groups (not counting HAL)
catch_data %>% 
  select(LongName,is_focal) %>% 
  distinct() %>% 
  filter(is_focal>0) %>% 
  arrange(LongName) %>%
  nrow() # OK - 12

# Catch proportion over time of 12 focal groups vs all groups
tmp <- catch_data %>%
  group_by(Decade, Year, is_focal) %>%
  summarise(tot_split = sum(Catch, na.rm = T)) %>%
  group_by(Decade, Year) %>%
  mutate(tot = sum(tot_split)) %>%
  ungroup() %>%
  mutate(prop = tot_split / tot)

# produce data frame to turn into table for the appendix
# At the bottom:
# aggregated all OY FMP species
# halibut (need to get it from elsewhere, and parse out GOA - maybe ask Cheryl)
# percentage of Atlantis grps?
# TO DO: double-check that all species are indeed summed up against OY (e.g. octopus etc)

catch_by_decade_spp <- catch_data %>%
  group_by(Decade, Catch_sp_clean, LongName, is_focal) %>%
  summarise(Catch_mean = mean(Catch, na.rm = T),
            Catch_SD = sd(Catch, na.rm = T)) %>%
  ungroup()

catch_by_decade_focal <- catch_data %>%
  group_by(Decade, Year, is_focal) %>%
  summarise(Catch_all_sp = sum(Catch, na.rm = T)) %>%
  group_by(Decade, is_focal) %>%
  summarise(Catch_mean = mean(Catch_all_sp, na.rm = T),
            Catch_SD = sd(Catch_all_sp, na.rm = T)) %>%
  ungroup() %>%
  filter(is_focal == 1)

catch_by_decade_spp_tot <- catch_data %>%
  group_by(Decade, Year) %>%
  summarise(Catch_all_sp = sum(Catch, na.rm = T)) %>%
  group_by(Decade) %>%
  summarise(Catch_mean = mean(Catch_all_sp, na.rm = T),
            Catch_SD = sd(Catch_all_sp, na.rm = T)) %>%
  ungroup()

# rearrange the table for word
final_table <- catch_by_decade_spp %>%
  pivot_wider(
    id_cols = c(Catch_sp_clean, LongName, is_focal),
    names_from = Decade,
    values_from = c(Catch_mean, Catch_SD)
  ) %>%
  mutate(
    Catch1990 = paste0(
      trimws(format(round(Catch_mean_1990, 0), big.mark = ",")),
      " ± ",
      trimws(format(round(Catch_SD_1990, 0), big.mark = ","))
    ),
    Catch2000 = paste0(
      trimws(format(round(Catch_mean_2000, 0), big.mark = ",")),
      " ± ",
      trimws(format(round(Catch_SD_2000, 0), big.mark = ","))
    ),
    Catch2010 = paste0(
      trimws(format(round(Catch_mean_2010, 0), big.mark = ",")),
      " ± ",
      trimws(format(round(Catch_SD_2010, 0), big.mark = ","))
    ),
    Catch2020 = paste0(
      trimws(format(round(Catch_mean_2020, 0), big.mark = ",")),
      " ± ",
      trimws(format(round(Catch_SD_2020, 0), big.mark = ","))
    )
  ) %>%
  select(
    LongName,
    Catch_sp_clean,
    is_focal,
    Catch1990,
    Catch2000,
    Catch2010,
    Catch2020
  )


# arrange
final_table <- final_table %>%
  arrange(-is_focal, LongName)

# If you haven't installed writexl yet, uncomment and run this line:
library(writexl)

write_xlsx(final_table, "catch_summary_table.xlsx")

# get total gf catch over time
catch_data %>%
  group_by(Year) %>%
  summarise(tot = sum(Catch)) %>%
  pull(tot) %>%
  summary()

catch_data %>%
  group_by(Year) %>%
  summarise(tot = sum(Catch)) %>%
  filter(Year > 1999) %>%
  ggplot()+
  geom_line(aes(x = Year, y = tot))+
  scale_y_continuous(limits = c(0,NA))
