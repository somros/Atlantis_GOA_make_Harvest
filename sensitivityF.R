# Alberto Rovellini
# 02/03/2023
# Test model sensitivity to F
library(tidyverse)
library(data.table)
library(here)
library(readxl)

# fg
atlantis_fg <- read.csv('data/GOA_Groups.csv')
# selex
selex <- read.csv('data/age_at_selex.csv')

# read biom and text. Use last time step?

f_vec <- data.frame('Fval' = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.25, 1.5, 2),
                    'dir' = c(831:836, 845:850))

# do you want to plot with perceived F calculated on exploitable age classes?
selectivity <- TRUE

dir.list <- list.dirs('C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/')

# make data frame with 12 rows (F vals) and as many cols as fished groups you have in the model

fished_grp <- atlantis_fg %>% filter(IsImpacted == 1) %>% pull(Code)
all_cols <- c(paste0(fished_grp, '_Biomass'), paste0(fished_grp, '_Catch'))

# get codes of age structured groups
age_grp <- atlantis_fg %>% filter(NumCohorts > 1) %>% pull(Code)

bio_catch_list <- list()
realized_f_list <- list()

# read in selectivity patterns for IsImpacted
# These are knife-edge curves 
# Note: given the issue with high contribution to total biomass of older age classes, this may result into a small difference from using 
# read in biomass at age for age-structured groups
# TODO: turn this into a function

for(i in 1:nrow(f_vec)){
  
  this_fval <- f_vec[i, 1]
  this_run <- f_vec[which(f_vec$Fval==this_fval),2]
  this_dir <- dir.list[grepl(this_run, dir.list)]
  
  # total biomass
  biodat_tmp <- read.table(paste0(this_dir, '/outputGOA0', this_run, '_testBiomIndx.txt'), sep = ' ', header = T)
  biodat <- biodat_tmp %>% 
    select(Time, any_of(fished_grp)) %>% 
    slice_tail(n = 5) %>% 
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE)))
  
  names(biodat) <- paste0(names(biodat), '_Biomass')
  
  # use GOA biomass by age as an option
  biodat_age_tmp <- read.table(paste0(this_dir, '/outputGOA0', this_run, '_testAgeBiomIndx.txt'), sep = ' ', header = T)
  biodat_age <- biodat_age_tmp %>% 
    slice_tail(n = 5) %>% 
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = 'Code.Age', values_to = 'biomass_mt') %>%
    separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    left_join(selex, by = 'Code') %>%
    mutate(idx = as.numeric(Age) - as.numeric(age_class)) %>%
    filter(is.na(idx) | idx >= 0) %>%
    group_by(Code) %>%
    summarise(biomass_mt = sum(biomass_mt)) %>%
    ungroup() %>%
    pivot_wider(names_from = Code, values_from = biomass_mt) %>% 
    select(any_of(fished_grp))
  
  names(biodat_age) <- paste0(names(biodat_age), '_Biomass')
  
  if(selectivity) biodat <- biodat_age

  # total catch
  catchdat_tmp <- read.table(paste0(this_dir, '/outputGOA0', this_run, '_testCatch.txt'), sep = ' ', header = T)
  catchdat <- catchdat_tmp %>% 
    select(Time, any_of(fished_grp)) %>% 
    slice_tail(n = 5) %>% 
    summarise(across(-"Time", ~ mean(.x, na.rm = TRUE)))
    
  names(catchdat) <- paste0(names(catchdat), '_Catch')
  
  this_df <- data.frame(this_fval, biodat, catchdat)
  
  bio_catch_list[[i]] <- this_df
  
  # get realized F at t1
  biom_t1 <- biodat_tmp %>% 
    select(Time, any_of(fished_grp)) %>% 
    filter(Time == 0) %>%
    #filter(between(Time, 365, 1825)) %>% 
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(-Time, names_to = 'Code', values_to = 'biomass') %>%
    select(-Time)
  
  # and for exploited age classes
  biom_age_t1 <- biodat_age_tmp %>% 
    filter(Time == 0) %>% 
    pivot_longer(-Time, names_to = 'Code.Age', values_to = 'biomass') %>%
    separate_wider_delim(Code.Age, delim = '.', names = c('Code', 'Age')) %>%
    left_join(selex, by = 'Code') %>%
    mutate(idx = as.numeric(Age) - as.numeric(age_class)) %>%
    filter(is.na(idx) | idx >= 0) %>%
    group_by(Code) %>%
    summarise(biomass = sum(biomass)) %>%
    ungroup() %>% 
    filter(Code %in% fished_grp)
  
  # catch
  catch_t1 <- catchdat_tmp %>% 
    select(Time, any_of(fished_grp)) %>% 
    filter(Time == 365) %>%
    #filter(between(Time, 365, 1825)) %>% 
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(-Time, names_to = 'Code', values_to = 'catch') %>%
    select(-Time)
  
  # calc realized f
  if(selectivity)biom_t1 <- biom_age_t1
  
  f_t1 <- biom_t1 %>% left_join(catch_t1) %>%
    mutate(exp_rate = catch/biomass,
           f = -log(1-exp_rate),
           this_fval = this_fval) %>%
    select(Code, this_fval, f)
  
  realized_f_list[[i]] <- f_t1
  
}

bio_catch <- rbindlist(bio_catch_list)
realized_f <- rbindlist(realized_f_list)
realized_f[is.nan(realized_f$f)]$f <- NA # change NaN to NA

to_plot <- atlantis_fg %>% filter(IsImpacted == 1 & GroupType == 'FISH') %>% pull(Code)

# plot
for_plot <- bio_catch %>%
  pivot_longer(-this_fval) %>%
  separate_wider_delim(name, delim ='_', names = c('Code', 'Type')) %>%
  left_join(atlantis_fg %>% select(Code, Name), by = 'Code') %>%
  filter(Code %in% to_plot) %>%
  left_join(realized_f, by = c('Code','this_fval'))

# # make set with Code and F where biomass drops below 50%
# f50_frame <- list()
# 
# for(i in 1:length(to_plot)){
#   
#   code <- to_plot[i]
#   B0 <- for_plot %>% filter(Code == code & this_fval == 0 & Type == 'Biomass') %>% pull(value)
#   f50 <- for_plot %>% filter(Code == code & Type == 'Biomass' & value < B0 / 2) %>% slice_head() %>% pull(this_fval)
#   
#   if(length(f50) == 0) f50 <- NA
#   
#   f50_frame[[i]] <- data.frame('Code' = code, f50)
#   
# }
# 
# f50_frame <- data.table::rbindlist(f50_frame)
# f50_frame <- f50_frame %>% left_join(atlantis_fg %>% select(Code, Name))

# make set with Code and F where catch is highest
fmsy_atlantis <- list()

for(i in 1:length(to_plot)){
  
  code <- to_plot[i]
  this_df <- for_plot %>% filter(Code == code)
  
  # imposed and realized f corresponding to highest catch
  fmsy_atlantis[[i]] <-  this_df %>% filter(Type == 'Catch') %>% arrange(desc(value)) %>% slice_head(n = 1) %>% select(Code, this_fval, f)
  
}

fmsy_atlantis <- data.table::rbindlist(fmsy_atlantis)
fmsy_atlantis <- fmsy_atlantis %>% left_join(atlantis_fg %>% select(Code, Name))

# FMSY from FMP
# read in F and M
tier3 <- read_xlsx('data/GOA MSY estimates tables.xlsx', sheet = 1, range = 'A3:J19') %>%
  select(Stock, FOFL) %>%
  set_names(c('Stock', 'FMSY'))

tier4_5 <- read_xlsx('data/GOA MSY estimates tables.xlsx', sheet = 2, range = 'A3:I10') %>%
  select(`Stock/Stock complex`, `M or FMSY`)%>%
  set_names(c('Stock', 'FMSY'))

tier_3_4_5 <- rbind(tier3, tier4_5)

# make key
tier_3_4_5 <- tier_3_4_5 %>%
  mutate(Code = c('POL','COD','SBF','FFS','FFS','FFS','FFS','FFD',
                  'REX','REX','ATF','FHS','POP','RFS','RFS','RFP',
                  'FFS','RFD','RFD','RFD','RFD','THO','DOG')) %>%
  group_by(Code) %>%
  summarise(FMSY = mean(FMSY)) %>%
  set_names('Code','FMSY')

# get M for other groups from parameters
other_M <- read.csv('data/life_history_parameters.csv') %>%
  select(Code, M_FUNC) %>%
  filter(Code %in% setdiff(atlantis_fg$Code, tier_3_4_5$Code)) %>%
  set_names('Code','FMSY')

all_f <- rbind(tier_3_4_5, other_M)

fmsy <- data.frame('Code' = to_plot) %>%
  left_join(all_f) %>% 
  left_join(atlantis_fg %>% select(Code, Name)) %>%
  mutate(FMSY_25 = FMSY/4)

# make plots
  
p <- for_plot %>%
  filter(this_fval <= 1) %>%
  ggplot(aes(x = f, y = value / 1000))+
  geom_line()+
  geom_point(aes(color = as.factor(this_fval)), size = 2)+
  geom_vline(data = fmsy_atlantis, aes(xintercept = f, group = Name), linetype = 'dashed', color = 'black')+
  geom_vline(data = fmsy, aes(xintercept = FMSY, group = Name), linetype = 'dashed', color = 'orange')+
  geom_vline(data = fmsy, aes(xintercept = FMSY_25, group = Name), linetype = 'dashed', color = 'blue')+
  theme_bw()+
  labs(x = 'F as perceived by the model', y = '1000\'s of tons', color = 'F as model input')+
  ggh4x::facet_grid2(Name~Type, scales = 'free', independent = 'all')

p
ggsave('sensitivity_stagger_SELEX.png', p, width = 8, height = 40)

# broken down in smaller groups
# Tier 3
# split in 2 for slides
t3_1 <- c('POL','COD','SBF','POP','RFD','RFP')

p1_1 <- for_plot %>%
  filter(this_fval <= 1) %>%
  filter(Code %in% t3_1) %>%
  ggplot(aes(x = f, y = value / 1000))+
  geom_line()+
  geom_point(aes(color = as.factor(this_fval)), size = 2)+
  geom_vline(data = fmsy_atlantis %>% filter(Code %in% t3_1), aes(xintercept = f, group = Name), linetype = 'dashed', color = 'black')+
  geom_vline(data = fmsy %>% filter(Code %in% t3_1), aes(xintercept = FMSY, group = Name), linetype = 'dashed', color = 'orange')+
  geom_vline(data = fmsy %>% filter(Code %in% t3_1), aes(xintercept = FMSY_25, group = Name), linetype = 'dashed', color = 'blue')+
  theme_bw()+
  labs(x = 'F as perceived by the model', y = '1000\'s of tons', color = 'F as model input')+
  ggh4x::facet_grid2(Name~Type, scales = 'free', independent = 'all')

p1_1
ggsave('sensitivity_tier3_1_1.png', p1_1, width = 5, height = 7)

t3_2 <- c('FFS','FFD','REX','ATF','FHS')

p1_2 <- for_plot %>%
  filter(this_fval <= 1) %>%
  filter(Code %in% t3_2) %>%
  ggplot(aes(x = f, y = value / 1000))+
  geom_line()+
  geom_point(aes(color = as.factor(this_fval)), size = 2)+
  geom_vline(data = fmsy_atlantis %>% filter(Code %in% t3_2), aes(xintercept = f, group = Name), linetype = 'dashed', color = 'black')+
  geom_vline(data = fmsy %>% filter(Code %in% t3_2), aes(xintercept = FMSY, group = Name), linetype = 'dashed', color = 'orange')+
  geom_vline(data = fmsy %>% filter(Code %in% t3_2), aes(xintercept = FMSY_25, group = Name), linetype = 'dashed', color = 'blue')+
  theme_bw()+
  labs(x = 'F as perceived by the model', y = '1000\'s of tons', color = 'F as model input')+
  ggh4x::facet_grid2(Name~Type, scales = 'free', independent = 'all')

p1_2
ggsave('sensitivity_tier3_2_1.png', p1_2, width = 5, height = 7)

# for methods
p1_3 <- for_plot %>%
  filter(this_fval <= 1) %>%
  filter(Code %in% c(t3_1, t3_2)) %>%
  ggplot(aes(x = f, y = value / 1000))+
  geom_line()+
  geom_point(aes(color = as.factor(this_fval)), size = 2)+
  geom_vline(data = fmsy_atlantis %>% filter(Code %in% c(t3_1, t3_2)), aes(xintercept = f, group = Name), linetype = 'dashed', color = 'black')+
  geom_vline(data = fmsy %>% filter(Code %in% c(t3_1, t3_2)), aes(xintercept = FMSY, group = Name), linetype = 'dashed', color = 'orange')+
  geom_vline(data = fmsy %>% filter(Code %in% c(t3_1, t3_2)), aes(xintercept = FMSY_25, group = Name), linetype = 'dashed', color = 'blue')+
  theme_bw()+
  labs(x = 'F as perceived by the model', y = '1000\'s of tons', color = 'F as model input')+
  ggh4x::facet_grid2(Name~Type, scales = 'free', independent = 'all')

p1_3
ggsave('sensitivity_tier3_all.png', p1_3, width = 7.5, height = 10)


# Tier 4 - 5
t4and5 <- c('RFS','THO','DOG','SKB','SKL','SKO')

p2 <- for_plot %>%
  filter(this_fval <= 1) %>%
  filter(Code %in% t4and5) %>%
  ggplot(aes(x = f, y = value / 1000))+
  geom_line()+
  geom_point(aes(color = as.factor(this_fval)), size = 2)+
  geom_vline(data = fmsy_atlantis %>% filter(Code %in% t4and5), aes(xintercept = f, group = Name), linetype = 'dashed', color = 'black')+
  geom_vline(data = fmsy %>% filter(Code %in% t4and5), aes(xintercept = FMSY, group = Name), linetype = 'dashed', color = 'orange')+
  geom_vline(data = fmsy %>% filter(Code %in% t4and5), aes(xintercept = FMSY_25, group = Name), linetype = 'dashed', color = 'blue')+
  theme_bw()+
  labs(x = 'F as perceived by the model', y = '1000\'s of tons', color = 'F as model input')+
  ggh4x::facet_grid2(Name~Type, scales = 'free', independent = 'all')

p2
ggsave('sensitivity_tier4and5_1.png', p2, width = 5, height = 7)

# forage fish
ff <- c('CAP','EUL','SAN','HER','FOS')

p3 <- for_plot %>%
  filter(this_fval <= 1) %>%
  filter(Code %in% ff) %>%
  ggplot(aes(x = f, y = value / 1000))+
  geom_line()+
  geom_point(aes(color = as.factor(this_fval)), size = 2)+
  geom_vline(data = fmsy_atlantis %>% filter(Code %in% ff), aes(xintercept = f, group = Name), linetype = 'dashed', color = 'black')+
  geom_vline(data = fmsy %>% filter(Code %in% ff), aes(xintercept = FMSY, group = Name), linetype = 'dashed', color = 'orange')+
  geom_vline(data = fmsy %>% filter(Code %in% ff), aes(xintercept = FMSY_25, group = Name), linetype = 'dashed', color = 'blue')+
  theme_bw()+
  labs(x = 'F as perceived by the model', y = '1000\'s of tons', color = 'F as model input')+
  ggh4x::facet_grid2(Name~Type, scales = 'free', independent = 'all')

p3
ggsave('sensitivity_forage_1.png', p3, width = 5, height = 7)

# migrating fish species
mig <- c('SCH','SCM','SCO','SPI','SSO','HAK')

p4 <- for_plot %>%
  filter(this_fval <= 1) %>%
  filter(Code %in% mig) %>%
  ggplot(aes(x = f, y = value / 1000))+
  geom_line()+
  geom_point(aes(color = as.factor(this_fval)), size = 2)+
  geom_vline(data = fmsy_atlantis %>% filter(Code %in% mig), aes(xintercept = f, group = Name), linetype = 'dashed', color = 'black')+
  geom_vline(data = fmsy %>% filter(Code %in% mig), aes(xintercept = FMSY, group = Name), linetype = 'dashed', color = 'orange')+
  geom_vline(data = fmsy %>% filter(Code %in% mig), aes(xintercept = FMSY_25, group = Name), linetype = 'dashed', color = 'blue')+
  theme_bw()+
  labs(x = 'F as perceived by the model', y = '1000\'s of tons', color = 'F as model input')+
  ggh4x::facet_grid2(Name~Type, scales = 'free', independent = 'all')

p4
ggsave('sensitivity_migration_1.png', p4, width = 5, height = 7)
