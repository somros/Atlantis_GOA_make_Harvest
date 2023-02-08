# Alberto Rovellini
# 02/03/2023
# Test model sensitivity to F
library(tidyverse)
library(data.table)
library(here)
library(readxl)

atlantis_fg <- read.csv('data/GOA_Groups.csv')

# read biom and text. Use last time step?

f_vec <- data.frame('Fval' = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.25, 1.5, 2),
                    'dir' = c(831:836, 845:850))



dir.list <- list.dirs('C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/')

# make data frame with 12 rows (F vals) and as many cols as fished groups you have in the model

fished_grp <- atlantis_fg %>% filter(IsImpacted == 1) %>% pull(Code)
all_cols <- c(paste0(fished_grp, '_Biomass'), paste0(fished_grp, '_Catch'))

bio_catch_list <- list()

for(i in 1:nrow(f_vec)){
  
  this_fval <- f_vec[i, 1]
  this_run <- f_vec[which(f_vec$Fval==this_fval),2]
  this_dir <- dir.list[grepl(this_run, dir.list)]
  
  biodat <- read.table(paste0(this_dir, '/outputGOA0', this_run, '_testBiomIndx.txt'), sep = ' ', header = T)
  biodat <- biodat %>% select(Time, any_of(fished_grp)) %>% slice_tail()
  names(biodat) <- c('Time', paste0(names(biodat)[-1], '_Biomass'))
  
  catchdat <- read.table(paste0(this_dir, '/outputGOA0', this_run, '_testCatch.txt'), sep = ' ', header = T)
  catchdat <- catchdat %>% select(Time, any_of(fished_grp)) %>% slice_tail()
  names(catchdat) <- c('Time', paste0(names(catchdat)[-1], '_Catch'))
  
  this_df <- data.frame(this_fval, biodat[,-1], catchdat[,-1])
  
  bio_catch_list[[i]] <- this_df
  
}

bio_catch <- rbindlist(bio_catch_list)

to_plot <- atlantis_fg %>% filter(IsImpacted == 1 & GroupType == 'FISH') %>% pull(Code)

# plot
for_plot <- bio_catch %>%
  filter(this_fval <= 1) %>%
  pivot_longer(-this_fval) %>%
  separate_wider_delim(name, delim ='_', names = c('Code', 'Type')) %>%
  left_join(atlantis_fg %>% select(Code, Name), by = 'Code') %>%
  filter(Code %in% to_plot) 

# make set with Code and F where biomass drops below 50%
f50_frame <- list()

for(i in 1:length(to_plot)){
  
  code <- to_plot[i]
  B0 <- for_plot %>% filter(Code == code & this_fval == 0 & Type == 'Biomass') %>% pull(value)
  f50 <- for_plot %>% filter(Code == code & Type == 'Biomass' & value < B0 / 2) %>% slice_head() %>% pull(this_fval)
  
  if(length(f50) == 0) f50 <- NA
  
  f50_frame[[i]] <- data.frame('Code' = code, f50)
  
}

f50_frame <- data.table::rbindlist(f50_frame)
f50_frame <- f50_frame %>% left_join(atlantis_fg %>% select(Code, Name))

# FMSY
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

fmsy <- data.frame('Code' = to_plot) %>%
  left_join(tier_3_4_5) %>% 
  left_join(atlantis_fg %>% select(Code, Name))
  
p <- for_plot %>%
  ggplot(aes(x = this_fval, y = value / 1000))+
  geom_line()+
  geom_point()+
  geom_vline(data = f50_frame, aes(xintercept = f50, group = Name), linetype = 'dashed', color = 'red')+
  geom_vline(data = fmsy, aes(xintercept = FMSY, group = Name), linetype = 'dashed', color = 'darkgreen')+
  theme_bw()+
  labs(x = 'F', y = '1000\'s of tons')+
  ggh4x::facet_grid2(Name~Type, scales = 'free', independent = 'all')

p

ggsave('sensitivity.png', p, width = 5, height = 35)
  


