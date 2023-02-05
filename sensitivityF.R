# Alberto Rovellini
# 02/03/2023
# Test model sensitivity to F
library(tidyverse)
library(data.table)
library(here)

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

# plot
p <- bio_catch %>%
  filter(this_fval < 1) %>%
  pivot_longer(-this_fval) %>%
  separate_wider_delim(name, delim ='_', names = c('Code', 'Type')) %>%
  left_join(atlantis_fg %>% select(Code, Name), by = 'Code') %>%
  ggplot(aes(x = this_fval, y = value / 1000, color = Type))+
  geom_line()+
  geom_point()+
  scale_color_manual(values = c('red','blue'))+
  theme_bw()+
  labs(x = 'F', y = '1000\'s of tons')+
  facet_wrap(~Name, scales = 'free', ncol = 3)

p

ggsave('sensitivity.png', p, width = 8, height = 22)
  


