# alberto rovellini
# move to functions later
# compute F reduction factor from OY

library(tidyverse)

all_runs <- c(2060, 2061, 2064, 2065, 2066, 2067)

oy_species <- c("COD","POL","POP","FFS","ATF","FHS","REX","FFD","RFS","RFP","RFD","SBF","THO","SCU","DFS","DFD","SKB","SKO","SKL")
oy_fleets <- "background"

this_cap <- rep(c(400000,200000,600000),2)

label_key <- data.frame("run" = all_runs,
                        "cap" = c(this_cap))

get_oy_scaling <- function(this_run){
  
  rundir <- paste0("../../Parametrization/output_files/data/out_",this_run)
  
  ecocap_report <- read.delim(paste0(rundir, "/outputGOA0", this_run, "_test_EcosystemCapResult.txt"),
                              sep = " ")
  
  # remove col in excess
  ecocap_report <- ecocap_report %>%
    select(Time, SpeciesName, FisheryName, SpBased_Frescale, PostSystCap_Frescale)
  
  # filter and compute rescaling
  ecocap_report <- ecocap_report %>%
    filter(SpeciesName %in% oy_species) %>%
    mutate(Time = Time / 365,
           oy_rescale = PostSystCap_Frescale / SpBased_Frescale)
  
  # add run
  ecocap_report <- ecocap_report %>%
    mutate(run = this_run)
  
  return(ecocap_report)
  
}

scale_df <- bind_rows(lapply(all_runs, get_oy_scaling)) %>%
  left_join(label_key, by = "run")

key_wgts <- data.frame("wgts" = c(rep("EQUAL",3),
                                  rep("BINARY",3)),#,
                       #rep("10,10,10,1,1",4)), # this was supposed to be the opposite but this shows us what happens
                       #rep("1000,100,10,0.1,0.01", 4)),
                       "run" = all_runs)

scale_df <- scale_df %>%
  left_join(key_wgts, by = "run")

# view
scale_df %>%
  filter(Time > 5) %>%
  ggplot(aes(x = Time, y = oy_rescale, color = wgts))+
  geom_line()+
  theme_bw()+
  facet_grid(SpeciesName~cap)

# none of these visualizations are great
# if you do not need time and are interetsed in EOR, can always simplify
# otherwise solution is to print out in a pdf and stretch it vertically

