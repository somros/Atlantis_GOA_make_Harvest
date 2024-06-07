# some checks for different methods of reading catch from a run
# base
catch_nc_file_base <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1517/outputGOA01517_testCATCH.nc"
bio_nc_file_base <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1517/outputGOA01517_test.nc"
# mpa
catch_nc_file_mpa <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1553/outputGOA01553_testCATCH.nc"
bio_nc_file_mpa <- "C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1553/outputGOA01553_test.nc"
# groups
grps <- read.csv("data/GOA_Groups.csv")
all_fg <- grps %>% filter(GroupType %in% c("FISH", "SHARK")) %>% pull(Name)

# geometry
goa_bgm <- read_bgm("data/GOA_WGS84_V4_final.bgm")
goa_sf <- goa_bgm %>% box_sf()

# catch reconstruction
fleets <- readRDS("fleets/fleet_total_catch_atl.RDS")
fleet_key <- read.csv("data/GOA_fisheries.csv")

# process fishery data

# load function to extract catch from netcdf
source("build_catch_output.R")

# Compare between methods of extracting catch from NECDF -----------------
# the results of these two should be the same - or very similar
check1 <- build_catch_output(catch_nc = catch_nc_file_mpa, 
                             bio_nc = bio_nc_file_mpa, 
                             age_struc = F,
                             relative = F,
                             run = 1553)

check2 <- build_catch_output_v2(catch_nc = catch_nc_file_mpa, 
                                fleet_struc = F,
                                relative = F,
                                run = 1553,
                                key = fleet_key)

check_both <- check1 %>%
  st_set_geometry(NULL) %>%
  select(box_id, Name, mt_tot) %>%
  rename(mt_1 = mt_tot) %>%
  left_join(check2 %>%
              st_set_geometry(NULL) %>%
              select(box_id, Name, mt_tot) %>%
              rename(mt_2 = mt_tot)) 

check_both %>%
  ggplot(aes(x = mt_1, y = mt_2, color = box_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")  # Add 1:1 line


# not identical. Method 1 overestimates the catch - which is unsurprising as it basically does an estimate based on WAA and numbers
# Use method 2 - though method one could be useful to look at catch at age
# TODO: clean this up and keep only method 2 to view catch

# Compare method 2 to CSV -------------------------------------------------
check2_agg <- check2 %>%
  group_by(ts, Name) %>%
  summarise(mt = sum(mt)) %>%
  mutate(check = "nc")

check3 <- read.csv("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1553/outputGOA01553_testCatch.txt", sep = " ", header = T)
check3 <- check3 %>%
  pivot_longer(-Time, names_to = "Code", values_to = "mt") %>%
  filter(Time > 0) %>%
  mutate(Time = Time / 365) %>%
  left_join(grps %>% select(Code,Name)) %>%
  drop_na() %>%
  select(Time, Name, mt) %>%
  rename(ts = Time) %>%
  mutate(check = "flat")

check4 <- rbind(check2_agg, check3)

tt <- check4 %>%
  filter(Name == "Pollock") %>%
  ggplot(aes(x = ts, y = mt, color = check))+
  geom_line()+
  facet_wrap(~Name, scales = "free")

# NC and flat seem identical