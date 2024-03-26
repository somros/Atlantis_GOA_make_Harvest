library(tidyverse)
library(rbgm)
library(sf)
library(maps)
library(mapdata)
library(viridis)
library(ggmap)
library(ggrepel)

# read in geometry
fl <- 'data/GOA_WGS84_V4_final.bgm'
bgm <- read_bgm(fl)
goa_sf <- box_sf(bgm)
st_crs(goa_sf) <- st_crs(attr(goa_sf$geometry, "crs")$proj)
boundary_boxes <- goa_sf %>% filter(boundary == TRUE) %>% pull(box_id) # get boundary boxes

# #read in ports
# ports <- read.csv("fleets/ports.csv", fileEncoding = "UTF-8-BOM")
# #ports <- ports[ports$port != "Seattle",]
# ports_sf <- ports %>% st_as_sf(coords = c(x = "lon", y = "lat"), crs = "WGS84") %>% st_transform(crs = bgm$extra$projection)

# read in port codes taken from the eLandings website
port_codes <- read.csv("fleets/port_codes_eLandings.csv")
# extract coordinates
register_google(key = 'AIzaSyDUVXbvQOdPZlPDEFocLJs7ildlJyIQxh0')  # Replace with your actual API key
port_codes <- port_codes %>%
  filter(State == "ak") %>% # for now only keep AK ports, but think of Seattle
  mutate(loc = paste(Port.Name, State, sep = ", "))

port_coords <- geocode(port_codes$loc)
ports <- cbind(port_codes, port_coords) %>%
  st_as_sf(coords = c(x = "lon", y = "lat"), crs = "WGS84") %>% 
  st_transform(crs = bgm$extra$projection)

# coast shapefile
coast <- maps::map(database = 'worldHires', regions = c('USA','Canada'), plot = FALSE, fill=TRUE)

coast_sf <- coast %>% 
  st_as_sf(crs = 4326) %>% 
  st_transform(crs = st_crs(goa_sf)) %>% 
  st_combine() %>%
  st_crop(goa_sf %>% st_bbox())
  
# ## little extra step to map labels correctly
# port_coords <- do.call(rbind, st_geometry(ports_sf)) %>% 
#   as_tibble() %>% setNames(c("lon","lat"))
# ports_sf <- cbind(ports_sf, port_coords)

# read fleets and fleet keys
fleets <- readRDS("fleets/fleet_total_catch_atl.RDS")
fleet_key <- read.csv("fleets/fleet_Atlantis.csv")

# for the purposes here, remove "Discards" from the data
fleets <- fleets %>%
  filter(port != "Discard")

# plot pollock trawl and halibut hook and line as examples
toplot <- fleets %>%
  filter(fleet %in% c("NSE_TRW_POL_NROA","ALL_HAL_HAL_ALL"), year == 2018) %>%
  group_by(fleet, box_id) %>%
  summarize(tot_mt = sum(weight_mton)) %>% 
  ungroup()

# join fleet information
toplot <- toplot %>%
  left_join(fleet_key, by = c("fleet"="Fleet"))

# add label for plotting
toplot <- toplot %>% mutate(label = paste(Gear, Primary.Spp, sep = ", "))

# make this dataframe lighter
toplot_lite <- toplot %>%
  dplyr::select(label, box_id, tot_mt)

# fill the data frame with empty catches so that we can map all boxes
filler <- expand.grid(label = unique(toplot_lite$label), box_id = 0:108) %>%
  mutate(tot_mt = 0)

empty_boxes <- anti_join(filler, toplot_lite, by = c("label", "box_id"))
toplot_lite <- rbind(toplot_lite, empty_boxes)

# set catch to NA in boundary boxes
toplot_lite <- toplot_lite %>%
  mutate(tot_mt = ifelse(box_id %in% boundary_boxes, NA, tot_mt))

# make it spatial
toplot_sf <- goa_sf %>%
  dplyr::select(box_id) %>%
  left_join(toplot_lite, by = "box_id")

# joint port information
# this will require making a different data frame
# label and ports to map (say top 4)
topports <- fleets %>%
  filter(fleet %in% c("NSE_TRW_POL_NROA","ALL_HAL_HAL_ALL"), year == 2018) %>%
  group_by(fleet, port) %>%
  summarize(tot_mt = sum(weight_mton)) %>% 
  group_by(fleet) %>%
  slice_max(tot_mt, n = 4) %>%
  ungroup() 

# join to ports and make it spatial
topports_sf <- ports %>%
  left_join(topports, by = c("Port.Code"="port")) %>%
  drop_na() %>%
  left_join(fleet_key, by = c("fleet"="Fleet")) %>%
  mutate(label = paste(Gear, Primary.Spp, sep = ", ")) %>%
  dplyr::select(label, Port.Name)

# plot, color by total catch
p_map <- toplot_sf %>%
  filter(box_id < 92) %>%
  split(.$label) %>%
  purrr::map(~ ggplot()+
               geom_sf(data = ., aes(fill = tot_mt), color = NA)+
               scale_fill_viridis()+
               geom_sf(data = coast_sf)+
               geom_sf(data = topports_sf %>% filter(label == .$label))+
               theme_bw()+
               facet_wrap(~label, nrow = 2)) %>%
  cowplot::plot_grid(plotlist = ., ncol = 1)
p_map

# need this for geom_label_repel
ports_labels <- cbind(topports_sf, st_coordinates(topports_sf))
  
# Define a custom plotting function
plot_function <- function(data, this_label) {
  ggplot(data) +
    geom_sf(aes(fill = tot_mt)) +
    #scale_fill_viridis() +
    geom_sf(data = coast_sf) +
    geom_sf(data = topports_sf %>% filter(label == this_label), size = 2, color = "orange") +
    geom_label_repel(data = ports_labels %>% filter(label == this_label), 
                  aes(x = X, y = Y, label = Port.Name),
                  size = 3,
                  fill = "orange",
                  alpha = 0.8) +
    theme_bw() +
    labs(x = "", y = "")+
    facet_wrap(~label, nrow = 2)
}

# Apply the function to each subset
plots <- toplot_sf %>%
  filter(box_id < 92) %>%
  split(.$label) %>%
  imap(~ plot_function(.x, .y))

# Combine plots
combined_plot <- cowplot::plot_grid(plotlist = plots, ncol = 1)
ggsave("footprint_example.pdf", combined_plot, width = 10, height = 4)
