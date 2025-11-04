# Set working directory
setwd("c:/Users/megal/Desktop/mantas")

library(SpatialKDE)
library(sp)
library(sf)
library(dplyr)
library(tmap)
print("started")
aerial_data <- read.csv("mantas2020.csv", na.strings = c("", "NA", "N/A", "NULL", "NaN"))

#disregard any lines that have a blank value for longitude or latitude
aerial_data <- aerial_data %>%
  filter(!is.na(Longitude) & !is.na(Latitude))

aerial_data <- aerial_data %>%
  mutate(
    Longitude = as.numeric(Longitude),
    Latitude  = as.numeric(Latitude),
    manta_pres = as.numeric(manta_pres)
  ) %>%
  filter(manta_pres > 0) %>%  
  st_as_sf(coords = c("Longitude", "Latitude"), dim = "XY") %>%
  st_set_crs(4326) %>%
  st_transform(32617) %>%
  select(manta_pres)

cell_size <- 50   
band_width <- 500 


grid_aerial_data <- aerial_data %>%
  create_grid_rectangular(cell_size = cell_size, side_offset = band_width)

#print(grid_aerial_data)

kde <- aerial_data %>%
  kde(band_width = band_width, kernel = "epanechnikov", grid = grid_aerial_data, weights = aerial_data$manta_pres)


kde <- kde %>%
  filter(!is.na(kde_value) & kde_value > 0)


print(paste("KDE value range:", min(kde$kde_value), "to", max(kde$kde_value)))



tmap_mode("plot")
p <- tm_shape(kde) +
  tm_polygons(
    col = "kde_value", 
    palette = "YlOrRd",  # Yellow-Orange-Red heatmap colors
    style = "cont",  # Continuous color scale
    border.col = NA,  # Remove cell borders
    title = "Manta Density",
    alpha = 0.8  # Some transparency
  ) +
  tm_shape(aerial_data) +
  tm_bubbles(size = 0.05, col = "black", alpha = 0.5)  # Smaller, subtle points


print(p)  

tmap_save(p, "kde.png", dpi = 200, width = 8, height = 8, units = "in")
tmap_mode("view")
print(p)


st_crs(aerial_data)
st_crs(grid_aerial_data)
print("WHYYY")
names(kde)
summary(sf::st_drop_geometry(kde))
sum(is.na(kde$kde_value))
range(kde$kde_value, na.rm = TRUE)

st_intersects(grid_aerial_data, aerial_data, sparse = FALSE)
