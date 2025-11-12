library(SpatialKDE)
library(sp)
library(sf)
library(dplyr)
library(tmap)
print("started")
aerial_data <- read.csv("Trimmed Aerial Data_backup file_08302025.xlsx - 2020.csv", na.strings = c("", "NA", "N/A", "NULL", "NaN"))

head(aerial_data)

print(aerial_data)

#disregard any lines that have a blank value for longitude or latitude
aerial_data <- aerial_data %>%
  filter(!is.na(Longitude) & !is.na(Latitude))

aerial_data <- aerial_data %>%
  mutate(
    Longitude = as.numeric(Longitude),
    Latitude  = as.numeric(Latitude),
    mantas = as.numeric("Number of Mantas")
  ) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), dim = "XY") %>%
  st_set_crs(4326) %>%
  st_transform(32617) %>%
  select()

cell_size <- 25
band_width <- 150
print("test")


grid_aerial_data <- aerial_data %>%
  create_grid_rectangular(cell_size = cell_size, side_offset = band_width)

#print(grid_aerial_data)

kde <- aerial_data %>%
  kde(band_width = band_width, kernel = "quartic", grid = grid_aerial_data)

#print(kde)

#plot the kde
#plot(kde)

tmap_mode("plot")
p <- tm_shape(kde) +
  tm_polygons(col = "kde_value", palette = "viridis", title = "KDE Estimate") +
  tm_shape(aerial_data) +
  tm_bubbles(size = 0.1, col = "red")


print(p)  # show in plot window

# Save to file (correct usage)
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