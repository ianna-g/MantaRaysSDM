# Set working directory
#setwd("c:/Users/iannag/Desktop/Desktop - Ianna’s MacBook Air/GradSchool/Internship/MantaRaysSDM")

library(SpatialKDE)
library(sp)
library(sf)
library(dplyr)
library(tmap)
print("started")

process_dataset <- function(csv_path,
                            out_stem = NULL,
                            cell_size = 100,
                            band_width = 500) {
  dat <- read.csv(
    csv_path,
    na.strings = c("", "NA", "N/A", "NULL", "NaN")
  )

  dat <- dat %>% 
    mutate(
      Longitude = as.numeric(Longitude),
      Latitude  = as.numeric(Latitude),
      `On.Off.Effort` = as.numeric(`On.Off.Effort`),
      manta_pres = as.numeric(`Number.of.Mantas`)
    ) %>% 
    filter(!is.na(Longitude) & !is.na(Latitude)) %>% 
    filter(!is.na(`On.Off.Effort`) & `On.Off.Effort` == 0) %>% 
    filter(!is.na(manta_pres)) %>% 
    st_as_sf(coords = c("Longitude", "Latitude"), dim = "XY") %>% 
    st_set_crs(4326) %>% 
    st_transform(32617) %>% 
    select(manta_pres)

  if (nrow(dat) == 0) {
    message(
      "No rows after filtering (coords present, on-effort == 1, and ",
      "manta_pres > 0). Skipping file: ",
      csv_path
    )
    return(invisible(NULL))
  }

  grid <- dat %>% 
    create_grid_rectangular(
      cell_size = cell_size,
      side_offset = band_width
    )

  kde <- dat %>% 
    kde(
      band_width = band_width,
      kernel = "epanechnikov",
      grid = grid,
      weights = dat$manta_pres
    ) %>% 
    filter(!is.na(kde_value) & kde_value > 0)

  message(
    paste(
      "KDE value range:",
      min(kde$kde_value),
      "to",
      max(kde$kde_value)
    )
  )

  p <- tm_shape(kde) +
    tm_polygons(
      col = "kde_value",
      palette = "YlOrRd",
      style = "cont",
      border.col = NA,
      title = "Manta Density",
      lwd = 0,
      alpha = 1
    )

  dir.create("outputs", showWarnings = FALSE)
  if (is.null(out_stem)) {
    out_stem <- tools::file_path_sans_ext(basename(csv_path))
  }
  png_path <- file.path("outputs", paste0(out_stem, ".png"))

  tmap_mode("plot")
  tmap_save(
    p,
    png_path,
    dpi = 200,
    width = 8,
    height = 8,
    units = "in"
  )

  tmap_mode("view")
  p_view <- tm_basemap("Esri.WorldTopoMap") + p
  html_path <- file.path("outputs", paste0(out_stem, ".html"))
  tmap_save(p_view, html_path)

  invisible(list(
    p_plot = p,
    p_view = p_view,
    kde = kde,
    grid = grid
  ))
}

files_to_run <- c(
  "Trimmed Aerial Data_backup file_08302025.xlsx - 2020.csv",
  "Trimmed Aerial Data_backup file_10152025.xlsx - 2021.csv",
  "Trimmed Aerial Data_backup file_10152025.xlsx - 2022.csv",
  "Trimmed Aerial Data_backup file_10152025.xlsx - 2023.csv",
  "Trimmed Aerial Data_backup file_10152025.xlsx - 2024.csv",
  "Trimmed Aerial Data_backup file_10152025.xlsx - 2025.csv"
)
# we are going to run the same function for each of the files so that we dont have to use 5 different r scripts for the same code
invisible(lapply(files_to_run, function(f) {
  process_dataset(f)
}))

# st_crs(aerial_data)
# st_crs(grid_aerial_data)
# names(kde)
# summary(sf::st_drop_geometry(kde))
# sum(is.na(kde$kde_value))
# range(kde$kde_value, na.rm = TRUE)

# st_intersects(grid_aerial_data, aerial_data, sparse = FALSE)
