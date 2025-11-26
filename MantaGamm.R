# --- Setup ------------------------------------------------------------------
# Install only if needed
setwd("/Users/iannag/Desktop/Desktop - Ianna’s MacBook Air/GradSchool/Internship/MantaRaysSDM")

# Install required packages
# install.packages(c(
#   "tidyverse",   # Data manipulation and visualization
#   "mgcv",        # Generalized Additive Models (GAMs)
#   "janitor",     # Data cleaning
#   "lubridate",   # Date handling
#   "DHARMa",      # Residual diagnostics
#   "gratia",      # GAM visualization
#   "performance"  # Model performance metrics
# ))


library(tidyverse)
library(mgcv)        # mgcv::bam
library(janitor)     # clean_names
library(lubridate)   # date parsing
library(DHARMa)      # residual diagnostics
library(gratia)      # plotting smooths
library(performance) # concurvity, other checks

# --- 1. Read & clean data ---------------------------------------------------
# Replace the filename below with your CSV file path
df <- read_csv("manta_all_years_filtered.csv", show_col_types = FALSE)

#add the na handling...###


# Make names safe and consistent (snake_case)
df <- df %>% janitor::clean_names()

# Remove duplicate columns (keeping first occurrence)
df <- df[!duplicated(names(df))]

# Print head to inspect
print(colnames(df))
head(df)

# --- 2. Filter on-effort rows and create presence ---------------------------
# NOTE: you told me "on effort is 0" so we keep on_off_effort == 0.
df <- df %>% filter(on_off_effort == 0)

# Create binary presence from number_of_mantas
df <- df %>%
  mutate(presence = if_else(number_of_mantas > 0, 1L, 0L))

# --- 3. Parse dates and derive day-of-year (cyclic) -------------------------
# Parse dates and ensure consistent format
df <- df %>%
  mutate(
    # First parse all dates to Date type, handling both formats
    date_parsed = case_when(
      str_detect(date, "^\\d{4}-\\d{1,2}-\\d{1,2}$") ~ as.Date(date, format = "%Y-%m-%d"),
      str_detect(date, "^\\d{1,2}/\\d{1,2}/\\d{4}$") ~ as.Date(date, format = "%m/%d/%Y"),
      TRUE ~ as.Date(NA)
    ),
    # Update the original date column to MM/DD/YYYY format
    date = format(date_parsed, "%m/%d/%Y"),
    # Then create day of year and radians
    day_of_year = yday(date_parsed),
    doy_rad = 2 * pi * day_of_year / 365
  )

# Verify the date format is consistent
cat("First few dates after formatting:\n")
head(df$date_2)

# If your date format is not YYYY-MM-DD (ymd), try:
# df$date_parsed <- lubridate::parse_date_time(df$date, orders = c("mdy", "dmy", "ymd HMS", "ymd"))

# --- 4. Clean variable types -------------------------------------------------
# Convert factors / factor-like columns
df <- df %>%
  mutate(
    transect = as.factor(transect),
    tide_stage = as.factor(tide_stage),   # if you prefer numeric for tide_stage, drop as.factor
    lunar_phase = as.factor(lunar_phase)  # if lunar_phase is categorical (new, full, etc.)
  )

# Detection covariates (ensure numeric)
# names curated from your header: visibility_1_5, glare_observer_1, beaufort_wind, sea_state, cloud_cover
# After clean_names(), pick exact names from colnames(df) if they differ.
# Example: visibility_1_5 may become visibility_1_5
glimpse(df %>% select(contains("visibility"), contains("glare"), contains("sea_state"), contains("beaufort")))

# --- 5. Compose GAMM formula -----------------------------------------------
# We'll include:
# - smoothers for continuous habitat variables
# - factor for tide_stage
# - linear detection covariates
# - random effect for transect (s(transect, bs="re"))
# - cyclic smoother for day_of_year (bs="cc")

# Example formula - adapt variable names if your cleaned names differ

# Start with a simple model
simple_formula <- presence ~ 
  s(tavg_c) + 
  s(nearest_inlet_km) + 
  s(tide_height_m) +
  tide_stage

mod_simple <- bam(
  formula = simple_formula,
  data = df,
  family = binomial(link = "logit"),
  method = "fREML"
)

# Check if this works
summary(mod_simple)

# the results of this model show that distance to nearest 
# inlet(p<2e-16) and tide height(p=5.86e-06) are the mose 
# influential variables. tide stage(p=0.549) and temperature(p=0.822)
# did not have much of an impact. since there is low deviance this indicates 
# the model is missing important variables. 


# Check unique values for each variable
var_check <- data.frame(
  variable = c("tavg_c", "tide_height_m", "nearest_inlet_km", "prcp_mm_prev_day",
               "tide_intensity", "lunar_illum", "day_of_year", "visibility_1_5",
               "sea_state", "glare_observer_1_0_3", "beaufort_wind"),
  unique_values = sapply(df[, c("tavg_c", "tide_height_m", "nearest_inlet_km",
                              "prcp_mm_prev_day", "tide_intensity", "lunar_illum",
                              "day_of_year", "visibility_1_5", "sea_state",
                              "glare_observer_1_0_3", "beaufort_wind")], 
                        function(x) length(unique(na.omit(x)))),
  suggested_k = NA
)

# Suggest k (maximum k = unique values - 1, but we'll cap at 5 for safety)
var_check$suggested_k <- pmin(floor(var_check$unique_values * 0.8), 5)
var_check$suggested_k <- ifelse(var_check$suggested_k < 3, 3, var_check$suggested_k)  # Minimum k=3

print(var_check)


sapply(df[, c("tavg_c",
              "tide_height_m",
              "nearest_inlet_km",
              "lunar_illum",
              "chlor_a",
              "tide_stage",
              "prcp_mm_prev_day",
              "tide_intensity",
              "sea_state",
              "transect")],
       function(x) sum(!is.na(x)))

sapply(df[, c("tavg_c",
              "tide_height_m",
              "nearest_inlet_km",
              "lunar_illum",
              "chlor_a",
              "prcp_mm_prev_day")],
       function(x) length(unique(x)))


# Fit the full GAMM model with robust settings
mod_full <- bam(
  presence ~ 
    # Continuous smooths
    s(tavg_c, k = 5) +
    s(tide_height_m, k = 5) +
    s(nearest_inlet_km, k = 5) +
    s(lunar_illum, k = 5, bs = "cc") +
    s(chlor_a, k=5)+
    s(prcp_mm_prev_day, k=3)+

    # Categorical fixed effects
    as.factor(tide_intensity) +
    as.factor(tide_stage)+
    as.factor(sea_state),

  data = df,
  family = binomial(link = "logit"),
  method = "fREML",
  discrete = TRUE,
  select = TRUE
)

print("its working")

# Check model summary
summary(mod_full)

# Model diagnostics
gam.check(mod_full)

# Visualize smooth terms (uncomment to use)
plot(mod_full, pages = 1, scale = 0)

# Check concurvity
concurvity(mod_full)

print(mod_full)

# --- 6. Fit GAMM with bam() -------------------------------------------------
# Use bam for larger datasets; family = binomial for presence/absence
mod_gamm <- bam(
  formula = mod_full,
  data = df,
  family = binomial(link = "logit"),
  method = "fREML"   # fast REML, stable with random effects
  #discrete = TRUE     # speeds up big datasets; comment out if small dataset
)

summary(mod_gamm)

# --- 7. Diagnostics ---------------------------------------------------------
# 7a: check residuals (DHARMa)
sim_res <- simulateResiduals(mod_gamm, plot = FALSE)
plot(sim_res)   # DHARMa diagnostic plots. outputs residuas. residuals from the model output
#run tests on resiuals to make sure output from the gamm is reliable. 
# acf(residuals(m))
# moran.test(residuals(m) ,lw)

# the residuals will allow us to identify obvious outliers.
# we are also going to plot them. [model evaluation process]
#generate a qq plot. we are going to want to look for some shape. nice scattered residual plot. 
#ideally no distinctive shape. else it is over/under fit. 

# 7b: concurvity checks (nonlinear collinearity)
concur <- mgcv::concurvity(mod_gamm, full = TRUE)
print(concur)

# 7c: check k-index / basis dimension (gives warnings if k too low)
gam.check(mod_gamm)

#if we gwt k values too low we may need to refit....
#s(x, k=10)
#s(x, k=20)

# --- 8. Visualize smooths and partial effects -------------------------------
# using gratia (nicer plotting)
draw(mod_gamm, select = NULL)  # will plot all smooths; or select = 1 etc.
# Or for single smoother:
# draw(mod_gamm, select = "s(tavg_c)")

# For effect tables and p-values
anova(mod_gamm)   # approximate significance of terms
plot(mod_gamm, pages = 1, rug = TRUE, shade = TRUE)
#### more tests!!!! Model evaluation steps!!!!

#Look for residual auto correlation###
# use this as a guide: s(x, k=10). s(x, k=20)

#And spatial auto correlation
#these will output values and we will look at them as they relate to our data 

# --- 9. Interpret important outputs ----------------------------------------
# - summary(mod_gamm) shows estimated degrees of freedom (edf) for smooths and p-values
# - s(variable) with edf > 1 means a non-linear shape was fitted
# - inspect p-values and plots to decide which variables are important

# --- 10. Export predictions (optional) -------------------------------------
# If you later prepare rasters of predictors for prediction, you can use predict(..., type="response")
# Example: predict on the original data frame
df$pred_prob <- predict(mod_full, newdata = df, type = "response")

# Save model
# saveRDS(mod_gamm, "gamm_manta_model.rds")

#step 1:
df$pred <- predict(mod_full, newdata = df, type = "response")

#step 2:
library(ggplot2)

ggplot(df, aes(x = longitude_6, y = latitude_5)) +
  geom_point(aes(color = pred), size = 2, alpha = 0.9) +
  scale_color_viridis_c(option = "plasma", name = "Predicted\nProbability") +
  theme_minimal() +
  coord_fixed() +
  labs(
    title = "Predicted Manta Presence Along Transects",
    x = "Longitude",
    y = "Latitude"
  )

#step 3:
library(dplyr)

# resolution of the grid (adjust if needed)
res <- 0.002   # ~200 m depending on latitude

lon_seq <- seq(min(df$longitude_6), max(df$longitude_6), by = res)
lat_seq <- seq(min(df$latitude_5),  max(df$latitude_5),  by = res)

grid <- expand.grid(
  longitude_6 = lon_seq,
  latitude_5 = lat_seq
)


# Fill grid with average/typical conditions
grid$tavg_c            <- mean(df$tavg_c, na.rm = TRUE)
grid$tide_height_m     <- mean(df$tide_height_m, na.rm = TRUE)
grid$nearest_inlet_km  <- mean(df$nearest_inlet_km, na.rm = TRUE)
grid$lunar_illum       <- mean(df$lunar_illum, na.rm = TRUE)
grid$chlor_a           <- mean(df$chlor_a, na.rm = TRUE)
grid$prcp_mm_prev_day  <- mean(df$prcp_mm_prev_day, na.rm = TRUE)
grid$tide_intensity    <- as.factor(names(sort(table(df$tide_intensity), decreasing=TRUE))[1])
grid$sea_state         <- as.factor(names(sort(table(df$sea_state), decreasing=TRUE))[1])
grid$transect          <- NA  # random effect is ignored for SDM mapping


grid$pred <- predict(mod_full, newdata = grid, type = "response")

#plot heat map
ggplot(grid, aes(longitude_6, latitude_5, fill = pred)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Predicted\nProbability") +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "Species Distribution Model (SDM) — Interpolated Heatmap",
    x = "Longitude",
    y = "Latitude"
  )



