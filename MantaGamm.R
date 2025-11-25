# --- Setup ------------------------------------------------------------------
# Install only if needed


# Install required packages
install.packages(c(
  "tidyverse",   # Data manipulation and visualization
  "mgcv",        # Generalized Additive Models (GAMs)
  "janitor",     # Data cleaning
  "lubridate",   # Date handling
  "DHARMa",      # Residual diagnostics
  "gratia",      # GAM visualization
  "performance"  # Model performance metrics
))


library(tidyverse)
library(mgcv)        # mgcv::bam
library(janitor)     # clean_names
library(lubridate)   # date parsing
library(DHARMa)      # residual diagnostics
library(gratia)      # plotting smooths
library(performance) # concurvity, other checks

# --- 1. Read & clean data ---------------------------------------------------
# Replace the filename below with your CSV file path
df <- read_csv("Trimmed_Aerial_Data_draft.csv", show_col_types = FALSE)

#add the na handling...###


# Make names safe and consistent (snake_case)
df <- df %>% janitor::clean_names()

# Print head to inspect
print(colnames(df))
head(df)

# --- 2. Filter on-effort rows and create presence ---------------------------
# NOTE: you told me "on effort is 0" so we keep on_off_effort == 0.
# If your file uses 1 for on-effort, change the filter to on_off_effort == 1
df <- df %>% filter(on_off_effort == 0)

# Create binary presence from number_of_mantas
df <- df %>%
  mutate(presence = if_else(number_of_mantas > 0, 1L, 0L))

# --- 3. Parse dates and derive day-of-year (cyclic) -------------------------
# Try ymd, if fails adjust (see comments below)
df <- df %>%
  mutate(
    date_parsed = lubridate::ymd(date),          # change if your date format differs
    day_of_year = yday(date_parsed),
    doy_rad = 2 * pi * day_of_year / 365        # optional numeric cyclic transform
  )

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


gamm_formula <- as.formula(
  "presence ~
    s(tavg_c, k=10) +
    s(curr_north_mps, k=8) +
    s(tide_height_m, k=8) +
    s(nearest_inlet_km, k=8) +
    s(prcp_mm_prev_day, k=6) +
    s(tide_intensity, k=6) +
    s(lunar_illum, k=6, bs='cc') +       # cyclic smoother for lunar illumination
    tide_stage +                         # factor
    visibility_1_5 +                     # detection linear covariate
    sea_state +                          # detection linear covariate
    glare_observer_1_0_3+                   # detection linear covariate
    beaufort_wind +
    cloud_cover +
    s(day_of_year, bs='cc', k=12) +      # seasonality
    s(transect, bs='re')"                # random effect for repeated points
)

print(gamm_formula)

# --- 6. Fit GAMM with bam() -------------------------------------------------
# Use bam for larger datasets; family = binomial for presence/absence
mod_gamm <- bam(
  formula = gamm_formula,
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
df$pred_prob <- predict(mod_gamm, newdata = df, type = "response")

# Save model
# saveRDS(mod_gamm, "gamm_manta_model.rds")
