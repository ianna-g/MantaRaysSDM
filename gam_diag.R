## =========================================================
## GAM DIAGNOSTIC TOOLKIT (mgcv)
## Save as: gam_diag.R
## Usage:
##   source("gam_diagnostics.R")
##   gam_diagnostics(m, data = dat,
##                   response = "y",
##                   time_var = "year",
##                   lon_var = "lon",
##                   lat_var = "lat")
## =========================================================

# ---- Required packages ----
pkgs <- c("mgcv", "gratia", "ggplot2")
new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
lapply(pkgs, require, character.only = TRUE)

# spdep is optional, for Moran's I (spatial autocorrelation of residuals)
if (!requireNamespace("spdep", quietly = TRUE)) {
  message("Package 'spdep' not installed; spatial Moran's I will be skipped.")
}

.gam_residual_df <- function(m, data, response) {
  df <- data

  # Extract core quantities
  df$.fitted_link  <- as.numeric(m$linear.predictors)
  df$.fitted_resp  <- as.numeric(m$fitted.values)
  df$.resid_resp   <- residuals(m, type = "response")
  df$.resid_pear   <- residuals(m, type = "pearson")
  df$.resid_dev    <- residuals(m, type = "deviance")
  df$.resid_std    <- rstandard(m)  # standardized residuals

  if (!is.null(response) && response %in% names(df)) {
    df$.y <- df[[response]]
  } else {
    df$.y <- NA
  }

  df
}

.basic_residual_plots <- function(df, model_label = "GAM") {
  p1 <- ggplot(df, aes(x = .fitted_resp, y = .resid_dev)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.4) +
    labs(x = "Fitted values (response scale)",
         y = "Deviance residuals",
         title = paste(model_label, "- Deviance residuals vs fitted"))

  p2 <- ggplot(df, aes(sample = .resid_dev)) +
    stat_qq(alpha = 0.4) +
    stat_qq_line() +
    labs(title = paste(model_label, "- QQ plot (Deviance residuals)"))

  p3 <- ggplot(df, aes(x = .fitted_link, y = .resid_dev)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.4) +
    labs(x = "Linear predictor (eta)",
         y = "Deviance residuals",
         title = paste(model_label, "- Residuals vs linear predictor"))

  list(p1 = p1, p2 = p2, p3 = p3)
}

.residuals_vs_time <- function(df, time_var, model_label = "GAM") {
  if (is.null(time_var) || !(time_var %in% names(df))) {
    return(NULL)
  }

  p <- ggplot(df, aes_string(x = time_var, y = ".resid_dev")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "loess", se = FALSE, span = 0.8) +
    labs(x = time_var,
         y = "Deviance residuals",
         title = paste(model_label, "- Residuals vs time"))

  p
}

.residuals_map <- function(df, lon_var, lat_var, model_label = "GAM") {
  if (is.null(lon_var) || is.null(lat_var) ||
      !(lon_var %in% names(df)) || !(lat_var %in% names(df))) {
    return(NULL)
  }

  p <- ggplot(df, aes_string(x = lon_var, y = lat_var, color = ".resid_dev")) +
    geom_point() +
    coord_equal() +
    scale_color_gradient2(mid = "grey80") +
    labs(x = lon_var, y = lat_var,
         color = "Dev. resid.",
         title = paste(model_label, "- Spatial pattern of residuals"))

  p
}

.moran_test <- function(df, lon_var, lat_var) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    return(NULL)
  }
  if (is.null(lon_var) || is.null(lat_var) ||
      !(lon_var %in% names(df)) || !(lat_var %in% names(df))) {
    return(NULL)
  }

  coords <- as.matrix(df[, c(lon_var, lat_var)])
  # simple distance-based neighbors (0–100 km, tune as needed)
  nb <- spdep::dnearneigh(coords, d1 = 0, d2 = 100000) # distance in m if coords are projected
  if (any(sapply(nb, length) == 0)) {
    message("Some points have no neighbors within 100 km for Moran's I.")
  }
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

  spdep::moran.test(df$.resid_dev, lw, zero.policy = TRUE)
}

.smooth_and_partial_plots <- function(m) {
  # Smooths
  print("=== Smooth effect plots ===")
  print(gratia::draw(m))

  # Partial residuals
  print("=== Partial residual plots (by smooth) ===")
  pr <- gratia::partial_residuals(m)
  print(gratia::draw(pr))
}

gam_diagnostics <- function(m,
                            data,
                            response,
                            time_var = NULL,
                            lon_var  = NULL,
                            lat_var  = NULL,
                            model_label = "GAM") {
  if (!inherits(m, "gam")) stop("m must be a 'gam' or 'bam' model from mgcv.")
  if (missing(data)) stop("You must provide the data frame used to fit the model.")
  if (missing(response)) stop("You must provide the response column name (string).")

  cat("\n==============================\n")
  cat("GAM DIAGNOSTICS\n")
  cat("==============================\n\n")

  # 1) Model summary
  cat(">>> Model summary:\n")
  print(summary(m))

  # 2) gam.check
  cat("\n>>> gam.check (basis size / k-index / basic residual diagnostics):\n")
  mgcv::gam.check(m)

  # 3) Residual data frame
  df <- .gam_residual_df(m, data, response)

  # 4) Basic residual plots
  cat("\n>>> Basic residual plots (fitted vs residuals, QQ, etc.)\n")
  rb <- .basic_residual_plots(df, model_label = model_label)
  print(rb$p1)
  print(rb$p2)
  print(rb$p3)

  # 5) ACF of residuals (temporal autocorrelation)
  cat("\n>>> ACF of deviance residuals (temporal autocorrelation):\n")
  acf(df$.resid_dev, main = paste(model_label, "- ACF of deviance residuals"))

  # 6) Residuals vs time (if available)
  if (!is.null(time_var) && time_var %in% names(df)) {
    cat("\n>>> Residuals vs time:\n")
    print(.residuals_vs_time(df, time_var, model_label = model_label))
  } else {
    cat("\n(No time_var provided or not in data; skipping residuals vs time.)\n")
  }

  # 7) Residuals vs space (if lon/lat provided)
  if (!is.null(lon_var) && !is.null(lat_var) &&
      lon_var %in% names(df) && lat_var %in% names(df)) {
    cat("\n>>> Spatial residual map:\n")
    print(.residuals_map(df, lon_var, lat_var, model_label = model_label))

    # Optional Moran's I
    if (requireNamespace("spdep", quietly = TRUE)) {
      cat("\n>>> Moran's I test on deviance residuals:\n")
      print(.moran_test(df, lon_var, lat_var))
    } else {
      cat("\n(spdep not installed; skipping Moran's I test.)\n")
    }
  } else {
    cat("\n(No lon/lat provided or missing from data; skipping spatial residual checks.)\n")
  }

  # 8) Smooths and partial residuals
  cat("\n>>> Smooth effect plots and partial residual plots (gratia):\n")
  .smooth_and_partial_plots(m)

  cat("\n=== Diagnostics complete. ===\n")
  invisible(df)
}
