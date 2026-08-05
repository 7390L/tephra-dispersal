
#===============================================================
# ~~~ Tephra deposition analysis from Monteath et al. (2025) ~~~
#===============================================================

#    Monteath, A.J., Jensen, B.J.L., Davies, L.J., Bolton, M.S.M., Hughes,
#    P.D.M., Mackay, H., Edwards, M.E., Finkenbinder, M., Booth, R.K.,
#    Cwynar, L.C., Harvey, J., Pyne‐O’Donnell, S., Papp, C.N., Froese, D.G.,
#    Mallon, G., Amesbury, M.J., & Mayfield, R.J. (2025). Increasing Tephra
#    Deposition in Northeastern North America Points to Atmospheric Circulation
#    Changes at the Early Mid Holocene Transition. Journal of Geophysical
#    Research: Atmospheres, 130(1). https://doi.org/10.1029/2024JD042135

# The script originally presented in the above publication is reproduced here
# with adaptations as needed to suit the Mediterranean dataset

# Contents:
# ------------
#   1. Packages
#   2. Load and clean data
#        Distal filter (>100km)
#        Volcanic source filter (Italian vs all)
#        Min vs max tephra count scenarios
#   3. Age binning
#        Count tephras in 500 yr bins
#        Site coverage/ intensity correction
#   4. Regime shift model
#   5. Slope metrics
#   6. Competing models
#         Total sites model
#         Site type model
#   7. Cache results

#==================
# 1. Load packages
#==================

library(readxl)
library(dplyr)
library(tidyr)
library(segmented)
library(minpack.lm)
library(propagate)


#========================
# 2. Load and clean data
#========================

DATA_FILE  <- "Mediterranean_tephra_data.xlsx"
DEPOSIT_SHEET <- "Tephra data"    # The sheet containing the list of tephra deposits
SITE_SHEET <- "Site spans"    # The sheet containing the list of site spans

MIN_DISTANCE_KM <- 100   # Minimum source-to-deposit distance (km); excludes proximal deposits

# Filter by volcanic source
# "italian" includes only Italian-sourced tephra
# "all" includes all tephras in the dataset regardless of source
ANALYSIS_SCOPE <- "italian"    # Change me: "italian", "all"
ANALYSIS_SCOPE <- match.arg(ANALYSIS_SCOPE, c("italian", "all"))

# Used to tag the cache filename by scope, matching the convention in SPD_final.R
SUFFIX <- paste0("_", ANALYSIS_SCOPE)

ITALIAN_SOURCES_EXACT <- c("Campi Flegrei", "Ischia", "Somma-Vesuvius", "Procida-Vivara",
                           "Campanian volcanic field", "Etna", "Lipari", "Vulcano",
                           "Palinuro Seamount", "Aeolian")    # Exact-match Italian sources
ITALIAN_SOURCE_PATTERN <- paste("Campi Flegrei|Ischia|Somma.Vesuvius|Procida|Campanian",
                                "Phlegraean|Etna|Lipari|Vulcano|Aeolian|Palinuro", sep = "|")    # Loose match, includes uncertain attributions


# Define function to account for deposits with multiple populations
# 'Source volcano' in dataset may contain two or more comma-separated sources
# Allow a multi-source deposit to contribute more than one observation
# However note that these observations are not independent
expand_sources <- function(data) {
  parts <- strsplit(as.character(data$`Source volcano`), ",", fixed = TRUE)    # Split cells after comma if present
  n1 <- trimws(vapply(parts, function(p) p[1], character(1)))    # n1 takes the first name from the cell
  n2 <- trimws(vapply(parts,    # n2 takes the second name
                      function(p) if (length(p) >= 2) p[2] else NA_character_,
                      character(1)))

  prim <- data    # Create a copy of the full dataset that only contains primary populations
  prim$`Source volcano` <- n1
  prim$source_rank <- "primary"

  dist_sec <- suppressWarnings(as.numeric(
    data$`Great circle distance between source and deposit (km) - secondary`))
  sec_ok <- !is.na(n2) & !is.na(dist_sec)      # Identify rows where a secondary population is present
                                               # and has a valid distance
  if (any(sec_ok)) {
    sec <- data[sec_ok, ]    # Create a subset of data only containing secondary population rows
    sec$`Source volcano` <- n2[sec_ok]
    sec$`Great circle distance between source and deposit (km)` <- dist_sec[sec_ok]
    sec$source_rank <- "secondary"
    out <- rbind(prim, sec[, names(prim)])    # Combine primary and secondary populations into one data frame
  } else {
    out <- prim
  }
  out
}

deposits <- expand_sources(read_excel(DATA_FILE, sheet = DEPOSIT_SHEET)) %>%
  mutate(
    `Median age (cal yr BP)`  = suppressWarnings(as.numeric(`Median age (cal yr BP)`)),
    `Minimum age (cal yr BP)` = suppressWarnings(as.numeric(`Minimum age (cal yr BP)`)),
    `Maximum age (cal yr BP)` = suppressWarnings(as.numeric(`Maximum age (cal yr BP)`)),
    distance_km = suppressWarnings(as.numeric(`Great circle distance between source and deposit (km)`))
  ) %>%
  filter(!is.na(`Median age (cal yr BP)`), !is.na(`Minimum age (cal yr BP)`),
         !is.na(`Maximum age (cal yr BP)`), !is.na(distance_km)) %>%   # Essential fields present
  filter(distance_km >= MIN_DISTANCE_KM)                              # Exclude proximal deposits

if (ANALYSIS_SCOPE == "italian") {
  # Minimum scenario = Certain (exact-match) Italian attributions, primary populations only
  tephra_data_min <- deposits %>%
    filter(`Source volcano` %in% ITALIAN_SOURCES_EXACT & source_rank == "primary")

  # Maximum scenario = Also includes uncertain "(?)" attributions and secondary populations
  tephra_data_max <- deposits %>%
    filter(grepl(ITALIAN_SOURCE_PATTERN, `Source volcano`, ignore.case = TRUE))
} else {
  # Minimum scenario = All deposits, primary populations only
  tephra_data_min <- deposits %>%
    filter(source_rank == "primary")

  # Maximum scenario = All deposits, including secondary populations
  tephra_data_max <- deposits
}

site_data <- read_excel(DATA_FILE, sheet = SITE_SHEET) %>%
  mutate(
    `Minimum site coverage (cal yr BP)` = suppressWarnings(
      as.numeric(`Minimum site coverage (cal yr BP)`)
    ),
    `Maximum site coverage (cal yr BP)` = suppressWarnings(
      as.numeric(`Maximum site coverage (cal yr BP)`)
    )
  )


#================
# 3. Age binning
#================

# Define age bins, spanning 0-25 ka in 500 yr increments
bins <- data.frame(bin_center = seq(0, 25000, by = 500)) %>%
  mutate(
    bin_start = bin_center - 250,
    bin_end   = bin_center + 250 - 0.1
  )

# Function to count sites in each bin by type, adapted for marine vs lake vs terrestrial
count_sites_in_bin <- function(start, end, site_data) {
  marine_count <-
    sum(
      site_data$`Minimum site coverage (cal yr BP)` <= end &
        site_data$`Maximum site coverage (cal yr BP)` >= start &
        grepl("Marine", site_data$Type, ignore.case = TRUE),
      na.rm = TRUE
    )
  lake_count <-
    sum(
      site_data$`Minimum site coverage (cal yr BP)` <= end &
        site_data$`Maximum site coverage (cal yr BP)` >= start &
        grepl("Lake", site_data$Type, ignore.case = TRUE),
      na.rm = TRUE
    )
  terrestrial_count <-
    sum(
      site_data$`Minimum site coverage (cal yr BP)` <= end &
        site_data$`Maximum site coverage (cal yr BP)` >= start &
        grepl("Terrestrial", site_data$Type, ignore.case = TRUE),
      na.rm = TRUE
    )
  return(c(marine_count, lake_count, terrestrial_count))
}

# Function to count tephras in each bin. Counts rows (one row = one dated
# deposit, after expand_sources() has already split any multi-source deposit
# into its own row)
count_tephras_in_bin <- function(start, end, tephra_data_min, tephra_data_max) {
  min_bin <- tephra_data_min %>%
    filter(`Median age (cal yr BP)` >= start & `Median age (cal yr BP)` <= end)
  max_bin <- tephra_data_max %>%
    filter(`Median age (cal yr BP)` >= start & `Median age (cal yr BP)` <= end)
  min_tephras <- nrow(min_bin)
  max_tephras <- nrow(max_bin)
  return(c(min_tephras, max_tephras))
}

# Calculate counts for each bin
bins <- bins %>%
  rowwise() %>%
  mutate(
    site_counts = list(count_sites_in_bin(bin_start, bin_end, site_data)),
    # site_counts_1 = number of marine sites; site_counts_2 = number of lake sites; site_counts_3 = number of terrestrial sites
    tephras_count = list(count_tephras_in_bin(bin_start, bin_end, tephra_data_min, tephra_data_max)) # tepras_count_1 = minimum number; tephras_count_2 = maximum number
  ) %>%
  unnest_wider(site_counts, names_sep = "_") %>%
  unnest_wider(tephras_count, names_sep = "_")

# Calculate intensity corrected tephra counts, updated to guard against total_sites == 0
bins <- bins %>%
  mutate(
    total_sites = site_counts_1 + site_counts_2 + site_counts_3,
    intensity_corrected_min_tephras = ifelse(total_sites > 0,
                                             tephras_count_1 / total_sites, 0),
    intensity_corrected_max_tephras = ifelse(total_sites > 0,
                                             tephras_count_2 / total_sites, 0)
  )

# Calculate cumulative intensity corrected tephra counts
bins <- bins %>%
  arrange(desc(bin_center)) %>%
  mutate(
    cumulative_intensity_corrected_min_tephras = cumsum(intensity_corrected_min_tephras),
    cumulative_intensity_corrected_max_tephras = cumsum(intensity_corrected_max_tephras)
  )


#=======================
# 4. Regime shift model
#=======================

# ~ AKA "segmented regression", "piecewise regression", "breakpoint analysis"

# Perform segmented regression for cumulative min tephras
# First fit a simple linear model, then let selgmented() search for the
# best-supported number and location of breakpoints (up to Kmax = 20 candidates),
# selecting the optimal number of breaks by BIC
mylm_min <-
  lm(cumulative_intensity_corrected_min_tephras ~ bin_center,
     data = bins
  )
myselg_min <-
  selgmented(
    mylm_min,
    Kmax = 20,
    type = "bic",
    stop.if = 20,
    plot.ic = TRUE
  )

# Extract psi values for min tephras, i.e. the breakpoints and their uncertainty
psi_min <- myselg_min$psi %>%
  as.data.frame() %>%
  mutate(type = "Min Tephras")

# Perform segmented regression for cumulative max tephras
# Same approach as above, applied to the max-tephras scenario
mylm_max <-
  lm(cumulative_intensity_corrected_max_tephras ~ bin_center,
     data = bins
  )
myselg_max <-
  selgmented(
    mylm_max,
    Kmax = 20,
    type = "bic",
    stop.if = 20,
    plot.ic = TRUE
  )

# Extract psi values for max tephras
psi_max <- myselg_max$psi %>%
  as.data.frame() %>%
  mutate(type = "Max Tephras")

# Combine psi values
all_psi <- bind_rows(psi_min, psi_max)

# Define the x values (ages) at which to generate predictions for a smooth plotted curve:
# the original bin centers, the breakpoint locations, and a dense 50-yr grid across 0-25 ka
high_res_points <- seq(from = 25000, to = 0, by = -50)
all_points_min <- c(bins$bin_center, psi_min$Est., high_res_points)
all_points_max <- c(bins$bin_center, psi_max$Est., high_res_points)

# Get the predicted cumulative tephra count (fit) and 95% confidence interval
# at each x value defined above
preds_min <-
  predict(
    myselg_min,
    newdata = data.frame(bin_center = all_points_min),
    interval = "confidence",
    level = 0.95
  )
pred_df_min <-
  data.frame(
    bin_center = all_points_min,
    pred = preds_min,
    type = "Min Tephras"
  )

preds_max <-
  predict(
    myselg_max,
    newdata = data.frame(bin_center = all_points_max),
    interval = "confidence",
    level = 0.95
  )
pred_df_max <-
  data.frame(
    bin_center = all_points_max,
    pred = preds_max,
    type = "Max Tephras"
  )

# Combine prediction data frames
pred_df_all <- bind_rows(pred_df_min, pred_df_max)

# Attach the fitted y-value to each breakpoint, by matching its x position
# against the predictions data frame, within the same min/max type
all_psi <- all_psi %>%
  left_join(
    dplyr::select(pred_df_all, bin_center, pred.fit, type),
    by = c("Est." = "bin_center", "type" = "type")
  )


#==================
# 5. Slope metrics
#==================

# Pull the slope of each fitted linear segment from the segmented regression models
# Slopes can be used to compare rate of change in deposition between segments

# Extract slopes for cumulative min tephras
slopes_min <- slope(myselg_min)$bin_center

# Extract slopes for cumulative max tephras
slopes_max <- slope(myselg_max)$bin_center

# Ensure both data frames have the same number of rows
n_min <- nrow(slopes_min)
n_max <- nrow(slopes_max)
max_rows <- max(n_min, n_max)

# Pad the shorter data frame with NA values
if (n_min < max_rows) {
  slopes_min <-
    rbind(slopes_min, matrix(NA, nrow = max_rows - n_min, ncol = ncol(slopes_min)))
}
if (n_max < max_rows) {
  slopes_max <-
    rbind(slopes_max, matrix(NA, nrow = max_rows - n_max, ncol = ncol(slopes_max)))
}

# Combine slopes into a single data frame and reverse the sign of the estimates
slopes_df <- data.frame(
  Segment = c(paste0("slope", 1:max_rows), paste0("slope", 1:max_rows)),
  Est. = -c(slopes_min[, "Est."], slopes_max[, "Est."]),
  St.Err. = c(slopes_min[, "St.Err."], slopes_max[, "St.Err."]),
  t_value = c(slopes_min[, "t value"], slopes_max[, "t value"]),
  CI_95_l = -c(slopes_min[, "CI(95%).l"], slopes_max[, "CI(95%).l"]),
  CI_95_u = -c(slopes_min[, "CI(95%).u"], slopes_max[, "CI(95%).u"]),
  Model = rep(c("Min Tephras", "Max Tephras"), each = max_rows)
)

# Print the slopes data frame
print(slopes_df)


#=====================
# 6. Competing models
#=====================

# Fit non-linear least squares models for competing hypotheses 
# (i.e., not the regime-shift model)


# Total sites model
# ------------------
# Define the model function that predicts tephra counts based on the count of all sites
myfunction_all <- function(total_sites, all) {
  (total_sites * all)
}

# Starting values for coefficients
start_params_all <- list(all = 1)

# Fit the NLS model using the Levenberg-Marquardt algorithm: here, for maximum tephras
maxfit_T <- nlsLM(
  tephras_count_2 ~ myfunction_all(total_sites, all),
  data = bins,
  start = start_params_all,
  lower = c(0)
)

# Fit the NLS model using the Levenberg-Marquardt algorithm: here, for minimum tephras
minfit_T <- nlsLM(
  tephras_count_1 ~ myfunction_all(total_sites, all),
  data = bins,
  start = start_params_all,
  lower = c(0)
)


# Site type model
# ----------------
# Define the model function that predicts tephra counts based on
# marine, lake and terrestrial sites

myfunction <- function(site_counts_1,
                       site_counts_2,
                       site_counts_3,
                       Mar,
                       Lake,
                       Terr) {

    site_counts_1 * Mar +
    site_counts_2 * Lake +
    site_counts_3 * Terr
}

# Starting values for coefficients
start_params <- list(Mar = 1, Lake = 1, Terr = 1)

# Fit the NLS model using the Levenberg-Marquardt algorithm: here, for minimum tephras
minfit <- nlsLM(
  tephras_count_1 ~ myfunction(
    site_counts_1,
    site_counts_2,
    site_counts_3,
    Mar,
    Lake,
    Terr
  ),
  data = bins,
  start = start_params,
  lower = c(0, 0, 0)
)

# Fit the NLS model using the Levenberg-Marquardt algorithm: here, for maximum tephras
maxfit <- nlsLM(
  tephras_count_2 ~ myfunction(
    site_counts_1,
    site_counts_2,
    site_counts_3,
    Mar,
    Lake,
    Terr
  ),
  data = bins,
  start = start_params,
  # weights = 1/sqrt(total_sites),
  lower = c(0, 0, 0)
)


# Get predictions from the fitted model
# Also predict confidence intervals

mypreds_max <- predictNLS(maxfit)$summary

mypreds_min <- predictNLS(minfit)$summary

mypreds_max_T <- predictNLS(maxfit_T)$summary

mypreds_min_T <- predictNLS(minfit_T)$summary

# Bind predictions to bins data
bins <- bins %>%
  mutate(
    predicted_max_tephras = mypreds_max$mean.1 / bins$total_sites,
    predicted_max_tephras_lower = mypreds_max$`2.5%` / bins$total_sites,
    predicted_max_tephras_upper = mypreds_max$`97.5%` / bins$total_sites,
    predicted_min_tephras = mypreds_min$mean.1 / bins$total_sites,
    predicted_min_tephras_lower = mypreds_min$`2.5%` / bins$total_sites,
    predicted_min_tephras_upper = mypreds_min$`97.5%` / bins$total_sites,
    predicted_min_tephras_T = mypreds_min_T$mean.1 / bins$total_sites,
    predicted_min_tephras_lower_T = mypreds_min_T$`2.5%` / bins$total_sites,
    predicted_min_tephras_upper_T = mypreds_min_T$`97.5%` / bins$total_sites,
    predicted_max_tephras_T = mypreds_max_T$mean.1 / bins$total_sites,
    predicted_max_tephras_lower_T = mypreds_max_T$`2.5%` / bins$total_sites,
    predicted_max_tephras_upper_T = mypreds_max_T$`97.5%` / bins$total_sites
  )

# Prepare for cumulative plot with type-based predictions

# Calculate cumulative intensity corrected tephra counts
bins <- bins %>%
  arrange(desc(bin_center)) %>%
  mutate(
    cumulative_predicted_min_tephras = cumsum(predicted_min_tephras),
    cumulative_predicted_max_tephras = cumsum(predicted_max_tephras),
    cumulative_predicted_min_tephras_lower = cumsum(predicted_min_tephras_lower),
    cumulative_predicted_min_tephras_upper = cumsum(predicted_min_tephras_upper),
    cumulative_predicted_max_tephras_lower = cumsum(predicted_max_tephras_lower),
    cumulative_predicted_max_tephras_upper = cumsum(predicted_max_tephras_upper),
    cumulative_predicted_min_tephras_T = cumsum(predicted_min_tephras_T),
    cumulative_predicted_max_tephras_T = cumsum(predicted_max_tephras_T),
    cumulative_predicted_min_tephras_lower_T = cumsum(predicted_min_tephras_lower_T),
    cumulative_predicted_min_tephras_upper_T = cumsum(predicted_min_tephras_upper_T),
    cumulative_predicted_max_tephras_lower_T = cumsum(predicted_max_tephras_lower_T),
    cumulative_predicted_max_tephras_upper_T = cumsum(predicted_max_tephras_upper_T)
  )

# Calculate cumulative error bars
bins <- bins %>%
  mutate(
    cumulative_predicted_min_tephras_error = sqrt(cumsum((predicted_min_tephras_upper - predicted_min_tephras_lower) / 2)^2),
    cumulative_predicted_max_tephras_error = sqrt(cumsum((predicted_max_tephras_upper - predicted_max_tephras_lower) / 2)^2),
    cumulative_predicted_min_tephras_error_T = sqrt(cumsum((
      predicted_min_tephras_upper_T - predicted_min_tephras_lower_T
    ) / 2)^2),
    cumulative_predicted_max_tephras_error_T = sqrt(cumsum((
      predicted_max_tephras_upper_T - predicted_max_tephras_lower_T
    ) / 2)^2)
  )


#==================
# 7. Cache results
#==================

cache_file <- paste0("Monteath_analysis_cache", SUFFIX, ".RData")
save(list = ls(), file = cache_file)
message("Saved analysis results to ", cache_file)
