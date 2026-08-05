
#===============================
# ~~~ SPD of tephra deposits ~~~
#===============================

# Calculates summed probability distribution (SPD) of Italian tephra
# deposits in the Mediterranean over the past 25 ka
# Fits continuous piecewise linear (CPL) models to the SPD

# Implements ADMUR package and associated methodology from:

#    Timpson, A., Barberena, R., Thomas, M. G., Méndez, C.,
#    & Manning, K. (2021). Directly modelling population dynamics
#    in the South American Arid Diagonal using 14 C dates. Philosophical
#    Transactions of the Royal Society B: Biological Sciences, 376(1816),
#    20190723. https://doi.org/10.1098/rstb.2019.0723)

#    ADMUR: Ancient Demographic Modelling Using Radiocarbon. Adrian
#    Timpson. University College London. Research Department of Genetics,
#    Environment and Evolution (GEE), Darwin Building, Gower Street,
#    London, WC1E 6BT. 2020 https://github.com/UCL/ADMUR


# Contents:
# ------------
#   1. Packages
#   2. Load and clean data
#   3. Summed probability distribution (SPD)
#   4. Continuous piecewise linear (CPL) models
#   5. MCMC credible interval for a chosen CPL model
#   6. Intensity correction by site availablity
#   6b. Site availability-adjusted CPL fits
#   7. Cumulative SPD (raw, intensity-corrected, and CPL fits)
#   8. Save analysis cache for SPD_figures.R


#==================
# 1. Load packages
#==================

library(readxl)
library(dplyr)
library(ADMUR)
library(DEoptimR)
library(parallel)


#========================
# 2. Load and clean data
#========================

DATA_FILE     <- "Mediterranean_tephra_data.xlsx"
DEPOSIT_SHEET <- "Tephra data"   # The sheet containing the list of tephra deposits
SITE_SHEET    <- "Site spans"    # The sheet containing the list of site spans

MIN_DISTANCE_KM <- 100   # Minimum source-to-deposit distance (km); excludes proximal deposits

# Filter by volcanic source
# "italian" includes only Italian-sourced tephra
# "all" includes all tephras in the dataset regardless of source
ANALYSIS_SCOPE <- "italian"    # Change me: "italian", "all"
ANALYSIS_SCOPE <- match.arg(ANALYSIS_SCOPE, c("italian", "all"))

# Attribution scenario:
# "min" = certain (exact-match) attributions, primary populations only
# "max" = also includes uncertain "(?)" attributions and secondary populations
SOURCE_SCENARIO <- "max"    # Change me: "min", "max"

# Phasing scheme (what ADMUR treats as one independent "phase"):
# "deposit"  = each dated deposit is its own phase (effectively unphased)
# "eruption" = deposits attributed to the same named "Source eruption" are
#              grouped into one phase; deposits with an unknown/blank source
#              eruption are each kept as their own singleton phase, since
#              they are not known to represent the same eruptive event
PHASE_TYPE <- "deposit"    # Change me: "deposit", "eruption"
PHASE_TYPE <- match.arg(PHASE_TYPE, c("deposit", "eruption"))

ITALIAN_SOURCES_EXACT <- c("Campi Flegrei", "Ischia", "Somma-Vesuvius", "Procida-Vivara",
                           "Campanian volcanic field", "Etna", "Lipari", "Vulcano",
                           "Palinuro Seamount", "Aeolian")
ITALIAN_SOURCE_PATTERN <- paste("Campi Flegrei|Ischia|Somma.Vesuvius|Procida|Campanian",
                                "Phlegraean|Etna|Lipari|Vulcano|Aeolian|Palinuro", sep = "|")

SUFFIX <- paste0("_", ANALYSIS_SCOPE, "_", SOURCE_SCENARIO,
                 if (PHASE_TYPE == "eruption") "_by_eruption" else "_by_deposit")

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
  sec_ok <- !is.na(n2) & !is.na(dist_sec)    # Identify rows where a secondary population is present 
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

tephra <- expand_sources(read_excel(DATA_FILE, sheet = DEPOSIT_SHEET)) %>%
  mutate(
    med = suppressWarnings(as.numeric(`Median age (cal yr BP)`)),
    lo  = suppressWarnings(as.numeric(`Minimum age (cal yr BP)`)),
    hi  = suppressWarnings(as.numeric(`Maximum age (cal yr BP)`)),
    distance_km = suppressWarnings(as.numeric(`Great circle distance between source and deposit (km)`))
  ) %>%
  filter(!is.na(med), !is.na(lo), !is.na(hi), !is.na(distance_km)) %>%   # Essential fields present
  filter(distance_km >= MIN_DISTANCE_KM) %>%                            # Exclude proximal deposits
  filter(
    if (ANALYSIS_SCOPE == "italian") {
      if (SOURCE_SCENARIO == "min") {
        `Source volcano` %in% ITALIAN_SOURCES_EXACT & source_rank == "primary"
      } else {
        grepl(ITALIAN_SOURCE_PATTERN, `Source volcano`, ignore.case = TRUE)
      }
    } else {
      if (SOURCE_SCENARIO == "min") {
        source_rank == "primary"
      } else {
        TRUE
      }
    }
  ) %>%
  mutate(
    # Sd from 95% (2-sigma) calibrated range 
    # Precisely-dated eruptions (e.g., 79 CE Vesuvius) are given sd = 15 yr
    # to avoid dividing by zero in downstream calculations
    sd = ifelse((hi - lo) > 0, (hi - lo) / 3.92, 0),
    sd = pmax(sd, 15)         # No tephra treated as more precise than +/-15 yr
  )


#============================================
# 3. Summed probabability distribution (SPD)
#============================================

# Assign each deposit to a phase, according to PHASE_TYPE
# "deposit": every deposit is its own phase (essentially unphased)
# "eruption": deposits sharing a named 'Source eruption' are grouped into one
#             phase; unknown/blank eruption deposits each get a unique
#             placeholder phase so unrelated eruptions are never merged
build_phase <- function(tephra, phase_type) {
  if (phase_type == "deposit") return(seq_len(nrow(tephra)))

  erup    <- trimws(as.character(tephra$`Source eruption`))
  unknown <- is.na(erup) | erup == "" | erup == "Unknown"
  phase   <- erup
  phase[unknown] <- paste0(".unknown.", seq_len(nrow(tephra))[unknown])
  phase
}

# Build the ADMUR data frame
# Contain the dates to be calibrated and summed
# Requires 'age' and standard deviation, 'sd'
data <- data.frame(
  age        = tephra$med,
  sd         = tephra$sd,
  phase      = build_phase(tephra, PHASE_TYPE),
  datingType = "calendar"    # Non-14C dates are treated as calendar dates
)                            # anything != '14C' is used directly, no calibration

# Construct calibration curve
# Note that this step is needed even when 14C dates are not used
age_range <- c(0, 25000)    # The date range for which to construct calibration curve
grid_inc  <- 10    # The resolution of the curve, in this case set as 10 yrs

CalArray <- makeCalArray(intcal20, calrange = age_range, inc = grid_inc)
# Generates an array of probabilities mapping the calibration curve and its error ribbon

cal <- summedCalibrator(data, CalArray, normalise = "none")
# Calibrates radiocarbon dates through CalArray then sums all probability distributions
# When normalise = "none", the output SPD, the area of the curve is equal to number of samples
# When normalise = "standard", the output SPD is normalised so the curve integrates to 1

years        <- as.numeric(rownames(cal))
spd_vals_raw <- cal[, 1]   

PD <- phaseCalibrator(data, CalArray, remove.external = TRUE)

spd_vals <- rowSums(PD)


#=======================================
# 4. Continuous Piecewise Linear models
#=======================================

# Optional taphonomic correction. Adds ADMUR's 'power' model as a second
# component alongside the CPL model
APPLY_TAPHONOMY_CORRECTION <- FALSE   # Set TRUE to include taphonomic correction

CPL_ORDERS <- 1:5   # Change to try more/fewer model orders (e.g. 1:5, 1:10)

# Fit CPL models with 1 to 5 pieces to the SPD, using multi-restart optimisation
fit_cpl_batch <- function(n_pieces_list, PD, years = 0:25000,
                           n_restarts = 5,    # Always runs 5 restarts per model
                           n_cores = max(1, detectCores() - 1),
                           taphonomy = APPLY_TAPHONOMY_CORRECTION) {

  model_type <- if (taphonomy) c("CPL", "power") else "CPL"
  jobs <- expand.grid(n_pieces = n_pieces_list, seed = seq_len(n_restarts))

  run_one <- function(n_pieces, seed) {
    n_par_cpl <- 2 * n_pieces - 1

    if (taphonomy) {
      lower <- c(rep(0, n_par_cpl), 0, -3)
      upper <- c(rep(1, n_par_cpl), max(years), 0)
    } else {
      lower <- rep(0, n_par_cpl)
      upper <- rep(1, n_par_cpl)
    }

    n_par   <- length(lower)
    NP      <- 10 * n_par
    maxiter <- max(500, 200 * n_par)

    set.seed(seed)
    tryCatch(
      JDEoptim(lower   = lower,
               upper   = upper,
               fn      = objectiveFunction,
               PDarray = PD,
               type    = model_type,
               NP      = NP,
               maxiter = maxiter,
               tol     = 1e-6,
               trace   = FALSE),
      error = function(e) NULL
    )
  }

  results <- mcmapply(run_one, jobs$n_pieces, jobs$seed,
                       mc.cores = n_cores, SIMPLIFY = FALSE)

  # Group restarts back by model order and keep the best of each
  fits <- lapply(n_pieces_list, function(n_pieces) {
    runs <- Filter(Negate(is.null), results[jobs$n_pieces == n_pieces])
    if (length(runs) == 0) stop("All restarts failed for ", n_pieces, "-CPL")

    values <- vapply(runs, `[[`, numeric(1), "value")
    best   <- runs[[which.min(values)]]
r
    hinges <- convertPars(pars = best$par, years = years, type = model_type)

    list(n_pieces        = n_pieces,
         type            = model_type,
         par             = best$par,
         neglogL         = best$value,
         hinges          = hinges,
         n_restarts_used = length(runs),
         neglogL_range   = diff(range(values)))
  })
  setNames(fits, paste0("cpl", n_pieces_list))
}

fits_file <- paste0("fits", SUFFIX, ".rds")

if (file.exists(fits_file)) {
  fits <- readRDS(fits_file)
  message("Loaded saved CPL fits from ", fits_file)
} else {
  message("Fitting CPL models (", min(CPL_ORDERS), " to ", max(CPL_ORDERS), " pieces, parallel multi-start",
          if (APPLY_TAPHONOMY_CORRECTION) ", with taphonomic correction" else "", ")...")
  t0 <- Sys.time()
  fits <- fit_cpl_batch(CPL_ORDERS, PD = PD)
  for (f in fits) {
    message(f$n_pieces, "-CPL: neglogL = ", round(f$neglogL, 2),
            " (", f$n_restarts_used, " restart",
            if (f$n_restarts_used != 1) "s" else "",
            " used, range = ", round(f$neglogL_range, 3), ")")
  }
  message("CPL fitting took ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  saveRDS(fits, fits_file)
  message("Saved CPL fits to ", fits_file)
}

# Combine all CPLs into one data frame
all_cpl_df <- do.call(rbind, lapply(fits, function(f) {
  data.frame(age      = f$hinges$year,
             density  = f$hinges$pdf,
             model    = paste0(f$n_pieces, "-CPL"))
}))

# BIC comparison across CPL orders
# Lower BIC is preferred
N_phases <- ncol(PD)
cpl_bic_table <- do.call(rbind, lapply(fits, function(f) {
  data.frame(model         = paste0(f$n_pieces, "-CPL"),
             n_pieces      = f$n_pieces,
             K             = length(f$par),
             neglogL       = f$neglogL,
             BIC           = log(N_phases) * length(f$par) + 2 * f$neglogL,
             neglogL_range = f$neglogL_range)
}))
cpl_bic_table <- cpl_bic_table[order(cpl_bic_table$BIC), ]
print(cpl_bic_table)

CPL_CONVERGENCE_TOL <- 1.0
converged_table <- cpl_bic_table[cpl_bic_table$neglogL_range <= CPL_CONVERGENCE_TOL, ]
if (nrow(converged_table) == 0) {
  warning("No CPL model converged reliably (all neglogL_range > ", CPL_CONVERGENCE_TOL,
          "); falling back to the full BIC table")
  converged_table <- cpl_bic_table
}

best_cpl_model <- converged_table$model[1]
cat("\nBest CPL model by BIC (restricted to reliably-converged fits): ", best_cpl_model, "\n", sep = "")
if (!identical(best_cpl_model, cpl_bic_table$model[1])) {
  cat("Note: ", cpl_bic_table$model[1], " has lower BIC but did not converge reliably ",
      "(neglogL_range = ", round(cpl_bic_table$neglogL_range[1], 3), "); excluded.\n", sep = "")
}

# Use the best BIC-selected CPL model throughout rest of the analysis
BEST_CPL_KEY <- paste0("cpl", converged_table$n_pieces[1])

# Goodness of fit (GOF) test
GOF_N_SIM <- 10000

GOF_CALRANGE <- c(200, age_range[2])

gof_file <- paste0("gof_", BEST_CPL_KEY, SUFFIX, ".rds")

if (file.exists(gof_file)) {
  gof <- readRDS(gof_file)
  message("Loaded saved GOF test from ", gof_file)
} else {
  message("Running GOF simulation test for ", BEST_CPL_KEY,
          " (N = ", GOF_N_SIM, " simulations, this can take a while)...")
  t0 <- Sys.time()
  gof <- SPDsimulationTest(data, calcurve = intcal20, calrange = GOF_CALRANGE,
                            pars = fits[[BEST_CPL_KEY]]$par,
                            type = fits[[BEST_CPL_KEY]]$type,
                            N    = GOF_N_SIM)
  message("GOF test took ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  saveRDS(gof, gof_file)
  message("Saved GOF test to ", gof_file)
}
cat("GOF p-value (", BEST_CPL_KEY, "): ", gof$pvalue, "\n", sep = "")

#===============================================
# 5. MCMC credible interval for the best CPL model
#===============================================

RUN_MCMC <- TRUE

# Highlight the BIC-selected model in the single-model figure
MCMC_MODEL <- BEST_CPL_KEY

mcmc_years <- seq(0, 25000, by = 100)  

CPL_JUMPS <- c(cpl1 = 0.7, cpl2 = 0.05, cpl3 = 0.01, cpl4 = 0.006,
               cpl5 = if (SOURCE_SCENARIO == "max") 0.005 else 0.0007)

tune_jump <- function(start_jump, startPars, type, target = 0.45,
                       tol = 0.1, pilot_N = 3000, max_tries = 6) {
  jump <- start_jump
  for (i in seq_len(max_tries)) {
    pilot <- mcmc(PDarray = PD, startPars = startPars, type = type,
                  N = pilot_N, burn = floor(pilot_N * 0.1), thin = 1, jumps = jump)
    ar <- pilot$acceptance.ratio
    message("  jump tuning attempt ", i, ": jump = ", signif(jump, 3),
            ", acceptance ratio = ", round(ar, 3))
    if (abs(ar - target) <= tol) return(jump)
    jump <- jump * (ar / target)
  }
  message("  jump tuning did not fully converge after ", max_tries,
          " tries; using last value ", signif(jump, 3))
  jump
}

# Run (or load cached) MCMC for one model, and derive its credible interval
run_mcmc <- function(model_name, max_attempts = 3) {
  model     <- fits[[model_name]]
  mcmc_file <- paste0("mcmc_", model_name, SUFFIX, ".rds")

  if (file.exists(mcmc_file)) {
    chain <- readRDS(mcmc_file)
    message("Loaded saved MCMC chain from ", mcmc_file)
  } else {
    message("Tuning jump size for ", model_name, "...")
    start_jump <- if (!is.na(CPL_JUMPS[[model_name]])) CPL_JUMPS[[model_name]] else 0.05
    jump <- tune_jump(start_jump, startPars = model$par, type = model$type)

    message("Running MCMC for ", model_name, " (this can take up to an hour)...")
    t0 <- Sys.time()
    chain <- NULL
    attempt <- 1
    while (is.null(chain) && attempt <= max_attempts) {
      chain <- tryCatch(
        mcmc(PDarray   = PD,
             startPars = model$par,
             type      = model$type,
             jumps     = jump),
        error = function(e) {
          message(model_name, ": attempt ", attempt, " failed (", conditionMessage(e), "), retrying...")
          NULL
        }
      )
      attempt <- attempt + 1
    }
    if (is.null(chain)) stop("MCMC failed for ", model_name, " after ", max_attempts, " attempts")
    message(model_name, ": MCMC took ",
            round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min",
            ", acceptance ratio = ", round(chain$acceptance.ratio, 2))
    saveRDS(chain, mcmc_file)
    message("Saved MCMC chain to ", mcmc_file)
  }

  res_mat <- if (is.matrix(chain$res)) chain$res else matrix(chain$res, ncol = 1)
  draws   <- convertPars(pars = res_mat, years = mcmc_years, type = model$type)
  pdf_mat <- as.matrix(draws[, -1])
  ci <- data.frame(
    age    = draws$year,
    lower  = apply(pdf_mat, 1, quantile, probs = 0.025),
    median = apply(pdf_mat, 1, quantile, probs = 0.5),
    upper  = apply(pdf_mat, 1, quantile, probs = 0.975)
  )

  list(chain = chain, ci = ci)
}

MCMC_MODELS <- BEST_CPL_KEY

if (RUN_MCMC) {
  mcmc_results <- mclapply(setNames(MCMC_MODELS, MCMC_MODELS), run_mcmc,
                            mc.cores = max(1, detectCores() - 1))
  mcmc_model <- fits[[MCMC_MODEL]]
  chain      <- mcmc_results[[MCMC_MODEL]]$chain
  mcmc_ci    <- mcmc_results[[MCMC_MODEL]]$ci
}


#=========================
# 6. Intensity-correction
#=========================

# Normalise by site coverage 
# Intensity-corrected SPD = raw SPD / N(t)

# Load site-span data
sites <- read_excel(DATA_FILE, sheet = SITE_SHEET) %>%
  mutate(
    `Minimum site coverage (cal yr BP)` = suppressWarnings(
      as.numeric(`Minimum site coverage (cal yr BP)`)
    ),
    `Maximum site coverage (cal yr BP)` = suppressWarnings(
      as.numeric(`Maximum site coverage (cal yr BP)`)
    )
  ) %>%
  filter(!is.na(`Minimum site coverage (cal yr BP)`),
         !is.na(`Maximum site coverage (cal yr BP)`))

# Compute number of sites with coverage at year t, N(t)
# Follows the same logic as the count_sites_in_bin() function in Monteath_adapted.R
n_sites_t <- colSums(
  outer(sites$`Minimum site coverage (cal yr BP)`, years, "<=") &
    outer(sites$`Maximum site coverage (cal yr BP)`, years, ">=")
)

# Build a tidy data frame combining raw SPD, N(t), and intensity-corrected SPD
spd_df <- data.frame(
  age      = years,
  spd      = spd_vals,
  n_sites  = n_sites_t
) %>%
  mutate(
    spd_intensity = ifelse(n_sites > 0, spd / n_sites, NA_real_)   
  )


#=========================================
# 6b. Site-availability-adjusted CPL fits
#=========================================

siteavail_fits_file <- paste0("fits", SUFFIX, "_siteavail_pilot.rds")
HAS_SITEAVAIL <- file.exists(siteavail_fits_file)

if (HAS_SITEAVAIL) {
  fits_siteavail <- readRDS(siteavail_fits_file)

  siteavail_bic_table <- do.call(rbind, lapply(fits_siteavail, function(f) {
    data.frame(model         = paste0(f$n_pieces, "-CPL+ts"),
               n_pieces      = f$n_pieces,
               K             = length(f$par),
               r_fitted      = f$par[length(f$par)],
               neglogL       = f$neglogL,
               BIC           = log(N_phases) * length(f$par) + 2 * f$neglogL,
               neglogL_range = f$neglogL_range)
  }))
  siteavail_bic_table <- siteavail_bic_table[order(siteavail_bic_table$BIC), ]
  print(siteavail_bic_table)

  converged_siteavail_table <- siteavail_bic_table[siteavail_bic_table$neglogL_range <= CPL_CONVERGENCE_TOL, ]
  if (nrow(converged_siteavail_table) == 0) {
    warning("No CPL+timeseries model converged reliably (all neglogL_range > ", CPL_CONVERGENCE_TOL,
            "); falling back to the full BIC table")
    converged_siteavail_table <- siteavail_bic_table
  }
  BEST_SITEAVAIL_KEY <- paste0("cpl", converged_siteavail_table$n_pieces[1])
  cat("\nBest CPL+timeseries model by BIC (restricted to reliably-converged fits): ",
      BEST_SITEAVAIL_KEY, "\n", sep = "")
  if (!identical(BEST_SITEAVAIL_KEY, paste0("cpl", siteavail_bic_table$n_pieces[1]))) {
    cat("Note: ", siteavail_bic_table$model[1], " has lower BIC but did not converge reliably ",
        "(neglogL_range = ", round(siteavail_bic_table$neglogL_range[1], 3), "); excluded.\n", sep = "")
  }

  timeseries_df <- data.frame(x = years, y = n_sites_t)

  all_cpl_ts_df <- do.call(rbind, lapply(converged_siteavail_table$n_pieces, function(n_pieces) {
    f          <- fits_siteavail[[paste0("cpl", n_pieces)]]
    n_cpl_pars <- 2 * n_pieces - 1
    curve      <- convertPars(pars = f$par[1:n_cpl_pars], years = years, type = "CPL")
    data.frame(age = curve$year, density = curve$pdf, model = paste0(n_pieces, "-CPL+ts"))
  }))

  siteavail_mcmc_file <- paste0("mcmc_", BEST_SITEAVAIL_KEY, SUFFIX, "_siteavail.rds")
  if (file.exists(siteavail_mcmc_file)) {
    chain_siteavail   <- readRDS(siteavail_mcmc_file)
    best_n_cpl_pars   <- 2 * as.integer(sub("cpl", "", BEST_SITEAVAIL_KEY)) - 1
    chain_res_cpl     <- chain_siteavail$res[, 1:best_n_cpl_pars, drop = FALSE]
    draws_siteavail   <- convertPars(pars = chain_res_cpl, years = mcmc_years, type = "CPL")
    pdf_mat_siteavail <- as.matrix(draws_siteavail[, -1])
    mcmc_ci_siteavail <- data.frame(
      age    = draws_siteavail$year,
      lower  = apply(pdf_mat_siteavail, 1, quantile, probs = 0.025),
      median = apply(pdf_mat_siteavail, 1, quantile, probs = 0.5),
      upper  = apply(pdf_mat_siteavail, 1, quantile, probs = 0.975)
    )
  } else {
    chain_siteavail   <- NULL
    mcmc_ci_siteavail <- NULL
  }

  siteavail_gof_file <- paste0("gof_", BEST_SITEAVAIL_KEY, SUFFIX, "_siteavail.rds")
  gof_siteavail <- if (file.exists(siteavail_gof_file)) readRDS(siteavail_gof_file) else NULL

  cat("GOF p-value (", BEST_SITEAVAIL_KEY, "+ts): ",
      if (!is.null(gof_siteavail)) gof_siteavail$pvalue else NA, "\n", sep = "")
} else {
  siteavail_bic_table <- NULL
  all_cpl_ts_df        <- NULL
  mcmc_ci_siteavail     <- NULL
  BEST_SITEAVAIL_KEY    <- NA_character_
  gof_siteavail         <- NULL
}


#==============================================================
# 7. Cumulative raw SPD, intensity-corrected SPD, and CPL fits
#==============================================================

# Cumulate the raw and intensity-corrected SPD (oldest to youngest)
spd_df_ord <- spd_df %>%
  arrange(desc(age)) %>%
  mutate(
    cum_spd_raw       = cumsum(spd) * grid_inc,
    cum_spd_intensity = cumsum(replace(spd_intensity, is.na(spd_intensity), 0)) * grid_inc
  )

# Intensity-corrected cumulative CPL curves
# Re-evaluates already-optimised parameters (no re-fitting) on the same 10-yr grid 
# as n_sites_t, so it can be divided by site coverage the same way the raw SPD was above
N_total <- max(spd_df_ord$cum_spd_raw)  

# Builds the cumulative-overlay curve and hinge (breakpoint) credible
# intervals for one CPL model key (e.g. "cpl4"), reusing its already-fitted
# parameters and already-run MCMC chain 
build_cpl_overlay <- function(key) {
  hinges_grid <- convertPars(pars = fits[[key]]$par, years = years, type = "CPL")
  pdf       <- hinges_grid$pdf * N_total
  intensity <- ifelse(n_sites_t > 0, pdf / n_sites_t, NA_real_)

  cum_df <- data.frame(age = years, pdf = pdf, intensity = intensity) %>%
    arrange(desc(age)) %>%
    mutate(cum_int = cumsum(replace(intensity, is.na(intensity), 0)) * grid_inc)

  hinges_pt <- CPLparsToHinges(fits[[key]]$par, years = 0:25000)

  chain_res   <- mcmc_results[[key]]$chain$res
  chain_res   <- if (is.matrix(chain_res)) chain_res else matrix(chain_res, ncol = 1)
  hinges_post <- CPLparsToHinges(chain_res, years = 0:25000)

  interior <- if (nrow(hinges_pt) > 2) 2:(nrow(hinges_pt) - 1) else integer(0)
  yr_cols  <- paste0("yr", interior)

  hinge_df <- if (length(interior) == 0) {
    data.frame(age = numeric(0), xmin = numeric(0), xmax = numeric(0))
  } else {
    data.frame(
      age  = hinges_pt$year[interior],
      xmin = vapply(yr_cols, function(col) quantile(hinges_post[[col]], 0.025), numeric(1)),
      xmax = vapply(yr_cols, function(col) quantile(hinges_post[[col]], 0.975), numeric(1))
    )
  }
  hinge_df$cum_int <- approx(x = cum_df$age, y = cum_df$cum_int, xout = hinge_df$age)$y

  list(hinges_grid = hinges_grid, pdf = pdf, intensity = intensity,
       cum_df = cum_df, hinges_pt = hinges_pt, hinges_post = hinges_post, hinge_df = hinge_df)
}

bestcpl_overlay     <- build_cpl_overlay(BEST_CPL_KEY)
bestcpl_hinges_grid <- bestcpl_overlay$hinges_grid
bestcpl_pdf         <- bestcpl_overlay$pdf
bestcpl_intensity   <- bestcpl_overlay$intensity
bestcpl_cum_df      <- bestcpl_overlay$cum_df
bestcpl_hinges_pt   <- bestcpl_overlay$hinges_pt
bestcpl_hinges_post <- bestcpl_overlay$hinges_post
bestcpl_hinge_df    <- bestcpl_overlay$hinge_df

build_cplts_overlay <- function(key) {
  f          <- fits_siteavail[[key]]
  n_cpl_pars <- 2 * f$n_pieces - 1
  cpl_par    <- f$par[1:n_cpl_pars]

  hinges_grid <- convertPars(pars = cpl_par, years = years, type = "CPL")
  intensity   <- hinges_grid$pdf * N_total / mean(n_sites_t)

  cum_df <- data.frame(age = years, pdf = hinges_grid$pdf, intensity = intensity) %>%
    arrange(desc(age)) %>%
    mutate(cum_int = cumsum(intensity) * grid_inc)

  hinges_pt <- CPLparsToHinges(cpl_par, years = 0:25000)
  interior  <- if (nrow(hinges_pt) > 2) 2:(nrow(hinges_pt) - 1) else integer(0)

  hinge_df <- data.frame(age = numeric(0), xmin = numeric(0), xmax = numeric(0))
  if (length(interior) > 0 && !is.null(chain_siteavail)) {
    chain_res_cpl <- chain_siteavail$res[, 1:n_cpl_pars, drop = FALSE]
    hinges_post   <- CPLparsToHinges(chain_res_cpl, years = 0:25000)
    yr_cols       <- paste0("yr", interior)
    hinge_df <- data.frame(
      age  = hinges_pt$year[interior],
      xmin = vapply(yr_cols, function(col) quantile(hinges_post[[col]], 0.025), numeric(1)),
      xmax = vapply(yr_cols, function(col) quantile(hinges_post[[col]], 0.975), numeric(1))
    )
  }
  hinge_df$cum_int <- approx(x = cum_df$age, y = cum_df$cum_int, xout = hinge_df$age)$y

  list(pdf = hinges_grid$pdf, intensity = intensity, cum_df = cum_df, hinge_df = hinge_df)
}

if (HAS_SITEAVAIL) {
  bestcplts_overlay   <- build_cplts_overlay(BEST_SITEAVAIL_KEY)
  bestcplts_pdf       <- bestcplts_overlay$pdf
  bestcplts_intensity <- bestcplts_overlay$intensity
  bestcplts_cum_df    <- bestcplts_overlay$cum_df
  bestcplts_hinge_df  <- bestcplts_overlay$hinge_df
} else {
  bestcplts_cum_df   <- NULL
  bestcplts_hinge_df <- NULL
}


#=====================================
# 8. Save analysis cache for plotting
#=====================================

analysis_cache_file <- paste0("SPD_analysis_cache", SUFFIX, ".RData")
save(list = ls(), file = analysis_cache_file)
message("Saved analysis cache to ", analysis_cache_file)

# Plotting occurs in SPD_figures.R
