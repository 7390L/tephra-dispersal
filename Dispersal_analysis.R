
#============================================
# ~~~ Tephra Transport Direction Analysis ~~~
#============================================

# Analyses tephra dispersal using circular statistics

# Contents:
# ------------
#   1. Packages
#   2. Configuration
#   3. Volcano grouping schemes
#   4. Age binning schemes
#   5. Circular statistics: Tests of uniformity ~ Rayleigh test and Rao's Spacings test
#   6. Circular statistics: Watson U2 test
#   7. Permutation null test
#   8. Residual bearings 
#   9. Monte Carlo age uncertainty 
#   10. Eruption collapse     
#    
#   Analysis pipeline:
#   11. Load and clean data
#   12. Assign groups and bins
#   13. Run statistical tests
#   14. Write CSVs
#   15. Save results (analysis_results.rds)


#==================
# 1. Load packages
#==================

library(readxl)
library(dplyr)
library(circular)


#==================
# 2. Configuration
#==================

DATA_FILE  <- "Mediterranean_tephra_data.xlsx"
DATA_SHEET <- "Tephra data"    # The sheet containing the list of tephra deposits

MIN_N <- 5    # Minimum deposits per bin for a result to be reported

MIN_DISTANCE_KM <- 100    # Minimum source-to-deposit distance (km); 
                          # More proximal deposits are excluded

# Monte Carlo sampling for tephra age uncertainty
N_MC    <- 1000   # Number of Monte Carlo iterations
MC_SEED <- 42     # Random seed for reproducibility

# Permutation null test
N_PERM  <- 1000    # Number of permutations for Watson U2 temporal null test

# Residual bearing analysis
MIN_VOLCANO_N <- 5   # Minimum tephra deposits per volcano to compute a reliable mean bearing


#=====================
# 3. Volcano grouping
#=====================

# Select a volcano grouping scheme by defining GROUPING_SCHEME

# "refined"  = three geographically distinct groups: Campanian, Etna, Aeolian
# "broad"    = two groups: Campanian, all other Italian volcanoes
# "italian"      = all Italian volcanoes pooled together
# "mediterranean" = all Mediterranean volcanoes pooled together (includes Santorini)
# or define a new grouping scheme

GROUPING_SCHEME <- "refined"   # Change me: "refined", "broad", "italian", "mediterranean"

GROUPING_SCHEMES <- list(
  
  refined = list(
    groups = list(
      Campanian = list(
        exact   = c("Campi Flegrei", "Ischia", "Somma-Vesuvius",    # Dataset entries must match exactly (used for minimum scenario)
                    "Procida-Vivara", "Campanian volcanic field"),  
        pattern = "Campi Flegrei|Ischia|Somma.Vesuvius|Procida|Campanian|Phlegraean"    # Includes uncertain attributions demarked by (?) (used for maximum scenario) 
      ),
      Etna = list(
        exact   = c("Etna"),
        pattern = "Etna"
      ),
      Aeolian = list(
        exact   = c("Lipari", "Vulcano", "Palinuro Seamount", "Aeolian"),
        pattern = "Lipari|Vulcano|Aeolian|Palinuro"
      )
    )
  ),
  
  broad = list(
    groups = list(
      Campanian = list(
        exact   = c("Campi Flegrei", "Ischia", "Somma-Vesuvius",
                    "Procida-Vivara", "Campanian volcanic field"),
        pattern = "Campi Flegrei|Ischia|Somma.Vesuvius|Procida|Campanian|Phlegraean"
      ),
      "Non-Campanian" = list(
        exact   = c("Etna", "Lipari", "Vulcano", "Palinuro Seamount", "Aeolian"),
        pattern = "Etna|Lipari|Vulcano|Aeolian|Palinuro"
      )
    )
  ),
  
  italian = list(
    groups = list(
      "All Italian" = list(
        exact   = c("Campi Flegrei", "Ischia", "Somma-Vesuvius",
                    "Procida-Vivara", "Campanian volcanic field",
                    "Etna", "Lipari", "Vulcano", "Palinuro Seamount", "Aeolian"),
        pattern = paste("Campi Flegrei|Ischia|Somma.Vesuvius|Procida|Campanian",
                        "Phlegraean|Etna|Lipari|Vulcano|Aeolian|Palinuro", sep = "|")
      )
    )
  ),

  mediterranean = list(
    groups = list(
      "All Mediterranean" = list(
        exact   = c("Campi Flegrei", "Ischia", "Somma-Vesuvius",
                    "Procida-Vivara", "Campanian volcanic field",
                    "Etna", "Lipari", "Vulcano", "Palinuro Seamount", "Aeolian",
                    "Santorini"),
        pattern = paste("Campi Flegrei|Ischia|Somma.Vesuvius|Procida|Campanian",
                        "Phlegraean|Etna|Lipari|Vulcano|Aeolian|Palinuro|Santorini", sep = "|")
      )
    )
  )

)

# Activate the selected grouping scheme
VOLCANO_GROUPS <- GROUPING_SCHEMES[[GROUPING_SCHEME]]$groups

# Function to assign each tephra deposit in the dataset to a volcano group
assign_group <- function(volcano_names, scenario = c("min", "max")) {
  scenario <- match.arg(scenario)
  vapply(volcano_names, function(v) {
    for (gname in names(VOLCANO_GROUPS)) {
      g <- VOLCANO_GROUPS[[gname]]
      matched <- if (scenario == "min") {
        v %in% g$exact
      } else {
        grepl(g$pattern, v, ignore.case = TRUE)
      }
      if (matched) return(gname)
    }
    return(NA_character_)
  }, character(1))
}


#================
# 4. Age binning
#================

# Select a binning scheme by defining BINNING_SCHEME

BINNING_SCHEME <- "holocene_onset"   # Change me: "climate", "deglac_fine", "holocene_onset", "fixed_width"

BINNING_SCHEMES <- list(
  
  # Bins defined by major climatic phases
  climate = list(
    breaks = c(25000, 19000, 11700, 8200, 4200, 0),
    labels = c("Last Glacial Maximum (25\u201319 ka)", "Deglaciation (19\u201311.7 ka)",
               "Early Holocene (11.7\u20138.2 ka)", "Mid Holocene (8.2\u20134.2 ka)",
               "Late Holocene (4.2\u20130 ka)")
  ),

  # 1-kyr resolution across the ~14.5 ka jetstream transition (North American
  # saddle collapse hypothesis)
  # Deposits outside the window get dropped
  deglac_fine = list(
    from = 16500,   # Oldest edge (cal yr BP)
    to   = 13500,   # Youngest edge
    by   = 1000     # Bin width
  ),

  # Younger Dryas vs full Early Holocene, testing the onset of the
  # Holocene (11700 cal BP) as a key transition point
  holocene_onset = list(
    breaks = c(12900, 11700, 8200),
    labels = c("Younger Dryas (12.9–11.7 ka)", "Early Holocene (11.7–8.2 ka)")
  ),

  # Even-width bins over the whole record
  fixed_width = list(
    from = 25000,
    to   = 0,
    by   = 2500   # Bin width - change for finer or coarser bins
  )
)

# Function to build equally-spaced breaks and labels for a width-based scheme
make_width_bins <- function(from, to, by) {
  edges <- seq(from, to, by = -abs(by))
  labs  <- paste0(edges[-length(edges)] / 1000, "-", edges[-1] / 1000, " ka")
  list(breaks = edges, labels = labs)
}

# Activate the selected scheme
local({
  sc <- BINNING_SCHEMES[[BINNING_SCHEME]]
  if (!is.null(sc$from)) {
    wb <- make_width_bins(sc$from, sc$to, sc$by)
    BIN_BREAKS <<- wb$breaks
    BIN_LABELS <<- wb$labels
  } else {
    BIN_BREAKS <<- sc$breaks
    BIN_LABELS <<- sc$labels
  }
})

# Function to assign each tephra deposit to a bin based on its median age 
assign_bins <- function(data, age_col = "age_median") {
  ages       <- data[[age_col]]
  asc_breaks <- rev(BIN_BREAKS)   # ascending edges for cut()
  asc_labels <- rev(BIN_LABELS)
  
  data$bin_label <- cut(
    ages,
    breaks         = asc_breaks,
    labels         = asc_labels,
    include.lowest = TRUE,
    right          = FALSE
  )
  midpoints        <- (BIN_BREAKS[-length(BIN_BREAKS)] + BIN_BREAKS[-1]) / 2
  names(midpoints) <- BIN_LABELS
  data$bin_mid     <- midpoints[as.character(data$bin_label)]
  data
}


#==========================================
# 5. Circular statistics: Uniformity tests
#==========================================

# Rayleigh test
# ----------------

# Function to compute Rayleigh test statistics for a vector of bearings in degrees
rayleigh_stats <- function(bearings_deg) {
  b <- bearings_deg[!is.na(bearings_deg)]
  n <- length(b)
  if (n < 3) return(list(n = n, mean_dir = NA_real_, R = NA_real_, p_rayleigh = NA_real_))    # Needs at least 3 bearings to compute Rayleigh test

  rad <- b * pi / 180    # Convert degrees to radians
  
  # Each bearing is treated as a unit vector on the unit circle
  C   <- mean(cos(rad))    # Mean of the x (cosine) components across all n bearings
  S   <- mean(sin(rad))    # Mean of the y (sine) components across all n bearings
  R   <- sqrt(C^2 + S^2)   # Compute length of the mean resultant vector, R 
                           # High R indicates bearings are strongly clustered

  mean_dir <- (atan2(S, C) * 180 / pi + 360) %% 360    # Calculate mean direction and convert back to degrees

  Z <- n * R^2    # Calculate the p-value
  p <- exp(-Z) * (1 + (2*Z - Z^2) / (4*n) -
                    (24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4) / (288*n^2))
  p <- max(0, min(1, p))
  
  # Return a list containing number of bearings, mean direction, R, and p-value
  list(n = n, mean_dir = mean_dir, R = R, p_rayleigh = p)
}


# Rao's spacing test
# ---------------------

# U statistic is computed by rao.spacing.test() from the circular package

rao_stats <- function(bearings_deg) {
  b <- bearings_deg[!is.na(bearings_deg)]
  n <- length(b)
  if (n < 4) return(list(n = n, rao_U = NA_real_, p_rao = NA_real_))    # Needs at least 4 bearings to compute Rao's U statistic

  cb  <- circular(b, units = "degrees", template = "geographics")
  res <- rao.spacing.test(cb)
  U   <- res$statistic

  crit <- getFromNamespace("rao.table", "circular")    # Get the critical values table

  # Select row of critical values table based on sample size n
  table.row <- if      (n <= 30)  n - 3
               else if (n <= 32)  27
               else if (n <= 37)  28
               else if (n <= 42)  29
               else if (n <= 47)  30
               else if (n <= 62)  31
               else if (n <= 87)  32
               else if (n <= 125) 33
               else               34

  # Calculate p-value based on U statistic and critical values for the sample size n
  # Supports four significance levels: 0.001, 0.01, 0.05, 0.10
  p_rao <- if      (U > crit[table.row, 1]) 0.001
            else if (U > crit[table.row, 2]) 0.01
            else if (U > crit[table.row, 3]) 0.05
            else if (U > crit[table.row, 4]) 0.10
            else                             1.00

  list(n = n, rao_U = round(U, 2), p_rao = p_rao)
}


# Summarise results of uniformity tests
# ----------------------------------------

# Function to summarise the results of the Rayleigh and Rao's spacing tests
# Returns a data frame with one row per group x bin combination
# Contains the number of bearings, mean direction, R statistic, U statistic, p-value, and significance

summarise_by_group <- function(data, scenario_name, bearing_col = "bearing") {
  all_bins   <- levels(data$bin_label)
  all_groups <- names(VOLCANO_GROUPS)

  rows <- lapply(all_groups, function(grp) {
    d_grp <- data[data$group == grp, ]
    lapply(all_bins, function(b) {
      sub <- d_grp[!is.na(d_grp$bin_label) & d_grp$bin_label == b, ]
      cs  <- rayleigh_stats(sub[[bearing_col]])
      rs  <- rao_stats(sub[[bearing_col]])
      data.frame(
        scenario        = scenario_name,
        group           = grp,
        bin             = b,
        bin_mid         = if (nrow(sub) > 0) mean(sub$bin_mid, na.rm = TRUE) else NA_real_,
        n               = cs$n,
        mean_dir        = round(cs$mean_dir, 1),
        R               = round(cs$R, 3),
        p_rayleigh      = round(cs$p_rayleigh, 4),
        significant     = !is.na(cs$p_rayleigh) && cs$p_rayleigh < 0.05,
        rao_U           = rs$rao_U,
        p_rao           = rs$p_rao,
        significant_rao = !is.na(rs$p_rao) && rs$p_rao <= 0.05,
        mean_dist       = round(mean(sub$distance_km, na.rm = TRUE), 0),
        low_n           = cs$n < MIN_N,
        stringsAsFactors = FALSE
      )
    })
  })
  do.call(rbind, do.call(c, rows))
}


#========================================
# 6. Circular statistics: Watson U2 test
#========================================

# Tests whether transport direction changed between adjacent time bins for each volcano group

# Function to compute Watson U2 statistic and p-value
watson_temporal <- function(data, scenario_name, bearing_col = "bearing") {
  all_bins   <- levels(data$bin_label)
  groups     <- names(VOLCANO_GROUPS)

  rows <- lapply(groups, function(grp) {
    d_grp <- data[data$group == grp, ]

    lapply(seq_along(all_bins[-length(all_bins)]), function(i) {
      bin_young <- all_bins[i]
      bin_old   <- all_bins[i + 1]    # The older bin immediately follows the preceding bin in time, e.g. mid Holocene -> late Holocene

      b_young <- d_grp[[bearing_col]][!is.na(d_grp$bin_label) & d_grp$bin_label == bin_young]    # Extract bearings that fall within bin
      b_old   <- d_grp[[bearing_col]][!is.na(d_grp$bin_label) & d_grp$bin_label == bin_old]
      b_young <- b_young[!is.na(b_young)]
      b_old   <- b_old[!is.na(b_old)]

      n_young <- length(b_young)
      n_old   <- length(b_old)

      U2  <- NA_real_
      sig <- NA_character_

      if (n_young >= 5 && n_old >= 5) {    # Verifies that each bin contains at least 5 bearings
        c1  <- circular(b_young, units = "degrees", template = "geographics")
        c2  <- circular(b_old,   units = "degrees", template = "geographics")
        res <- tryCatch(watson.two.test(c1, c2), error = function(e) NULL)    # Compute Watson U2 statistic; returns NULL if the test fails (e.g. due to identical bearing distribution)
        if (!is.null(res)) {
          U2  <- round(res$statistic, 4)
          sig <- if      (U2 > 0.385) "< 0.001"    # Critical U2 values from print.watson.two.test() in the circular package
                 else if (U2 > 0.268) "< 0.01"     # These are large-value approximations; interpret cautiously when n < 10 in either group
                 else if (U2 > 0.187) "< 0.05"
                 else if (U2 > 0.152) "< 0.10"
                 else                 ">= 0.10"
        }
      }

      data.frame(
        scenario         = scenario_name,
        group            = grp,
        bin_young        = bin_young,
        bin_old          = bin_old,
        n_young          = n_young,
        n_old            = n_old,
        U2               = U2,
        p_approx         = sig,
        insufficient     = n_young < 5 | n_old < 5,
        stringsAsFactors = FALSE
      )
    })
  })
  do.call(rbind, do.call(c, rows))
}


#==========================
# 7. Permutation null test
#==========================

# Test the null hypothesis that the observed Watson U2 statistic could arise by chance
# Pool bearings from adjacent time bins and randomly shuffle N_PERM times (defined above)
# Compute the Watson U2 statistic for each shuffle
# Compare with the observed U2 statistic 
# Permutation p-value = the proportion of shuffled U2 values >= the observed U2 value
# Low perm_p value indicates the observed U2 is unlikely to arise by chance

permutation_test <- function(data, scenario_name, bearing_col = "bearing",
                                      n_perm = N_PERM) {
  set.seed(MC_SEED + 1L)    # Set seed for the random number generator for reproducibility 
  all_bins <- levels(data$bin_label)
  groups   <- names(VOLCANO_GROUPS)

  rows <- lapply(groups, function(grp) {
    d_grp <- data[data$group == grp, ]

    lapply(seq_along(all_bins[-length(all_bins)]), function(i) {
      bin_young <- all_bins[i]
      bin_old   <- all_bins[i + 1]

      b_young <- d_grp[[bearing_col]][!is.na(d_grp$bin_label) & d_grp$bin_label == bin_young]
      b_old   <- d_grp[[bearing_col]][!is.na(d_grp$bin_label) & d_grp$bin_label == bin_old]
      b_young <- b_young[!is.na(b_young)]
      b_old   <- b_old[!is.na(b_old)]

      n_young <- length(b_young)
      n_old   <- length(b_old)

      U2_obs <- NA_real_
      perm_p <- NA_real_

      if (n_young >= 5 && n_old >= 5) {
        c1  <- circular(b_young, units = "degrees", template = "geographics")
        c2  <- circular(b_old,   units = "degrees", template = "geographics")
        res <- tryCatch(watson.two.test(c1, c2), error = function(e) NULL)

        if (!is.null(res)) {
          U2_obs <- round(res$statistic, 4)    # Run the Watson test on observed values first
          b_all  <- c(b_young, b_old)

          U2_perm <- vapply(seq_len(n_perm), function(j) {    # Permutation loop
            perm <- sample(b_all)
            cp1  <- circular(perm[seq_len(n_young)],                    units = "degrees", template = "geographics")    # Original group sizes are preserved
            cp2  <- circular(perm[seq(n_young + 1L, length(perm))],     units = "degrees", template = "geographics")
            r    <- tryCatch(watson.two.test(cp1, cp2), error = function(e) NULL)    # Recompute Watson U2 on shuffled bearings
            if (is.null(r)) NA_real_ else r$statistic
          }, numeric(1))

          n_valid <- sum(!is.na(U2_perm))    # Compute p-value
          perm_p  <- if (n_valid > 0)
            round(mean(U2_perm >= U2_obs, na.rm = TRUE), 4) else NA_real_
        }
      }

      data.frame(    
        scenario     = scenario_name,
        group        = grp,
        bin_young    = bin_young,
        bin_old      = bin_old,
        n_young      = n_young,
        n_old        = n_old,
        U2           = U2_obs,
        perm_p       = perm_p,
        insufficient = n_young < 5 | n_old < 5,
        stringsAsFactors = FALSE
      )
    })
  })
  do.call(rbind, do.call(c, rows))
}


#======================
# 8. Residual bearings
#======================

# Account for geographic bias with residual bearings
# Deposits cluster where sites a) happen to exist and b) have been studied
# Subtract each deposit's source volcano mean bearing from its observed bearing
# Caveat: the bearings of deposits within the bin in question are not included in the 
# source volcano mean bearing to avoid circularity

# Compute circular mean bearing per source volcano
compute_volcano_means <- function(data) {
  valid <- !is.na(data$bearing) & !is.na(data$source_volcano)
  vols  <- split(data$bearing[valid], data$source_volcano[valid])
  vapply(vols, function(b) {
    if (length(b) < MIN_VOLCANO_N) return(NA_real_)    # Minimum of 5 deposits per volcano needed to compute a reliable mean otherwise NA
    r <- b * pi / 180    # Convert to radians
    (atan2(mean(sin(r)), mean(cos(r))) * 180 / pi + 360) %% 360
  }, numeric(1))
}

# Function to compute difference between a deposit's bearing and its volcano's mean
circ_diff <- function(bearing, ref) {
  ((bearing - ref + 180) %% 360) - 180
}

# Convert [0, 360] bearing to signed [−180°, +180°] convention
# 0 = typical direction
# +ive = clockwise shift 
# -ive = anticlockwise shift 
to_signed <- function(deg) {
  ifelse(is.na(deg), NA_real_, ((deg + 180) %% 360) - 180)
}

# Add a new column to the dataset containing residual bearings 
add_bearing_residuals <- function(data) {
  data$bearing_resid <- NA_real_    # Residual column; filled in bin by bin below
  bins <- levels(data$bin_label)

  for (b in bins) {
    deposits_in_bin <- !is.na(data$bin_label) & data$bin_label == b    # Logical index of deposits in this bin
    if (!any(deposits_in_bin)) next    # Skip if bin is empty

    # Leave-one-bin-out: compute volcano means from all bins except the current one
    other_bins    <- data[!is.na(data$bin_label) & data$bin_label != b, ]
    volcano_means <- compute_volcano_means(other_bins)

    # Look up the reference mean for each deposit's source volcano, then subtract it
    ref_mean <- volcano_means[data$source_volcano[deposits_in_bin]]
    data$bearing_resid[deposits_in_bin] <- circ_diff(data$bearing[deposits_in_bin], ref_mean)
  }
  data
}


#=============================================
# 9. Monte Carlo sampling for age uncertainty
#=============================================

# For each deposit, sample an age from a normal distribution centred on the median age,
# truncated to the reported [age_min, age_max] range, which is treated as the 95% (+/-2 sigma) interval
# This is repeated N_MC times (defined above) to generate a distribution of possible ages for each deposit
# Each iteration re-bins the data and re-runs all core analyses

sample_ages <- function(data) {
  data$age_sampled <- data$age_median
  valid <- !is.na(data$age_min) & !is.na(data$age_max) &
    is.finite(data$age_min) & is.finite(data$age_max) &
    data$age_max > data$age_min
  if (any(valid)) {
    mu    <- data$age_median[valid]
    sigma <- (data$age_max[valid] - data$age_min[valid]) / (2 * 1.96)    # age_min/age_max treated as the 95% CI (+/-1.96 sigma)

    # Truncated normal via inverse-CDF sampling, so ages never fall outside the reported range
    p_lo <- pnorm(data$age_min[valid], mean = mu, sd = sigma)
    p_hi <- pnorm(data$age_max[valid], mean = mu, sd = sigma)
    u    <- runif(sum(valid), min = p_lo, max = p_hi)

    data$age_sampled[valid] <- qnorm(u, mean = mu, sd = sigma)
  }
  data
}

circ_mean_deg <- function(x) {
  a <- x[!is.na(x)]
  if (length(a) == 0) return(NA_real_)
  r <- a * pi / 180
  (atan2(mean(sin(r)), mean(cos(r))) * 180 / pi + 360) %% 360
}

# Aggregation functions collapse the N_MC result tables into summaries per bin and group
# (% significant, mean direction, mean U2, etc.)

aggregate_mc_uniformity <- function(mc_list) {
  all_res <- do.call(rbind, lapply(seq_along(mc_list), function(i)
    cbind(mc_list[[i]], mc_iter = i)))
  keys <- unique(all_res[, c("group", "bin")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
    d <- all_res[all_res$group == keys$group[k] & all_res$bin == keys$bin[k], ]
    data.frame(
      group        = keys$group[k],
      bin          = keys$bin[k],
      bin_mid      = mean(d$bin_mid, na.rm = TRUE),
      n_median     = round(median(d$n, na.rm = TRUE)),
      mean_dir_mc  = round(circ_mean_deg(d$mean_dir), 1),
      R_mean       = round(mean(d$R, na.rm = TRUE), 3),
      R_sd         = round(sd(d$R, na.rm = TRUE), 3),
      pct_sig      = round(100 * mean(d$significant, na.rm = TRUE), 1),
      U_mean       = round(mean(d$rao_U, na.rm = TRUE), 2),
      U_sd         = round(sd(d$rao_U, na.rm = TRUE), 2),
      pct_sig_rao  = round(100 * mean(d$significant_rao, na.rm = TRUE), 1),
      stringsAsFactors = FALSE
    )
  }))
}

# Aggregate N_MC watson_temporal() tables into one row per group × transition
aggregate_mc_watson_temp <- function(mc_list) {
  all_res <- do.call(rbind, lapply(seq_along(mc_list), function(i)
    cbind(mc_list[[i]], mc_iter = i)))
  keys <- unique(all_res[, c("group", "bin_old", "bin_young")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
    d      <- all_res[all_res$group     == keys$group[k]     &
                        all_res$bin_old   == keys$bin_old[k]   &
                        all_res$bin_young == keys$bin_young[k], ]
    tested <- !d$insufficient & !is.na(d$U2)
    data.frame(
      group     = keys$group[k],
      bin_old   = keys$bin_old[k],
      bin_young = keys$bin_young[k],
      n_old_med = round(median(d$n_old)),
      n_yng_med = round(median(d$n_young)),
      U2_mean   = round(mean(d$U2[tested], na.rm = TRUE), 4),
      pct_p010  = round(100 * mean(tested & d$U2 > 0.152, na.rm = TRUE), 1),
      pct_p005  = round(100 * mean(tested & d$U2 >= 0.187, na.rm = TRUE), 1),
      pct_p001  = round(100 * mean(tested & d$U2 >= 0.385, na.rm = TRUE), 1),
      stringsAsFactors = FALSE
    )
  }))
}

# MC driver: re-samples ages, re-bins, and re-runs core analyses N_MC times
# Also runs residual bearing analysis on each iteration
run_mc <- function(data, scenario_name, collapse_eruptions = FALSE) {
  set.seed(MC_SEED)
  mc_r <- mc_wt <- mc_r_resid <- mc_wt_resid <- vector("list", N_MC)

  cat(sprintf("MC analysis — %s scenario (%d iterations)...\n", scenario_name, N_MC))
  for (i in seq_len(N_MC)) {
    if (i %% 250 == 0) cat(sprintf("  %d / %d\n", i, N_MC))
    d_mc       <- sample_ages(data)
    d_mc       <- assign_bins(d_mc, age_col = "age_sampled")
    d_mc_resid <- add_bearing_residuals(d_mc)    # residuals always computed on deposit-level data first
    if (collapse_eruptions) {
      d_mc       <- collapse_to_eruptions(d_mc)
      d_mc_resid <- collapse_to_eruptions(d_mc_resid, bearing_col = "bearing_resid")
    }
    mc_r[[i]]        <- summarise_by_group(d_mc,       scenario_name)
    mc_wt[[i]]       <- watson_temporal(d_mc,          scenario_name)
    mc_r_resid[[i]]  <- summarise_by_group(d_mc_resid, scenario_name, bearing_col = "bearing_resid")
    mc_wt_resid[[i]] <- watson_temporal(d_mc_resid,    scenario_name, bearing_col = "bearing_resid")
  }
  cat("  done.\n")

  list(
    uniformity        = aggregate_mc_uniformity(mc_r),
    watson_temp       = aggregate_mc_watson_temp(mc_wt),
    uniformity_resid  = aggregate_mc_uniformity(mc_r_resid),
    watson_temp_resid = aggregate_mc_watson_temp(mc_wt_resid)
  )
}


#=======================
# 10. Eruption collapse
#=======================

# Function to collapse data to one observation per eruption per bin, using the circular mean of all its deposits
# Deposits with no eruption label are kept as individual observations
collapse_to_eruptions <- function(data, bearing_col = "bearing") {
  erup_id <- data$source_eruption
  missing <- is.na(erup_id) | erup_id == ""
  erup_id[missing] <- paste0(".deposit.", seq_len(nrow(data))[missing])    # Unknown eruptions are given a unique placeholder ID
                                                                           # This allows all unknown eruptions to be counted separately
  bin_levels    <- levels(data$bin_label)    # Preserve bin levels so empty bins are not lost
  data$.bearing <- data[[bearing_col]]  

  # Group and summarise (by eruption ID and age bin)
  out <- data |>
    mutate(.erup_id = erup_id) |>
    group_by(.erup_id, source_volcano, group, bin_label, bin_mid) |>    # Grouping by source volcano allows a multi-source deposit to be
    summarise(                                                          # counted twice under max scenario even when eruption is unknown
      .bearing    = circ_mean_deg(.bearing),    # Calcular mean bearing of all deposits attributed to the same eruption in a given bin
      distance_km = mean(distance_km, na.rm = TRUE),    
      .groups     = "drop"
    ) |>
    select(-.erup_id) |>
    data.frame()

  names(out)[names(out) == ".bearing"] <- bearing_col
  out$bin_label <- factor(out$bin_label, levels = bin_levels)

  # Residual bearings use the signed [-180, 180] convention
  # circ_mean_deg() always returns [0, 360) so convert back here
  if (bearing_col == "bearing_resid") {
    out[[bearing_col]] <- to_signed(out[[bearing_col]])
  }

  out
}


#=========================
# -- ANALYSIS PIPELINE ---
#=========================

#=========================
# 11. Load and clean data
#=========================

raw <- read_excel(DATA_FILE, sheet = DATA_SHEET)

# Rename spreadsheet columns
raw <- rename(raw,
  source_volcano        = `Source volcano`,
  source_eruption       = `Source eruption`,
  age_min               = `Minimum age (cal yr BP)`,
  age_median            = `Median age (cal yr BP)`,
  age_max               = `Maximum age (cal yr BP)`,
  distance_km           = `Great circle distance between source and deposit (km)`,
  bearing               = `Angle between source and deposit (bearing in degrees)`,
  distance_km_secondary = `Great circle distance between source and deposit (km) - secondary`,
  bearing_secondary     = `Angle between source and deposit (bearing in degrees) - secondary`
)

# Work on a plain data.frame 
raw <- data.frame(raw, check.names = FALSE, stringsAsFactors = FALSE)

# Ensures key columns are numerics
numeric_cols <- c("age_median", "age_min", "age_max", "bearing", "distance_km",
                  "distance_km_secondary", "bearing_secondary")
raw <- suppressWarnings(mutate(raw, across(any_of(numeric_cols), as.numeric)))

# Fail loudly on an empty or mis-read data sheet
if (nrow(raw) == 0)
  stop("Check DATA_FILE and DATA_SHEET.")

# Account for deposits with multiple populations
# source_volcano may contain two or more comma-separated sources
# Allow a multi-source deposit to contribute more than one observation
# However note that these observations are not independent
expand_sources <- function(data) {
  parts <- strsplit(as.character(data$source_volcano), ",", fixed = TRUE)    # Split cells after comma if present
  n1 <- trimws(vapply(parts, function(p) p[1], character(1)))    # n1 takes the first name from the cell
  n2 <- trimws(vapply(parts,    # n2 takes the second name
                      function(p) if (length(p) >= 2) p[2] else NA_character_,
                      character(1)))
  
  prim <- data    # Create a copy of the full dataset that only contains primary populations
  prim$source_volcano <- n1
  prim$source_rank    <- "primary"
  
  sec_ok <- !is.na(n2) & !is.na(data$distance_km_secondary) & !is.na(data$bearing_secondary)
  
  if (any(sec_ok)) {
    sec <- data[sec_ok, ]    # Create a subset of data only containing secondary population rows
    sec$source_volcano <- n2[sec_ok]
    sec$distance_km    <- sec$distance_km_secondary    # Move secondary source geometry into standard column 
    sec$bearing        <- sec$bearing_secondary
    sec$source_rank    <- "secondary"
    out <- rbind(prim, sec[, names(prim)])
  } else {
    out <- prim
  }
  out
}

df <- expand_sources(raw)    # Call the function


df <- df[!is.na(df$age_median) & !is.na(df$age_min) & !is.na(df$age_max) &
           !is.na(df$bearing) & !is.na(df$distance_km), ]    # Filter out observations lacking essential information fields

df <- df[df$distance_km >= MIN_DISTANCE_KM, ]    # Filter out proximal deposits 

cat(sprintf("Observations after cleaning: %d (distance >= %d km)\n",
            nrow(df), MIN_DISTANCE_KM))
cat(sprintf("  of which secondary-source observations: %d\n",
            sum(df$source_rank == "secondary")))


#============================
# 12. Assign groups and bins
#============================

# Label each deposit with its group
df$group_min <- assign_group(df$source_volcano, "min")
df$group_max <- assign_group(df$source_volcano, "max")

# Minimum: certain attributions and primary sources only
df_min <- df |> filter(!is.na(group_min), source_rank == "primary") |> rename(group = group_min)
# Maximum: uncertain attributions and secondary sources included
df_max <- df |> filter(!is.na(group_max))                            |> rename(group = group_max)

cat(sprintf("\nMinimum scenario: %d deposits\n", nrow(df_min)))
cat(sprintf("Maximum scenario: %d deposits\n\n", nrow(df_max)))

cat("Group counts (minimum scenario):\n")
print(table(df_min$group))

cat("\nGroup counts (maximum scenario):\n")
print(table(df_max$group))

# Bin deposits using the active binning scheme
df_min <- assign_bins(df_min)
df_max <- assign_bins(df_max)

# Collapse to one observation per eruption per bin
df_min_erup <- collapse_to_eruptions(df_min)
df_max_erup <- collapse_to_eruptions(df_max)

cat(sprintf("\nEruption-level (minimum scenario): %d eruptions\n", nrow(df_min_erup)))
cat(sprintf("Eruption-level (maximum scenario): %d eruptions\n", nrow(df_max_erup)))


#===========================
# 13. Run statistical tests
#===========================

# Run Rayleigh and Rao's spacing tests 
# Produce summaries by group and bin
results_min <- summarise_by_group(df_min, "Minimum")
results_max <- summarise_by_group(df_max, "Maximum")

# Run Watson U2 temporal tests 
temporal_min <- watson_temporal(df_min, "Minimum")
temporal_max <- watson_temporal(df_max, "Maximum")

# Run permutation null test for Watson U2 temporal
cat("\nRunning permutation null tests (primary bearings, ~1 min) ...\n")
perm_min <- permutation_test(df_min, "Minimum")
perm_max <- permutation_test(df_max, "Maximum")

# Monte Carlo over age uncertainty
mc_min <- run_mc(df_min, "Minimum")
mc_max <- run_mc(df_max, "Maximum")


# Statistical tests on eruption-level bearings
# -----------------------------------------------
results_erup_min <- summarise_by_group(df_min_erup, "Minimum (eruption)")
results_erup_max <- summarise_by_group(df_max_erup, "Maximum (eruption)")

temporal_erup_min <- watson_temporal(df_min_erup, "Minimum (eruption)")
temporal_erup_max <- watson_temporal(df_max_erup, "Maximum (eruption)")

cat("\nRunning permutation null tests (eruption-level bearings) ...\n")
perm_erup_min <- permutation_test(df_min_erup, "Minimum (eruption)")
perm_erup_max <- permutation_test(df_max_erup, "Maximum (eruption)")

mc_erup_min <- run_mc(df_min, "Minimum (eruption)", collapse_eruptions = TRUE)
mc_erup_max <- run_mc(df_max, "Maximum (eruption)", collapse_eruptions = TRUE)


# Statistical tests on residual bearings
# -----------------------------------------
df_min_resid <- add_bearing_residuals(df_min)
df_max_resid <- add_bearing_residuals(df_max)

# Print n deposits that successfully received a residual bearing
# Not possible to compute the required reference mean for those with fewer than MIN_VOLCANO_N deposits 
cat(sprintf("\nResidual coverage — Minimum scenario: %d / %d deposits\n",    
            sum(!is.na(df_min_resid$bearing_resid)), nrow(df_min_resid)))    
cat(sprintf("Residual coverage — Maximum scenario: %d / %d deposits\n",       
            sum(!is.na(df_max_resid$bearing_resid)), nrow(df_max_resid)))    

# Run Rayleigh and Rao's spacing tests on residual bearings
results_resid_min  <- summarise_by_group(df_min_resid, "Minimum (residual)", bearing_col = "bearing_resid")
results_resid_max  <- summarise_by_group(df_max_resid, "Maximum (residual)", bearing_col = "bearing_resid")

# Fix sign conventions (-180,180)
results_resid_min$mean_dir <- to_signed(results_resid_min$mean_dir)
results_resid_max$mean_dir <- to_signed(results_resid_max$mean_dir)

# Run Watson U2 test on the residual bearings
temporal_resid_min <- watson_temporal(df_min_resid, "Minimum (residual)", bearing_col = "bearing_resid")
temporal_resid_max <- watson_temporal(df_max_resid, "Maximum (residual)", bearing_col = "bearing_resid")

# Run permutation null test for residual bearings
cat("\nRunning permutation null tests (residual bearings) ...\n")
perm_resid_min <- permutation_test(df_min_resid, "Minimum (residual)", bearing_col = "bearing_resid")
perm_resid_max <- permutation_test(df_max_resid, "Maximum (residual)", bearing_col = "bearing_resid")


# Statistical tests on eruption-level residual bearings
# --------------------------------------------------------
df_min_resid_erup <- collapse_to_eruptions(df_min_resid, bearing_col = "bearing_resid")
df_max_resid_erup <- collapse_to_eruptions(df_max_resid, bearing_col = "bearing_resid")

results_resid_erup_min <- summarise_by_group(df_min_resid_erup, "Minimum (eruption, residual)", bearing_col = "bearing_resid")
results_resid_erup_max <- summarise_by_group(df_max_resid_erup, "Maximum (eruption, residual)", bearing_col = "bearing_resid")
results_resid_erup_min$mean_dir <- to_signed(results_resid_erup_min$mean_dir)
results_resid_erup_max$mean_dir <- to_signed(results_resid_erup_max$mean_dir)

temporal_resid_erup_min <- watson_temporal(df_min_resid_erup, "Minimum (eruption, residual)", bearing_col = "bearing_resid")
temporal_resid_erup_max <- watson_temporal(df_max_resid_erup, "Maximum (eruption, residual)", bearing_col = "bearing_resid")

cat("\nRunning permutation null tests (eruption-level residual bearings) ...\n")
perm_resid_erup_min <- permutation_test(df_min_resid_erup, "Minimum (eruption, residual)", bearing_col = "bearing_resid")
perm_resid_erup_max <- permutation_test(df_max_resid_erup, "Maximum (eruption, residual)", bearing_col = "bearing_resid")


#===========================
# 14. Write results to CSVs
#===========================

combine_mc_uniformity <- function(mc_raw, mc_resid, mc_erup_resid) {
  rename_block <- function(d, suffix) {
    stat_cols <- setdiff(names(d), c("group", "bin", "bin_mid"))
    names(d)[names(d) %in% stat_cols] <- paste0(stat_cols, "_", suffix)
    d
  }
  raw   <- rename_block(mc_raw,        "raw")
  resid <- rename_block(mc_resid,      "resid")
  erup  <- rename_block(mc_erup_resid, "erup_resid")

  d <- merge(raw,   resid, by = c("group", "bin", "bin_mid"), all = TRUE)
  d <- merge(d,     erup,  by = c("group", "bin", "bin_mid"), all = TRUE)

  d$group <- factor(d$group, levels = names(VOLCANO_GROUPS))
  d$bin   <- factor(d$bin,   levels = BIN_LABELS)
  d[order(d$group, d$bin), ]
}

mc_uniformity_combined_min <- combine_mc_uniformity(mc_min$uniformity, mc_min$uniformity_resid,
                                                     mc_erup_min$uniformity_resid)
mc_uniformity_combined_max <- combine_mc_uniformity(mc_max$uniformity, mc_max$uniformity_resid,
                                                     mc_erup_max$uniformity_resid)
                                                     
write.csv(mc_uniformity_combined_min, "mc_uniformity_combined_min.csv", row.names = FALSE)
write.csv(mc_uniformity_combined_max, "mc_uniformity_combined_max.csv", row.names = FALSE)
write.csv(perm_min, "perm_watson_temp_min.csv", row.names = FALSE)
write.csv(perm_max, "perm_watson_temp_max.csv", row.names = FALSE)
write.csv(perm_resid_min, "perm_watson_resid_min.csv", row.names = FALSE)
write.csv(perm_resid_max, "perm_watson_resid_max.csv", row.names = FALSE)
write.csv(mc_min$watson_temp, "mc_watson_temp_min.csv", row.names = FALSE)
write.csv(mc_max$watson_temp, "mc_watson_temp_max.csv", row.names = FALSE)
write.csv(mc_min$watson_temp_resid, "mc_watson_temp_resid_min.csv", row.names = FALSE)
write.csv(mc_max$watson_temp_resid, "mc_watson_temp_resid_max.csv", row.names = FALSE)
write.csv(results_resid_min, "residual_uniformity_min.csv", row.names = FALSE)
write.csv(results_resid_max, "residual_uniformity_max.csv", row.names = FALSE)
write.csv(temporal_resid_min,"residual_watson_temp_min.csv",  row.names = FALSE)
write.csv(temporal_resid_max,"residual_watson_temp_max.csv",  row.names = FALSE)
write.csv(results_erup_min,         "eruption_uniformity_min.csv",    row.names = FALSE)
write.csv(results_erup_max,         "eruption_uniformity_max.csv",    row.names = FALSE)
write.csv(temporal_erup_min,        "eruption_watson_temp_min.csv",   row.names = FALSE)
write.csv(temporal_erup_max,        "eruption_watson_temp_max.csv",   row.names = FALSE)
write.csv(perm_erup_min,            "perm_watson_erup_min.csv",       row.names = FALSE)
write.csv(perm_erup_max,            "perm_watson_erup_max.csv",       row.names = FALSE)
write.csv(mc_erup_min$watson_temp,       "mc_watson_temp_erup_min.csv",          row.names = FALSE)
write.csv(mc_erup_max$watson_temp,       "mc_watson_temp_erup_max.csv",          row.names = FALSE)
write.csv(results_resid_erup_min,        "eruption_residual_uniformity_min.csv",  row.names = FALSE)
write.csv(results_resid_erup_max,        "eruption_residual_uniformity_max.csv",  row.names = FALSE)
write.csv(temporal_resid_erup_min,       "eruption_residual_watson_temp_min.csv", row.names = FALSE)
write.csv(temporal_resid_erup_max,       "eruption_residual_watson_temp_max.csv", row.names = FALSE)
write.csv(perm_resid_erup_min,           "perm_watson_erup_resid_min.csv",        row.names = FALSE)
write.csv(perm_resid_erup_max,           "perm_watson_erup_resid_max.csv",        row.names = FALSE)
write.csv(mc_erup_min$watson_temp_resid, "mc_watson_temp_erup_resid_min.csv",     row.names = FALSE)
write.csv(mc_erup_max$watson_temp_resid, "mc_watson_temp_erup_resid_max.csv",     row.names = FALSE)
cat("CSVs written.\n")


#===============================
# 15. Save results for plotting
#===============================

analysis_results <- list(
  grouping_scheme    = GROUPING_SCHEME,
  binning_scheme     = BINNING_SCHEME,
  bin_breaks         = BIN_BREAKS,
  bin_labels         = BIN_LABELS,
  volcano_groups     = VOLCANO_GROUPS,
  min_n              = MIN_N,
  n_mc               = N_MC,
  n_perm             = N_PERM,
  df_min             = df_min,
  df_max             = df_max,
  df_min_resid       = df_min_resid,
  df_max_resid       = df_max_resid,
  results_min        = results_min,
  results_max        = results_max,
  results_resid_min  = results_resid_min,
  results_resid_max  = results_resid_max,
  temporal_min       = temporal_min,
  temporal_max       = temporal_max,
  temporal_resid_min = temporal_resid_min,
  temporal_resid_max = temporal_resid_max,
  mc_min             = mc_min,
  mc_max             = mc_max,
  mc_min_resid       = mc_min$uniformity_resid,
  mc_max_resid       = mc_max$uniformity_resid,
  mc_wt_min_resid    = mc_min$watson_temp_resid,
  mc_wt_max_resid    = mc_max$watson_temp_resid,
  perm_min           = perm_min,
  perm_max           = perm_max,
  perm_resid_min     = perm_resid_min,
  perm_resid_max     = perm_resid_max,
  df_min_erup        = df_min_erup,
  df_max_erup        = df_max_erup,
  results_erup_min   = results_erup_min,
  results_erup_max   = results_erup_max,
  temporal_erup_min  = temporal_erup_min,
  temporal_erup_max  = temporal_erup_max,
  perm_erup_min      = perm_erup_min,
  perm_erup_max      = perm_erup_max,
  mc_erup_min          = mc_erup_min,
  mc_erup_max          = mc_erup_max,
  df_min_resid_erup    = df_min_resid_erup,
  df_max_resid_erup    = df_max_resid_erup,
  results_resid_erup_min  = results_resid_erup_min,
  results_resid_erup_max  = results_resid_erup_max,
  temporal_resid_erup_min = temporal_resid_erup_min,
  temporal_resid_erup_max = temporal_resid_erup_max,
  perm_resid_erup_min  = perm_resid_erup_min,
  perm_resid_erup_max  = perm_resid_erup_max
)

RESULTS_FILE <- sprintf("analysis_results_%s_%s.rds", GROUPING_SCHEME, BINNING_SCHEME)
saveRDS(analysis_results, RESULTS_FILE)
cat("\nAnalysis bundle written to", RESULTS_FILE, "\n")

# Plotting occurs in Dispersal_analysis_figures_final.R
