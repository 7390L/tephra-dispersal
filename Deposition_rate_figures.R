
#=================================
# ~~~ Deposition rate figures ~~~
#=================================

# Plotting only; analysis occurs in Monteath_adapted.R and SPD.R

# Contents:
# ------------
#   1. Packages
#   2. Configuration
#   3. Set global style theme
#   4. Monteath et al. (2025): Main figure
#        Plot cumulative deposition and all three competing models
#        ~ regime-shift, total-sites, site-type
#   5. SPD and CPL model fits figure
#   6. Simple SPD figure
#   7. Synthesis figure: Monteath-style vs SPD analysis
#   8. Save figures


# =================
# 1. Load packages
# =================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scico)
library(scales)
library(showtext)
library(sysfonts)
library(systemfonts)
library(patchwork)


#==================
# 2. Configuration
#==================

# Choose which cached results to use for plotting 

#  "italian" = Italian-sourced tephra only, "all" = every tephra in the dataset 

ANALYSIS_SCOPE <- "italian"    # Change me: "italian", "all"

# Cached results from the Monteath et al. (2025)-style analysis alway contain both the min and max scenarios 
MONTEATH_ANALYSIS_CACHE <- paste0("Monteath_analysis_cache_", ANALYSIS_SCOPE, ".RData")

# Cached results from the SPD analysis save min and max scenarios separately
SOURCE_SCENARIO <- "max"    # Change me: "min", "max"
                            # This toggle also selects which of the min/max scenarios from the
                            # Monteath analysis is plotted

# Phasing scheme the cached SPD_final.R run used:
# "deposit" = each dated deposit is its own phase
# "eruption" = deposits from the same named eruption are grouped into one phase
PHASE_TYPE <- "deposit"    # Change me: "deposit", "eruption"

SPD_ANALYSIS_CACHE <- paste0("SPD_analysis_cache_", ANALYSIS_SCOPE, "_", SOURCE_SCENARIO,
                             "_by_", PHASE_TYPE, ".RData")

# Label for figure titles
# Add adjective to figure title when scope is restricted to "italian"
scope_adj <- if (ANALYSIS_SCOPE == "italian") "Italian " else ""

# Show/hide the 95% CI error bars on each breakpoint's estimated location
SHOW_BREAKPOINT_CI <- TRUE

# Show/hide the site-availability-adjusted CPL fits
SHOW_SITEAVAIL_CPL <- TRUE


#===========================
# 3. Set global style theme
#===========================

font_add_google("Source Sans 3", regular.wt = 400, bold.wt = 500)
showtext_auto()
showtext_opts(dpi = 400)

thesis_theme <- theme_classic(base_size = 14, base_family = "Source Sans 3") +
  theme(
    plot.title         = element_text(size = 14, face = "bold", hjust = 0.5,
                                       margin = margin(b = 20)),   
    axis.text          = element_text(size = 11.5, colour = "grey20"),
    axis.text.x        = element_text(size = 10.5, colour = "grey20"),
    axis.text.y        = element_text(size = 10.5, colour = "grey20"),
    axis.title.x       = element_text(size = 11.5, colour = "grey20",margin = margin(t = 18)),
    axis.title.y       = element_text(size = 11.5, colour = "grey20", margin = margin(r = 18)),
    legend.title       = element_text(size = 10.5, face = "bold"),
    legend.text        = element_text(size = 10.5),
    legend.key.size    = unit(0.8, "lines"),   
    legend.key.spacing.y = unit(5, "pt"),      
    legend.margin      = margin(3, 5, 3, 5),   
    strip.text         = element_text(size = 13, colour = "grey20"),
    panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.3),
    axis.line          = element_line(colour = "grey20"),
    axis.ticks         = element_line(colour = "grey20")
  )

n_label_layer <- function(label, order) {
  list(
    geom_point(data = data.frame(x = NA_real_, y = NA_real_, n_grp = label),
               aes(x = x, y = y, shape = n_grp), alpha = 0, na.rm = TRUE, inherit.aes = FALSE),
    scale_shape_manual(name = NULL, values = setNames(NA, label)),
    guides(shape = guide_legend(order = order, override.aes = list(alpha = 0)))
  )
}


#=======================================
# 4. Monteath et al (2025): Main figure
#=======================================

# Builds one figure overlaying both cumulative deposition scenarios (min/max)
# and all three competing models (regime-shift, site-type, total-sites)

# Colour palette
col_min    <- "#314E5E"    # Minimum regime-shift line
col_max    <- "#863F37"   # Maximum regime-shift line
col_min_fill <- "#B7C1C7"  # Observed minimum ribbon fill 
col_max_fill <- "#D5BCB9"  # Observed maximum ribbon fill  
col_site_type     <- "#bd9145"   # Site-type model (line)
col_site_type_fill <- "#fce2b5"    # Site-type model CI ribbon fill
col_total_sites     <- "#4D4D4D"      # Total-sites model (line)
col_total_sites_fill <- "#d5d6d0"   # Total-sites model CI ribbon fill 

# Load cached analysis results
if (!file.exists(MONTEATH_ANALYSIS_CACHE)) {
  stop("No cache found at '", MONTEATH_ANALYSIS_CACHE, "'. Run Monteath_adapted_final.R ",
       "first to generate it.")
}
load(MONTEATH_ANALYSIS_CACHE)
message("Loaded cached analysis results from ", MONTEATH_ANALYSIS_CACHE)

# Number of deposits included in each attribution scenario, shown as a caption on the figures
n_min <- nrow(tephra_data_min)
n_max <- nrow(tephra_data_max)

# Observed cumulative deposition
obs_wide <- bins %>%
  arrange(desc(bin_center)) %>%
  transmute(bin_center,
            min_cum = cumulative_intensity_corrected_min_tephras,
            max_cum = cumulative_intensity_corrected_max_tephras)

obs_steps <- obs_wide %>%
  mutate(x_old = bin_center, x_new = bin_center - 500) %>%   # 500 = bin width
  pivot_longer(c(x_old, x_new), values_to = "x") %>%
  dplyr::select(-name) %>%
  arrange(desc(x))

# Regime-shift fit lines
fit_long <- bind_rows(pred_df_min, pred_df_max) %>%
  arrange(type, bin_center) %>%
  distinct(type, bin_center, .keep_all = TRUE)

# Competing-model cumulative series (site-type + total-sites) including their 95% CIs
build_model_df <- function(bins, scenario) {
  p <- function(stub) paste0("cumulative_predicted_", scenario, "_tephras", stub)
  b <- bins %>% arrange(desc(bin_center))
  bind_rows(
    data.frame(bin_center = b$bin_center, y = b[[p("")]],
               lower = b[[p("_lower")]], upper = b[[p("_upper")]],
               series = "Site-type"),
    data.frame(bin_center = b$bin_center, y = b[[p("_T")]],
               lower = b[[p("_lower_T")]], upper = b[[p("_upper_T")]],
               series = "Total-sites")
  )
}

# Site-type and total-sites models:
# Toggle SOURCE_SCENARIO above to switch which scenario's models are plotted
alt_scenario_label <- if (SOURCE_SCENARIO == "max") "Maximum" else "Minimum"
site_type_label      <- paste0("Site-type (", SOURCE_SCENARIO, ")")
total_sites_label    <- paste0("Total-sites (", SOURCE_SCENARIO, ")")
site_type_ci_label   <- "Site-type 95% CI"
total_sites_ci_label <- "Total-sites 95% CI"

alt_models <- build_model_df(bins, SOURCE_SCENARIO) %>%
  mutate(scenario = alt_scenario_label,
         colour_key = ifelse(series == "Site-type", site_type_label, total_sites_label))

# Regime-shift lines are coloured to match their scenario's observed ribbon
# Keeps both min and max scenarios
regime_lines <- fit_long %>%
  transmute(bin_center, y = pmax(pred.fit, 0),   
            series     = "Regime-shift",  
            scenario   = ifelse(type == "Min Tephras", "Minimum", "Maximum"),
            colour_key = ifelse(type == "Min Tephras",
                                "Regime-shift (min)",
                                "Regime-shift (max)"))

model_lines <- bind_rows(alt_models, regime_lines) %>%
  mutate(series   = factor(series, levels = c("Regime-shift",
                                              "Site-type",
                                              "Total-sites")),
         scenario = factor(scenario, levels = c("Minimum", "Maximum")))

breakpoints_min <- dplyr::filter(all_psi, type == "Min Tephras") %>%
  transmute(bin_center = Est., y = pred.fit, colour_key = "Regime-shift (min)",
            xmin = Est. - 1.96 * `St.Err`, xmax = Est. + 1.96 * `St.Err`)
breakpoints_max <- dplyr::filter(all_psi, type == "Max Tephras") %>%
  transmute(bin_center = Est., y = pred.fit, colour_key = "Regime-shift (max)",
            xmin = Est. - 1.96 * `St.Err`, xmax = Est. + 1.96 * `St.Err`)
breakpoints <- bind_rows(breakpoints_min, breakpoints_max)

# Breakpoint error bars
breakpoint_ci_layer <- if (SHOW_BREAKPOINT_CI) {
  geom_segment(data = breakpoints,
               aes(x = xmin, xend = xmax, y = y, yend = y, colour = colour_key),
               linewidth = 0.5, show.legend = FALSE)
} else {
  NULL
}

# Major climatic periods shown as bands behind the data
climate_periods <- data.frame(
  period = c("Last Glacial Maximum", "Heinrich Stadial 1", "Bølling-Allerød",
             "Younger Dryas", "Early Holocene", "Mid-Holocene", "Late Holocene"),
  start  = c(25000, 19000, 14700, 12900, 11700, 8200, 4200),
  end    = c(19000, 14700, 12900, 11700,  8200, 4200,    0)
) %>%
  mutate(mid   = (start + end) / 2,
         shade = rep(c("#fcf9f9", "#f5f1f1"), length.out = n()))

# Position period labels just above the highest line/ribbon in the plot
plot_ymax <- max(c(obs_wide$max_cum, model_lines$y, alt_models$upper), na.rm = TRUE)
climate_periods$label_y <- plot_ymax * 1.05

# Plot:
# ------
combined_models_plot <- ggplot() +
  geom_rect(data = climate_periods,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = I(shade)),
            colour = NA) +
  geom_ribbon(data = obs_steps,
              aes(x = x, ymin = 0, ymax = min_cum, fill = "Cumulative (min)"),
              alpha = 1, colour = NA, key_glyph = draw_key_rect) +
  geom_ribbon(data = obs_steps,
              aes(x = x, ymin = min_cum, ymax = max_cum, fill = "Cumulative (max)"),
              alpha = 1, colour = NA, key_glyph = draw_key_rect) +
  scale_fill_manual(name = "Observed tephra deposition",
                    values = c("Cumulative (min)" = col_min_fill,
                               "Cumulative (max)" = col_max_fill)) +
  geom_ribbon(data = dplyr::filter(alt_models, series == "Total-sites"),
              aes(x = bin_center, ymin = lower, ymax = upper,
                  colour = total_sites_ci_label),
              fill = col_total_sites_fill, alpha = 1,
              key_glyph = draw_key_rect) +
  geom_ribbon(data = dplyr::filter(alt_models, series == "Site-type"),
              aes(x = bin_center, ymin = lower, ymax = upper,
                  colour = site_type_ci_label),
              fill = col_site_type_fill, alpha = 0.7, linewidth = 0,
              key_glyph = draw_key_rect) +
  geom_line(data = model_lines,
            aes(bin_center, y, colour = colour_key, linetype = colour_key,
                alpha = scenario, group = interaction(series, scenario)),
            linewidth = 0.9) +
  geom_segment(data = breakpoints,
               aes(x = bin_center, xend = bin_center, y = 0, yend = y,
                   colour = colour_key),
               linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  breakpoint_ci_layer +
  geom_point(data = breakpoints,
             aes(x = bin_center, y = y, colour = colour_key),
             size = 2) +
  geom_text(data = climate_periods,
            aes(x = mid, y = label_y, label = period),
            inherit.aes = FALSE, angle = 90, hjust = 1, vjust = 0.5,
            size = 3.5, colour = "grey30") +
  scale_x_reverse(name = "Age (cal yrs BP)", breaks = seq(0, 25000, 5000),
                  minor_breaks = seq(0, 25000, 1000),
                  guide = guide_axis(minor.ticks = TRUE),
                  limits = c(25000, 0), expand = c(0, 0)) +
  scale_y_continuous(name = "Site-averaged cumulative tephra deposition",
                     minor_breaks = scales::breaks_width(1),
                     guide = guide_axis(minor.ticks = TRUE),
                     expand = expansion(mult = c(0, 0), add = c(0, 0.1))) +
  scale_colour_manual(name = "Competing models",
                      breaks = c("Regime-shift (max)",
                                 "Regime-shift (min)",
                                 total_sites_label,
                                 total_sites_ci_label,
                                 site_type_label,
                                 site_type_ci_label),
                      values = setNames(
                        c(col_max, col_min, col_site_type, col_total_sites,
                          col_site_type_fill, col_total_sites_fill),
                        c("Regime-shift (max)", "Regime-shift (min)",
                          site_type_label, total_sites_label,
                          site_type_ci_label, total_sites_ci_label))) +
  scale_linetype_manual(values = setNames(
                          c("solid", "solid", "solid", "dotted"),
                          c("Regime-shift (max)", "Regime-shift (min)",
                            site_type_label, total_sites_label))) +
  scale_alpha_manual(values = c("Minimum" = 0.9, "Maximum" = 0.9)) +
  guides(fill     = guide_legend(order = 1),
         colour   = guide_legend(order = 2,
                                 override.aes = list(
                                   linetype = c("solid", "solid", "dotted", "blank",
                                                "solid", "blank"),
                                   alpha    = c(0.9, 0.9, 1, 1, 1, 1))),
         linetype = "none",
         alpha    = "none") +
  labs(title = paste0("Monteath et al. (2025) analysis: Cumulative ", scope_adj,
                      "tephra deposition in the Mediterranean from LGM to present")) +
  n_label_layer(paste0("n = ", n_min, " (min), ", n_max, " (max)"), order = 3) +
  thesis_theme +
  theme(plot.subtitle          = element_text(size = 10, colour = "grey30", margin = margin(b = 12)),
        axis.line              = element_line(linewidth = 0.6),
        axis.ticks             = element_line(linewidth = 0.4),
        axis.minor.ticks.length = unit(0.08, "cm"),
        axis.ticks.length      = unit(0.16, "cm"),
        panel.grid             = element_blank(),
        panel.grid.major.y     = element_blank(),
        legend.position        = "right")

print(combined_models_plot)


#==================================
# 5. SPD and CPL model fits figure
#==================================

# Load cached SPD analysis results
if (!file.exists(SPD_ANALYSIS_CACHE)) {
  stop("No cache found at '", SPD_ANALYSIS_CACHE, "'. Run SPD_final.R first to generate it.")
}
load(SPD_ANALYSIS_CACHE)
message("Loaded cached SPD analysis results from ", SPD_ANALYSIS_CACHE)

# Number of deposits included in this attribution scenario, shown as a caption
n_tephra <- nrow(tephra)

# Colours for the CPL model orders 
pal <- c("1-CPL" = "#65724A", "2-CPL" = "#29295C", "3-CPL" = "#eb819f",
         "4-CPL" = "#C1811A", "5-CPL" = "#8E5572")

# Exclude models that didn't converge reliably across the 5 restarts
converged_models <- cpl_bic_table$model[cpl_bic_table$neglogL_range <= CPL_CONVERGENCE_TOL]
all_cpl_df_plot <- dplyr::filter(all_cpl_df, model %in% converged_models)

# Plot higher order CPL models above lower order models
model_levels <- paste0(sort(as.integer(sub("-CPL", "", unique(all_cpl_df_plot$model)))), "-CPL")
all_cpl_df_plot$model <- factor(all_cpl_df_plot$model, levels = model_levels)

# CPL fits overlaid on the SPD
# SPD rescaled to a probability density (area = 1) so it shares the y-axis with the CPL fits
spd_df_rate <- data.frame(age = years, density = spd_vals)
spd_df_rate$density <- spd_df_rate$density / (sum(spd_df_rate$density) * grid_inc)

# MCMC 95% credible interval for the BIC-selected reference model
best_label    <- paste0(sub("cpl", "", BEST_CPL_KEY), "-CPL")
best_ci_label <- paste0(best_label, " 95% credible interval")
mcmc_ci_best  <- mcmc_results[[BEST_CPL_KEY]]$ci

model_order      <- as.integer(sub("-CPL", "", model_levels))
alpha_values     <- setNames(ifelse(model_order <= 2, 0.8, 1), model_levels)   # Give lower-order models lower opacity
linetype_values  <- setNames(rep("solid", length(model_levels)), model_levels)
if ("3-CPL" %in% model_levels && ANALYSIS_SCOPE == "italian" && SOURCE_SCENARIO == "max") {
  linetype_values[["3-CPL"]] <- "dotted"
}

cpl_plot <- ggplot() +
  geom_area(data = spd_df_rate, aes(age, density), fill = "grey85") +
  geom_ribbon(data = mcmc_ci_best,
              aes(x = age, ymin = lower, ymax = upper,
                  fill = best_ci_label),
              inherit.aes = FALSE, alpha = 0.15) +
  geom_line(data = all_cpl_df_plot, aes(age, density, colour = model, alpha = model,
                                        linetype = model),
            linewidth = 0.8) +
  scale_colour_manual(values = pal, breaks = model_levels) +
  scale_alpha_manual(values = alpha_values, guide = "none") +
  scale_linetype_manual(values = linetype_values, guide = "none") +
  scale_fill_manual(name = NULL,
                    values = setNames(pal[[best_label]], best_ci_label)) +
  guides(colour = guide_legend(order = 1), fill = guide_legend(order = 2)) +
  scale_x_reverse(name = "Age (cal yrs BP)",
                  breaks = seq(0, 25000, 5000),
                  minor_breaks = seq(0, 25000, 1000),
                  guide = guide_axis(minor.ticks = TRUE),
                  expand = c(0, 0)) +
  scale_y_continuous(name = "Probability density",
                     labels = scales::label_number(),
                     breaks = seq(0, 0.0004, 0.0001),
                     minor_breaks = seq(0, 0.0004, 0.00005),
                     guide = guide_axis(minor.ticks = TRUE)) +
  # Crop the tallest SPD spikes out of view so the CPL fits are readable
  coord_cartesian(ylim = c(0, 0.0004), expand = FALSE) +
  thesis_theme +
  theme(panel.grid.major.y = element_blank(),
        axis.line = element_line(colour = "grey20"),
        axis.minor.ticks.length = unit(0.08, "cm"),
        axis.ticks.length      = unit(0.16, "cm"),
        legend.spacing.y      = unit(0, "pt")) +
  labs(title = paste0("Summed probability distribution of ", scope_adj,
                      "tephra deposits with continuous\npiecewise linear model fits from LGM to present (",
                      alt_scenario_label, " scenario)"),
       colour = "Model") +
  n_label_layer(paste0("n = ", n_tephra), order = 3)

print(cpl_plot)


#  Site-availability-corrected CPL fits figure 
#----------------------------------------------

# Same as above, but uses site-availability adjusted CPL model fits

cpl_siteavail_plot <- NULL

if (SHOW_SITEAVAIL_CPL && HAS_SITEAVAIL) {
  
  all_cpl_ts_df$model <- sub("-CPL\\+ts", "-CPL + TS", all_cpl_ts_df$model)

  best_ts_label    <- paste0(sub("cpl", "", BEST_SITEAVAIL_KEY), "-CPL + TS")
  best_ts_ci_label <- paste0(best_ts_label, " 95% credible interval")

  model_levels_ts <- paste0(sort(as.integer(sub("-CPL \\+ TS", "", unique(all_cpl_ts_df$model)))), "-CPL + TS")
  all_cpl_ts_df$model <- factor(all_cpl_ts_df$model, levels = model_levels_ts)

  pal_ts <- setNames(pal, sub("-CPL$", "-CPL + TS", names(pal)))

  model_order_ts     <- as.integer(sub("-CPL \\+ TS", "", model_levels_ts))
  alpha_values_ts    <- setNames(ifelse(model_order_ts <= 2, 0.8, 1), model_levels_ts)
  linetype_values_ts <- setNames(rep("solid", length(model_levels_ts)), model_levels_ts)

  mcmc_ci_ts_layer <- if (!is.null(mcmc_ci_siteavail)) {
    geom_ribbon(data = mcmc_ci_siteavail, aes(x = age, ymin = lower, ymax = upper, fill = best_ts_ci_label),
                inherit.aes = FALSE, alpha = 0.15)
  } else {
    NULL
  }

  cpl_siteavail_plot <- ggplot() +
    geom_area(data = spd_df_rate, aes(age, density), fill = "grey85") +
    mcmc_ci_ts_layer +
    geom_line(data = all_cpl_ts_df, aes(age, density, colour = model, alpha = model,
                                        linetype = model),
              linewidth = 0.8) +
    scale_colour_manual(values = pal_ts, breaks = model_levels_ts) +
    scale_alpha_manual(values = alpha_values_ts, guide = "none") +
    scale_linetype_manual(values = linetype_values_ts, guide = "none") +
    scale_fill_manual(name = NULL, values = setNames(pal_ts[[best_ts_label]], best_ts_ci_label)) +
    guides(colour = guide_legend(order = 1), fill = guide_legend(order = 2)) +
    scale_x_reverse(name = "Age (cal yrs BP)",
                    breaks = seq(0, 25000, 5000),
                    minor_breaks = seq(0, 25000, 1000),
                    guide = guide_axis(minor.ticks = TRUE),
                    expand = c(0, 0)) +
    scale_y_continuous(name = "Probability density",
                       labels = scales::label_number(),
                       breaks = seq(0, 0.0004, 0.0001),
                       minor_breaks = seq(0, 0.0004, 0.00005),
                       guide = guide_axis(minor.ticks = TRUE)) +

    coord_cartesian(ylim = c(0, 0.0004), expand = FALSE) +
    thesis_theme +
    theme(panel.grid.major.y = element_blank(),
          axis.line = element_line(colour = "grey20"),
          axis.minor.ticks.length = unit(0.08, "cm"),
          axis.ticks.length      = unit(0.16, "cm"),
          legend.spacing.y      = unit(0, "pt")) +
    labs(title = paste0("Site-availability-corrected (CPL+timeseries) model fits for ", scope_adj,
                        "tephra deposits\nfrom LGM to present (", alt_scenario_label, " scenario)"),
         colour = "Model") +
    n_label_layer(paste0("n = ", n_tephra), order = 3)

  print(cpl_siteavail_plot)
} else if (SHOW_SITEAVAIL_CPL) {
  message("Skipping site-availability-corrected figure: SPD_final.R's cache has no site-",
          "availability results for this scope/scenario/phase combination")
}


cpl_combined_plot <- if (!is.null(cpl_siteavail_plot)) {
  cpl_plot + cpl_siteavail_plot + plot_layout(ncol = 2, guides = "keep")
} else {
  NULL
}


#======================
# 6. Simple SPD figure
#======================

# Simple SPD, uncropped and without CPL fits overlaid 
raw_spd_plot <- ggplot(spd_df_rate, aes(age, density)) +
  geom_area(fill = "grey85") +
  scale_x_reverse(name = "Age (cal yrs BP)",
                  breaks = seq(0, 25000, 5000),
                  minor_breaks = seq(0, 25000, 1000),
                  guide = guide_axis(minor.ticks = TRUE),
                  expand = c(0, 0)) +
  scale_y_continuous(name = "Probability density",
                     labels = scales::label_number(),
                     breaks = seq(0, 0.0012, 0.0002),
                     minor_breaks = seq(0, 0.0012, 0.0001),
                     guide = guide_axis(minor.ticks = TRUE),
                     expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 0.0012), expand = FALSE) +
  thesis_theme +
  theme(panel.grid.major.y = element_blank(),
        axis.line = element_line(colour = "grey20"),
        axis.minor.ticks.length = unit(0.08, "cm"),
        axis.ticks.length      = unit(0.16, "cm")) +
  labs(title = paste0("Summed probability distribution of ", scope_adj,
                      "tephra deposits\n from LGM to present (",
                      alt_scenario_label, " scenario)")) +
  n_label_layer(paste0("n = ", n_tephra), order = 1)

print(raw_spd_plot)


#=====================================================
# 7. Synthesis figure: Monteath-style vs SPD analysis
#=====================================================

psi_type    <- if (SOURCE_SCENARIO == "max") "Max Tephras" else "Min Tephras"
obs_cum_col <- if (SOURCE_SCENARIO == "max") "max_cum" else "min_cum"

col_max_bin  <- if (SOURCE_SCENARIO == "max") col_max else col_min   
col_bin_fill <- if (SOURCE_SCENARIO == "max") col_max_fill else col_min_fill   
col_spd_fill <- "#E4D0B2"  

if (HAS_SITEAVAIL) {
  synth_pieces <- sub("cpl", "", BEST_SITEAVAIL_KEY)
  synth_label  <- paste0(synth_pieces, "-CPL + TS")
} else {
  synth_label  <- best_label
}
col_cpl_hinge <- pal[[best_label]]   

obs_max   <- obs_steps %>% mutate(y = .data[[obs_cum_col]])
bps_bin   <- dplyr::filter(all_psi, type == psi_type) %>%
  mutate(xmin = Est. - 1.96 * `St.Err`, xmax = Est. + 1.96 * `St.Err`)
spd_curve <- spd_df_ord %>%
  dplyr::select(age, cum_int = cum_spd_intensity) %>% filter(!is.na(cum_int))
bestcpl_curve <- (if (HAS_SITEAVAIL) bestcplts_cum_df else bestcpl_cum_df) %>%
  dplyr::select(age, cum_int)
bestcpl_hinge_df <- if (HAS_SITEAVAIL) bestcplts_hinge_df else bestcpl_hinge_df
bin_fit   <- dplyr::filter(fit_long, type == psi_type) %>%
  dplyr::select(age = bin_center, fit = pred.fit)

regime_bin_label <- paste0("Regime-shift (", SOURCE_SCENARIO, ")")
bestcpl_label    <- if (HAS_SITEAVAIL) synth_label else paste0("Intensity-corrected ", synth_label, " model")

bps_bin_ci_layer <- if (SHOW_BREAKPOINT_CI) {
  geom_segment(data = bps_bin,
               aes(x = xmin, xend = xmax, y = pred.fit, yend = pred.fit,
                   colour = regime_bin_label),
               linewidth = 0.5, show.legend = FALSE)
} else {
  NULL
}

# BIC-selected model's hinge (breakpoint) markers, with 95% credible intervals from the MCMC chain 
bestcpl_hinge_ci_layer <- if (SHOW_BREAKPOINT_CI) {
  geom_segment(data = bestcpl_hinge_df,
               aes(x = xmin, xend = xmax, y = cum_int, yend = cum_int),
               colour = col_cpl_hinge, linewidth = 0.5)
} else {
  NULL
}

# Major climatic periods shown as bands behind the data
climate_periods <- data.frame(
  period = c("Last Glacial Maximum", "Heinrich Stadial 1", "Bølling-Allerød",
             "Younger Dryas", "Early Holocene", "Mid-Holocene", "Late Holocene"),
  start  = c(25000, 19000, 14700, 12900, 11700, 8200, 4200),
  end    = c(19000, 14700, 12900, 11700,  8200, 4200,    0)
) %>%
  mutate(mid   = (start + end) / 2,
         shade = rep(c("#fcf9f9", "#f5f1f1"), length.out = n()))

plot_ymax <- max(c(obs_max$y, spd_curve$cum_int, bin_fit$fit, bestcpl_curve$cum_int), na.rm = TRUE)
climate_periods$label_y <- plot_ymax * 1.05

# Plot
combined_plot <- ggplot() +
  geom_rect(data = climate_periods,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = I(shade)),
            colour = NA) +
  geom_ribbon(data = obs_max,
              aes(x = x, ymin = 0, ymax = y, fill = "Median-age binning"),
              alpha = 1, colour = NA) +
  geom_area(data = spd_curve,
            aes(x = age, y = cum_int, fill = "Intensity-corrected SPD"),
            alpha = 0.5, colour = "#C1811A", linewidth = 0.3) +
  geom_line(data = bin_fit, aes(x = age, y = fit, colour = regime_bin_label),
            linewidth = 1.0, alpha = 0.9) +
  geom_line(data = bestcpl_curve, aes(x = age, y = cum_int, colour = bestcpl_label),
            linewidth = 1.0, alpha = 0.9) +
  geom_segment(data = bps_bin,
               aes(x = Est., xend = Est., y = 0, yend = pred.fit, colour = regime_bin_label),
               linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  bps_bin_ci_layer +
  geom_point(data = bps_bin, aes(Est., pred.fit, colour = regime_bin_label), size = 2) +
  geom_segment(data = bestcpl_hinge_df,
               aes(x = age, xend = age, y = 0, yend = cum_int),
               colour = col_cpl_hinge, linetype = "dashed", linewidth = 0.5) +
  bestcpl_hinge_ci_layer +
  geom_point(data = bestcpl_hinge_df, aes(x = age, y = cum_int, colour = bestcpl_label), size = 2) +
  geom_text(data = climate_periods,
            aes(x = mid, y = label_y, label = period),
            inherit.aes = FALSE, angle = 90, hjust = 1, vjust = 0.5,
            size = 3.5, colour = "grey30") +
  scale_x_reverse(name = "Age (cal yrs BP)", breaks = seq(0, 25000, 5000),
                  minor_breaks = seq(0, 25000, 1000),
                  guide = guide_axis(minor.ticks = TRUE),
                  limits = c(25000, 0), expand = c(0, 0)) +
  scale_y_continuous(name = "Site-averaged cumulative tephra deposition",
                     minor_breaks = scales::breaks_width(1),
                     guide = guide_axis(minor.ticks = TRUE),
                     limits = c(0, NA), expand = expansion(mult = c(0, 0), add = c(0, 0.1))) +
  scale_fill_manual(name = "Observed tephra deposition",
                    breaks = c("Intensity-corrected SPD", "Median-age binning"),
                    values = c("Median-age binning" = col_bin_fill,
                               "Intensity-corrected SPD"  = col_spd_fill)) +
  scale_colour_manual(name = "Model fits",
                      breaks = c(bestcpl_label, regime_bin_label),
                      values = setNames(c(col_cpl_hinge, col_max_bin), c(bestcpl_label, regime_bin_label))) +
  guides(fill   = guide_legend(order = 1, override.aes = list(alpha = c(0.5, 1))),
         colour = guide_legend(order = 2, override.aes = list(linewidth = 1.0, alpha = c(0.9, 0.9)))) +
  labs(title = paste0("Cumulative ", scope_adj,
                      "tephra deposition: Intensity-corrected SPD and ", synth_label, " model\n compared with ",
                      "Monteath et al. (2025) regime-shift model (",
                      alt_scenario_label, " scenario)")) +
  n_label_layer(paste0("n = ", n_tephra), order = 3) +
  thesis_theme +
  theme(panel.grid.major.y   = element_blank(),
        axis.line            = element_line(linewidth = 0.6),
        axis.ticks           = element_line(linewidth = 0.4),
        axis.minor.ticks.length = unit(0.08, "cm"),
        axis.ticks.length      = unit(0.16, "cm"),
        legend.position      = "right")

print(combined_plot)


#=================
# 8. Save figures
#=================

# Save figures to the same folder
OUTPUT_DIR <- "Deposition_rate_figures"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

ggsave(file.path(OUTPUT_DIR, paste0("combined_models_plot_", ANALYSIS_SCOPE, "_", SOURCE_SCENARIO, ".png")),
       combined_models_plot, width = 13, height = 8, dpi = 400)

spd_fig_suffix <- paste0("_", ANALYSIS_SCOPE, "_", SOURCE_SCENARIO, "_by_", PHASE_TYPE)

ggsave(file.path(OUTPUT_DIR, paste0("cpl_plot", spd_fig_suffix, ".png")),
       cpl_plot, width = 9.75, height = 6, dpi = 400)

if (!is.null(cpl_siteavail_plot)) {
  ggsave(file.path(OUTPUT_DIR, paste0("cpl_siteavail_plot", spd_fig_suffix, ".png")),
         cpl_siteavail_plot, width = 9.75, height = 6, dpi = 400)
}

if (!is.null(cpl_combined_plot)) {
  ggsave(file.path(OUTPUT_DIR, paste0("cpl_plot_and_siteavail_combined", spd_fig_suffix, ".png")),
         cpl_combined_plot, width = 19, height = 6, dpi = 400)
}

ggsave(file.path(OUTPUT_DIR, paste0("raw_spd_plot", spd_fig_suffix, ".png")),
       raw_spd_plot, width = 9.75, height = 6, dpi = 400)

ggsave(file.path(OUTPUT_DIR, paste0("combined_plot", spd_fig_suffix, ".png")),
       combined_plot, width = 13, height = 8, dpi = 400)
