
# ======================================================
# ~~~ Tephra Transport Direction Analysis (figures) ~~~
# ======================================================

# Contents:
# ------------
#   1. Packages
#   2. Configuration
#   3. Load analysis results
#   4. Set global style theme
#   5. Rose diagrams
#   6. Composite Watson U2 and permutation test figure


# =================
# 1. Load packages
# =================

library(dplyr)
library(ggplot2)
library(scico)
library(showtext)
library(systemfonts)
library(sysfonts)
library(patchwork)
library(cowplot) 


#==================
# 2. Configuration
#==================

# Grouping/binning schemes
# Toggle these to select which cached analysis_results bundle to load
# Generate a scheme's cache by setting these in Dispersal_analysis.R 

GROUPING_SCHEME <- "refined"   # Change me: "refined", "broad", "italian", "mediterranean"
BINNING_SCHEME  <- "deglac_fine"   # Change me: "climate", "deglac_fine", "fixed_width"

RESULTS_FILE <- sprintf("analysis_results_%s_%s.rds", GROUPING_SCHEME, BINNING_SCHEME)

# Every figure is built for both the min and max attribution scenarios, and
# for all three bearing-type/collapse combinations, in one run


# =========================
# 3. Load analysis results
# =========================

if (!file.exists(RESULTS_FILE)) {
  stop("No cache found at '", RESULTS_FILE, "'. Run Dispersal_analysis_final.R with ",
       "GROUPING_SCHEME = \"", GROUPING_SCHEME, "\" and BINNING_SCHEME = \"", BINNING_SCHEME,
       "\" first to generate it.")
}
analysis_results <- readRDS(RESULTS_FILE)
list2env(analysis_results, envir = environment())    # Copy every element of analysis_results into current environment
message("Loaded cached analysis results from ", RESULTS_FILE)

stopifnot(grouping_scheme == GROUPING_SCHEME, binning_scheme == BINNING_SCHEME)

select_df <- function(scenario, bearing_type, collapse) {
  suffix <- paste0(if (bearing_type == "residual") "_resid" else "",
                    if (collapse) "_erup" else "")
  get(paste0("df_", scenario, suffix))
}

select_results <- function(scenario, bearing_type, collapse) {
  suffix <- paste0(if (bearing_type == "residual") "_resid" else "",
                    if (collapse) "_erup" else "")
  get(paste0("results", suffix, "_", scenario))
}

select_mc_uniformity <- function(scenario, bearing_type, collapse) {
  base <- if (collapse) get(paste0("mc_erup_", scenario)) else get(paste0("mc_", scenario))
  if (bearing_type == "residual") base$uniformity_resid else base$uniformity
}

to_signed <- function(deg) {
  ifelse(is.na(deg), NA_real_, ((deg + 180) %% 360) - 180)
}

bearing_col <- function(bearing_type) if (bearing_type == "residual") "bearing_resid" else "bearing"

# Bearing axis differs by bearing_type: observed bearings are compass
# directions (0-360 deg); residuals are signed deviations from each source
# volcano's mean bearing (-180 to +180, 0 = baseline direction)
rose_axis <- function(bearing_type) {
  if (bearing_type == "residual") {
    list(limits = c(-180, 180),
         breaks = seq(-180, 135, by = 45),
         labels = c("±180°", "-135°", "-90°", "-45°",
                    "0° (baseline)", "+45°", "+90°", "+135°"),
         start  = pi)
  } else {
    list(limits = c(0, 360),
         breaks = seq(0, 315, by = 45),
         labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
         start  = 0)
  }
}

SCHEME_TAG <- paste(GROUPING_SCHEME, BINNING_SCHEME, sep = "_")


# ====================
# 4. Set global theme
# ====================

font_add_google("Source Sans 3", regular.wt = 400, bold.wt = 600)
showtext_auto()
showtext_opts(dpi = 400)

thesis_theme <- theme_classic(base_size = 14, base_family = "Source Sans 3") +
  theme(
    plot.title         = element_text(size = 16, face = "bold", hjust = 0.5,
                                      margin = margin(b = 6)),   # bold title
    plot.subtitle      = element_text(size = 12, colour = "grey40", hjust = 0.5,
                                      margin = margin(b = 16)),
    axis.text          = element_text(size = 14, colour = "grey20"),
    axis.text.x        = element_text(size = 14, colour = "grey20"),
    axis.text.y        = element_text(size = 14, colour = "grey20"),
    axis.title.x       = element_text(margin = margin(t = 18)),
    axis.title.y       = element_text(margin = margin(r = 18)),
    legend.title       = element_text(size = 14),
    legend.text        = element_text(size = 13),
    strip.text         = element_text(size = 13, colour = "grey20"),
    panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.3),
    axis.line          = element_line(colour = "grey20")
  )

FIG_DIR <- "Dispersal_analysis_figures"
dir.create(FIG_DIR, showWarnings = FALSE)

save_fig <- function(p, filename, width = 14, height = 8) {
  if (!is.null(p)) ggsave(file.path(FIG_DIR, filename), p, width = width, height = height, dpi = 400)
}

# Volcano group colours
# ------------------------

.sc <- scico(12, palette = "romaO")   # Colours drawn evenly from scico 'romaO' palette

GROUP_COLOUR_SCHEMES <- list(
  refined       = c(Campanian = .sc[7], Etna = .sc[11], Aeolian = .sc[3]),
  broad         = c(Campanian = .sc[2], "Non-Campanian" = .sc[8]),
  italian       = c("All Italian" = .sc[5]),
  mediterranean = c("All Mediterranean" = .sc[5])
)

group_colours <- GROUP_COLOUR_SCHEMES[[grouping_scheme]]


# =================
# 5. Rose diagrams
# =================

BOTTOM_BEARING <- 180


ROSE_ROWS <- list(
  list(bearing_type = "observed", collapse = FALSE),
  list(bearing_type = "residual", collapse = FALSE),
  list(bearing_type = "residual", collapse = TRUE)
)

rose_row_max_count <- function(group_name, scenario, bearing_type, collapse) {
  bcol <- bearing_col(bearing_type)
  data <- select_df(scenario, bearing_type, collapse)
  data <- data[data$group == group_name & !is.na(data$bin_label) & !is.na(data[[bcol]]), ]
  if (nrow(data) == 0) return(0)
  data$bin_label <- factor(data$bin_label, levels = bin_labels)
  p <- ggplot(data, aes(x = .data[[bcol]])) +
    geom_histogram(binwidth = 20, boundary = 0) +
    facet_grid(. ~ bin_label)
  max(ggplot_build(p)$data[[1]]$count, na.rm = TRUE)
}

build_rose_row <- function(group_name, scenario, bearing_type, collapse, r_scale, show_bin_labels = TRUE) {
  bcol <- bearing_col(bearing_type)
  axis <- rose_axis(bearing_type)
  data <- select_df(scenario, bearing_type, collapse)
  data <- data[data$group == group_name & !is.na(data$bin_label) & !is.na(data[[bcol]]), ]
  data$bin_label <- factor(data$bin_label, levels = bin_labels)

  p <- ggplot(data, aes(x = .data[[bcol]])) +
    geom_histogram(
      binwidth  = 20,
      boundary  = 0,
      colour    = "white",
      linewidth = 0.1,
      fill      = group_colours[[group_name]]
    ) +
    coord_polar(start = axis$start, direction = 1, clip = "off") +
    scale_x_continuous(limits = axis$limits, breaks = axis$breaks, labels = axis$labels) +
    facet_grid(. ~ bin_label) +
    labs(x = NULL, y = NULL) +
    thesis_theme +
    theme(
      axis.text.x       = element_text(size = 9, colour = "#676767"),
      axis.ticks.y      = element_blank(),
      axis.text.y       = element_blank(),
      panel.grid.major  = element_line(colour = "#cececeff", linewidth = 0.3),
      panel.spacing.x   = unit(1.2, "cm"),
      strip.background  = element_blank(),
      strip.text.x      = if (show_bin_labels) element_text(size = 13, colour = "grey20") else element_blank(),
      axis.line         = element_blank()
    )

  # Point-estimate statistics
  pt_stats <- select_results(scenario, bearing_type, collapse)
  pt_stats <- pt_stats[pt_stats$group == group_name, ]
  pt_stats <- data.frame(
    bin_label = factor(pt_stats$bin, levels = bin_labels),
    n         = pt_stats$n,
    mean_dir  = pt_stats$mean_dir,
    R         = pt_stats$R
  )
  pt_stats <- pt_stats[!is.na(pt_stats$bin_label), ]

  # MC-mean statistics - averaged across the 1000 age-resampled iterations
  mc_stats <- select_mc_uniformity(scenario, bearing_type, collapse)
  mc_stats <- mc_stats[mc_stats$group == group_name, ]
  mc_stats <- data.frame(
    bin_label   = factor(mc_stats$bin, levels = bin_labels),
    n_mc        = mc_stats$n_median,
    mean_dir_mc = if (bearing_type == "residual") to_signed(mc_stats$mean_dir_mc) else mc_stats$mean_dir_mc,
    R_mc        = mc_stats$R_mean,
    R_mc_sd     = mc_stats$R_sd,
    pct_sig_mc  = mc_stats$pct_sig
  )
  mc_stats <- mc_stats[!is.na(mc_stats$bin_label), ]

  pt_arrows <- pt_stats[!is.na(pt_stats$mean_dir) & !is.na(pt_stats$R), ]
  mc_arrows <- mc_stats[!is.na(mc_stats$mean_dir_mc) & !is.na(mc_stats$R_mc), ]
  arrow_data <- rbind(
    data.frame(bin_label = pt_arrows$bin_label,
               mean_dir = pt_arrows$mean_dir, r_len = pt_arrows$R * r_scale,
               source = rep("Point estimate", nrow(pt_arrows))),
    data.frame(bin_label = mc_arrows$bin_label,
               mean_dir = mc_arrows$mean_dir_mc, r_len = mc_arrows$R_mc * r_scale,
               source = rep("MC mean", nrow(mc_arrows)))
  )
  arrow_data$source <- factor(arrow_data$source, levels = c("Point estimate", "MC mean"))

  # Ring radius
  n_rings <- 4
  ring_breaks <- r_scale * (seq_len(n_rings - 1) / n_rings) * (0.45 / 0.4)
  p <- p + scale_y_continuous(limits = c(0, r_scale), breaks = ring_breaks,
                               expand = c(0, 0), oob = scales::oob_keep)

  # Labels below roses
  fmt_deg <- function(x) ifelse(is.na(x), "—", sprintf("%.0f°", x))
  fmt_r   <- function(x) ifelse(is.na(x), "—", sprintf("%.2f", x))
  fmt_int <- function(x) ifelse(is.na(x), "—", as.character(round(x)))

  label_data <- merge(pt_stats, mc_stats, by = "bin_label", all = TRUE)
  label_data$label <- ifelse(
    !is.na(label_data$n) & label_data$n == 0,
    "n = 0",
    sprintf(
      "n = %s, μ = %s, R = %s\nMC: n = %s, μ = %s, R = %s ± %s",
      fmt_int(label_data$n), fmt_deg(label_data$mean_dir), fmt_r(label_data$R),
      fmt_int(label_data$n_mc), fmt_deg(label_data$mean_dir_mc),
      fmt_r(label_data$R_mc), fmt_r(label_data$R_mc_sd)
    )
  )
  label_data$y <- r_scale * 1.35

  p +
    geom_segment(
      data          = arrow_data,
      mapping       = aes(x = mean_dir, xend = mean_dir, y = 0, yend = r_len,
                           colour = source, linetype = source, linewidth = source),
      inherit.aes   = FALSE,
      lineend       = "round",
      arrow         = arrow(length = unit(0.11, "cm"), type = "closed")
    ) +
    scale_colour_manual(values = c("Point estimate" = "grey10", "MC mean" = "#C97C99"), name = NULL) +
    scale_linetype_manual(values = c("Point estimate" = "solid", "MC mean" = "42"), name = NULL) +
    scale_linewidth_manual(values = c("Point estimate" = 0.7, "MC mean" = 0.45), name = NULL) +
    geom_text(
      data        = label_data,
      mapping     = aes(x = BOTTOM_BEARING, y = y, label = label),
      inherit.aes = FALSE,
      size        = 9,
      size.unit   = "pt",
      colour      = "black",
      lineheight  = 0.9,
      vjust       = 1
    ) +
    theme(legend.position    = "bottom",
          legend.background  = element_blank(),
          legend.key         = element_blank(),
          legend.box.spacing = unit(0.4, "cm"),
          strip.text         = element_text(size = 10, colour = "grey20"))
}

# Build the full 3-row composite for one volcano group x scenario
plot_rose_group <- function(group_name, scenario, scenario_name) {
  r_scale <- max(vapply(ROSE_ROWS, function(r)
    rose_row_max_count(group_name, scenario, r$bearing_type, r$collapse), numeric(1)))

  if (!is.finite(r_scale) || r_scale <= 0) {
    message("No data to plot for ", group_name, " (", scenario_name, ").")
    return(invisible(NULL))
  }

  rows <- Map(function(r, i)
    build_rose_row(group_name, scenario, r$bearing_type, r$collapse, r_scale,
                    show_bin_labels = (i == 1)),
    ROSE_ROWS, seq_along(ROSE_ROWS))

  subtitle_main <- sprintf("%s scenario · %s / %s", scenario_name, GROUPING_SCHEME, BINNING_SCHEME)
  subtitle_note <- paste(strwrap(
    paste("Rows top to bottom: observed, residual, residual eruption-collapsed bearings.",
          "For residual rows, +ve = clockwise deviation, -ve = anticlockwise."),
    width = 110), collapse = "\n")

  (rows[[1]] / rows[[2]] / rows[[3]]) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title    = group_name,
      subtitle = paste0(subtitle_main, "\n", subtitle_note),
      theme    = thesis_theme
    ) &
    theme(legend.position = "bottom")
}

for (group_name in names(volcano_groups)) {
  p_rose_min <- plot_rose_group(group_name, "min", "Minimum")
  p_rose_max <- plot_rose_group(group_name, "max", "Maximum")
  print(p_rose_min)
  print(p_rose_max)
  save_fig(p_rose_min, sprintf("rose_%s_%s_min.png", SCHEME_TAG, tolower(group_name)),
           width = 11.69, height = 9.5)
  save_fig(p_rose_max, sprintf("rose_%s_%s_max.png", SCHEME_TAG, tolower(group_name)),
           width = 11.69, height = 9.5)
}


# ===================================================
# 6. Composite Watson U2 and permutation test figure
# ===================================================

# Produces 6 subplots showing Watson U2 and permutation test results
# Subplots are saved individually and assembled externally

# Includes: 
# 1. Observed bearings, deposit-level
# 2. Residual bearings, deposit-level
# 3. Residual bearings, eruption-collapsed

# Watson U2 critical-value thresholds, shown as reference lines/labels on
# the lollipop panels below
CRIT <- data.frame(
  U2  = c(0.152, 0.187, 0.268, 0.385),
  lab = c("italic(p) < 0.10", "italic(p) < 0.05",
          "italic(p) < 0.01", "italic(p) < 0.001")
)

PANEL_BUNDLES <- list(
  refined = readRDS("analysis_results_refined_climate.rds"),   # Load cached results
  italian = readRDS("analysis_results_italian_climate.rds")
)

make_scheme_rows <- function(bundle) {
  list(
    list(bundle = bundle, bearing_type = "observed", collapse = FALSE,
         subheading = "Observed bearings"),
    list(bundle = bundle, bearing_type = "residual", collapse = FALSE,
         subheading = "Residual bearings"),
    list(bundle = bundle, bearing_type = "residual", collapse = TRUE,
         subheading = "Residual bearings\nEruption-collapsed")
  )
}

PANEL_SETS <- list(
  refined = make_scheme_rows("refined"),
  italian = make_scheme_rows("italian")
)


GLOBAL_SHAPES <- c(Campanian = 21, Etna = 22, Aeolian = 23, "All Italian" = 24)   # Set one shape per volcano group

LEGEND_BARHEIGHT <- unit(3, "cm")   # Ensure legend colour bars for MC-significance and permutation p-values 
                                    # are the same height

CRIT_LABEL_MARGIN_PT <- 130   # Give the top row bigger margin so critical p value labels have enough space
 
Y_GAP_BOTTOM <- 1.0   # y-axis spacing
Y_GAP_TOP    <- 0

# Ensure that no valid Watson U2 results are dropped due to insufficient space
# in the plot 
WATSON_DODGE_SLACK <- 0.25
Y_EXPAND_BOTTOM <- Y_GAP_BOTTOM - WATSON_DODGE_SLACK
Y_EXPAND_TOP    <- Y_GAP_TOP

save_square_panel <- function(p, filename, panel_size = 3.2) {
  if (is.null(p)) return(invisible(NULL))
  g <- ggplotGrob(p)
  panel_rows <- unique(g$layout$t[g$layout$name == "panel"])
  panel_cols <- unique(g$layout$l[g$layout$name == "panel"])

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  fixed_w <- sum(grid::convertWidth(g$widths[-panel_cols],   "in", valueOnly = TRUE))
  fixed_h <- sum(grid::convertHeight(g$heights[-panel_rows], "in", valueOnly = TRUE))

  ggsave(file.path(FIG_DIR, filename), p,
        width = fixed_w + panel_size, height = fixed_h + panel_size,
        units = "in", dpi = 400)
}

# File naming conventions
panel_select_bundle <- function(bundle, prefix, scenario, bearing_type, collapse) {
  suffix <- paste0(if (bearing_type == "residual") "_resid" else "",
                    if (collapse) "_erup" else "",
                    "_", scenario)
  bundle[[paste0(prefix, suffix)]]
}
panel_select_mc <- function(bundle, scenario, bearing_type, collapse) {
  base <- if (collapse) bundle[[paste0("mc_erup_", scenario)]] else bundle[[paste0("mc_", scenario)]]
  list(uniformity  = if (bearing_type == "residual") base$uniformity_resid  else base$uniformity,
       watson_temp = if (bearing_type == "residual") base$watson_temp_resid else base$watson_temp)
}

# Build the Watson U2 and MC-significance data frame for one panel row
panel_build_watson <- function(bundle, scenario, bearing_type, collapse) {
  temp <- panel_select_bundle(bundle, "temporal", scenario, bearing_type, collapse)
  mc   <- panel_select_mc(bundle, scenario, bearing_type, collapse)$watson_temp
  keys <- c("group", "bin_young", "bin_old")
  d <- merge(temp[, c(keys, "U2", "p_approx", "insufficient")],
             mc[,   c(keys, "pct_p005")], by = keys, all.x = TRUE)

# Ensure all phases and volcano groups are shown on axes even when data is unsufficient
# to compute results
  strip_age  <- function(x) sub("\\s*\\(.*$", "", x)
  phases     <- strip_age(bundle$bin_labels)             
  all_trans  <- paste0(phases[-length(phases)], " → ", phases[-1]) 

  d$transition <- paste0(strip_age(d$bin_old), " → ", strip_age(d$bin_young))
  d$transition <- factor(d$transition, levels = rev(all_trans))
  d$group <- factor(d$group, levels = names(bundle$volcano_groups))
  d
}

# Watson U2 lollipop
plot_watson_panel <- function(d, xmax, show_crit_labels, subheading, show_legend = FALSE) {
  tested <- d[!is.na(d$U2), ]
  tested$ti <- as.numeric(tested$transition)
  tested <- tested |>
    dplyr::group_by(transition) |>
    dplyr::mutate(ypos = ti + (dplyr::row_number(group) - (dplyr::n() + 1) / 2) * 0.20) |>
    dplyr::ungroup()

  shp <- GLOBAL_SHAPES[levels(d$group)]

  p <- ggplot(tested) +
    geom_vline(data = CRIT, aes(xintercept = U2),
               linetype = "dashed", colour = "grey78", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = U2, y = ypos, yend = ypos),
                 colour = "grey55", linewidth = 0.5) +
    geom_point(aes(x = U2, y = ypos, fill = pct_p005, shape = group),
               size = 3.8, stroke = 0.3, colour = "grey30") +
    geom_text(aes(x = U2, y = ypos, label = sprintf("%.3f", U2)),
              hjust = 0, nudge_x = 0.016, size = 16, size.unit = "pt", colour = "grey20")

  if (show_crit_labels) {
    p <- p + geom_text(data = CRIT, aes(x = U2, y = Inf, label = lab),
                        parse = TRUE, angle = 90, vjust = 0.5, hjust = -0.35,
                        size = 14, size.unit = "pt", colour = "grey55")
  }

  p +
    scale_y_continuous(breaks = seq_along(levels(d$transition)),
                       labels = levels(d$transition),
                       limits = c(1 - WATSON_DODGE_SLACK,
                                 length(levels(d$transition))),
                       expand = expansion(add = c(Y_EXPAND_BOTTOM, Y_EXPAND_TOP))) +
    scale_x_continuous(limits = c(0, xmax), breaks = seq(0, 0.5, 0.1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_shape_manual(values = shp, name = "Volcano group", drop = FALSE) +
    scico::scale_fill_scico(palette = "lajolla", direction = -1, limits = c(0, 100),
                            name = "MC iterations\nsignificant (%)") +
    guides(
      shape = guide_legend(order = 1, override.aes = list(fill = "grey70")),
      fill  = guide_colourbar(order = 2, barheight = LEGEND_BARHEIGHT, reverse = TRUE)
    ) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL, title = subheading) +
    thesis_theme +
    theme(
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      plot.title    = element_text(size = 12.5, face = "plain", colour = "grey20",
                                    hjust = 0.5, lineheight = 1.15,
                                    margin = margin(b = 10 + if (show_crit_labels) CRIT_LABEL_MARGIN_PT else 0)),
      plot.subtitle = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      plot.margin   = margin(t = 6, r = 50, b = 6, l = 8)
    )
}

# Permutation heatmap
# Every group defined in the grouping scheme gets its own column even when 
# a group has too little data to test 

plot_permutation_panel <- function(bundle, scenario, bearing_type, collapse,
                                    subheading, show_legend = FALSE) {
  perm_results <- panel_select_bundle(bundle, "perm", scenario, bearing_type, collapse)
  d <- perm_results

  strip_age <- function(x) sub("\\s*\\(.*$", "", x)
  phases    <- strip_age(bundle$bin_labels)
  all_trans <- paste0(phases[-length(phases)], " → ", phases[-1])

  d$transition <- paste0(strip_age(d$bin_old), " → ", strip_age(d$bin_young))
  d$transition <- factor(d$transition, levels = rev(all_trans))
  d$group      <- factor(d$group, levels = names(bundle$volcano_groups))

  is_insufficient <- d$insufficient | is.na(d$U2)
  
  p_label <- ifelse(
    is.na(d$perm_p),          "p = —",
    ifelse(d$perm_p < 0.001,  "p < 0.001",
                              sprintf("p = %.3f", d$perm_p))
  )

  d$tile_label <- ifelse(
    is_insufficient,
    sprintf("n = %d, %d", d$n_old, d$n_young),
    sprintf("%s\n(n = %d, %d)", p_label, d$n_old, d$n_young)
  )

  ggplot(d, aes(x = group, y = transition, fill = perm_p)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = tile_label),
              colour = ifelse(is_insufficient, "grey45", "grey10"),
              size = 16, size.unit = "pt", lineheight = 1.1) +
    scale_y_discrete(expand = expansion(add = c(Y_GAP_BOTTOM, Y_GAP_TOP))) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(colours = colorRampPalette(c("#155A94", "#EFF6FC"))(5),     # Use non-linear colour scale to help distinguish different p-values 
                        values  = seq(0, 1, length.out = 5)^3,                           # i.e. using a linear scale, p = 0.1 and p = 0.01 would appear very similar
                        limits = c(0, 1), na.value = "grey85",                           # despite being very different in terms of significance
                        name = "Permutation\np-value",
                        guide = guide_colourbar(barheight = LEGEND_BARHEIGHT)) +
    labs(
      x     = NULL,
      y     = NULL,
      title = subheading
    ) +
    thesis_theme +
    theme(
      panel.grid          = element_blank(),
      panel.grid.major.y  = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      plot.title      = element_text(size = 12.5, face = "plain", colour = "grey20",
                                      hjust = 0.5, lineheight = 1.15, margin = margin(b = 10)),
      plot.subtitle   = element_blank(),
      axis.line       = element_blank(),      # No axis line or tick marks
      axis.ticks      = element_blank(),
      axis.ticks.length = unit(0, "pt"),
      axis.text.x     = element_text(margin = margin(t = 4)),
      axis.text.y     = element_blank(),
      plot.margin     = margin(t = 6, r = 10, b = 6, l = 8)
    )
}

# Build and save one full set of panels 
build_panel_set <- function(set_name, rows, scenario) {
  watson_data <- lapply(rows, function(r) {
    panel_build_watson(PANEL_BUNDLES[[r$bundle]], scenario, r$bearing_type, r$collapse)
  })
  xmax <- max(vapply(watson_data, function(d) max(d$U2, na.rm = TRUE), numeric(1)),
              CRIT$U2, na.rm = TRUE) * 1.08

  for (i in seq_along(rows)) {
    r <- rows[[i]]
    is_top <- (i == 1)
    row_tag <- sprintf("%s_row%d_%s_%s_%s", set_name, i, scenario, r$bearing_type,
                       if (r$collapse) "eruption" else "deposit")

    p_watson <- plot_watson_panel(watson_data[[i]], xmax,
                                   show_crit_labels = is_top,
                                   subheading       = r$subheading)
    save_square_panel(p_watson, sprintf("panel_%s_watson.png", row_tag))

    p_perm <- plot_permutation_panel(PANEL_BUNDLES[[r$bundle]], scenario, r$bearing_type, r$collapse,
                                      subheading = r$subheading)
    save_square_panel(p_perm, sprintf("panel_%s_perm.png", row_tag))
  }
}

for (set_name in names(PANEL_SETS)) {
  build_panel_set(set_name, PANEL_SETS[[set_name]], "max")   # Build set of panels for max scenario
  build_panel_set(set_name, PANEL_SETS[[set_name]], "min")   # Min scenario
}

# Legends
legend_dummy_watson <- data.frame(
  x = 0.1, y = 1, pct = 50,
  group = factor(names(GLOBAL_SHAPES), levels = names(GLOBAL_SHAPES))
)
legend_ref_watson <- ggplot(legend_dummy_watson, aes(x = x, y = y, shape = group, fill = pct)) +
  geom_point(size = 3.8, stroke = 0.3, colour = "grey30") +
  scale_shape_manual(values = GLOBAL_SHAPES, name = "Volcano group", drop = FALSE) +
  scico::scale_fill_scico(palette = "lajolla", direction = -1, limits = c(0, 100),
                          name = "MC iterations\nsignificant (%)") +
  guides(
    shape = guide_legend(order = 1, override.aes = list(fill = "grey70")),
    fill  = guide_colourbar(order = 2, barheight = LEGEND_BARHEIGHT, reverse = TRUE)
  ) +
  thesis_theme + theme(legend.position = "right")
save_fig(cowplot::ggdraw(cowplot::get_legend(legend_ref_watson)),
         "legend_watson.png", width = 3, height = 4)

legend_ref_perm <- plot_permutation_panel(PANEL_BUNDLES[["refined"]], "max", "observed", FALSE,
                                          subheading = NULL, show_legend = TRUE)
save_fig(cowplot::ggdraw(cowplot::get_legend(legend_ref_perm)),
         "legend_perm.png", width = 2.6, height = 2.9)

PANEL_BUNDLES$refined_deglac_fine <- readRDS("analysis_results_refined_deglac_fine.rds")
PANEL_SETS$deglac_fine <- make_scheme_rows("refined_deglac_fine")

build_panel_set("deglac_fine", PANEL_SETS$deglac_fine, "max")
build_panel_set("deglac_fine", PANEL_SETS$deglac_fine, "min")

# Composite panels for the refine key transitions tests
PANEL_BUNDLES$refined_holocene_onset <- readRDS("analysis_results_refined_holocene_onset.rds")

KEY_BOUNDARIES <- c(
  "~14.5 ka (saddle collapse)"    = "refined_deglac_fine",
  "~11.7 ka (Holocene onset)"     = "refined_holocene_onset"
)

strip_age <- function(x) sub("\\s*\\(.*$", "", x)

panel_build_watson_multi <- function(scenario, bearing_type, collapse) {
  parts <- lapply(names(KEY_BOUNDARIES), function(boundary_label) {
    bundle <- PANEL_BUNDLES[[KEY_BOUNDARIES[[boundary_label]]]]
    temp <- panel_select_bundle(bundle, "temporal", scenario, bearing_type, collapse)
    mc   <- panel_select_mc(bundle, scenario, bearing_type, collapse)$watson_temp
    keys <- c("group", "bin_young", "bin_old")
    d <- merge(temp[, c(keys, "U2", "p_approx", "insufficient")],
               mc[,   c(keys, "pct_p005")], by = keys, all.x = TRUE)

    phases    <- strip_age(bundle$bin_labels)
    all_trans <- paste0(phases[-length(phases)], " → ", phases[-1])
    d$transition <- paste0(strip_age(d$bin_old), " → ", strip_age(d$bin_young))
    d$transition <- factor(d$transition, levels = rev(all_trans))
    d$group    <- factor(d$group, levels = names(bundle$volcano_groups))
    d$boundary <- boundary_label
    list(d = d, levels = all_trans)
  })

  combined <- do.call(rbind, lapply(parts, function(x) {
    x$d$transition <- as.character(x$d$transition)
    x$d
  }))
  all_levels <- unlist(lapply(parts, `[[`, "levels"))
  combined$transition <- factor(combined$transition, levels = rev(all_levels))
  combined$boundary <- factor(combined$boundary, levels = names(KEY_BOUNDARIES))
  combined
}

panel_build_perm_multi <- function(scenario, bearing_type, collapse) {
  parts <- lapply(names(KEY_BOUNDARIES), function(boundary_label) {
    bundle <- PANEL_BUNDLES[[KEY_BOUNDARIES[[boundary_label]]]]
    d <- panel_select_bundle(bundle, "perm", scenario, bearing_type, collapse)

    phases    <- strip_age(bundle$bin_labels)
    all_trans <- paste0(phases[-length(phases)], " → ", phases[-1])
    d$transition <- paste0(strip_age(d$bin_old), " → ", strip_age(d$bin_young))
    d$transition <- factor(d$transition, levels = rev(all_trans))
    d$group    <- factor(d$group, levels = names(bundle$volcano_groups))
    d$boundary <- boundary_label
    list(d = d, levels = all_trans)
  })

  combined <- do.call(rbind, lapply(parts, function(x) {
    x$d$transition <- as.character(x$d$transition)
    x$d
  }))
  all_levels <- unlist(lapply(parts, `[[`, "levels"))
  combined$transition <- factor(combined$transition, levels = rev(all_levels))
  combined$boundary <- factor(combined$boundary, levels = names(KEY_BOUNDARIES))
  combined
}

plot_permutation_panel_from_data <- function(d, subheading, show_legend = FALSE) {
  is_insufficient <- d$insufficient | is.na(d$U2)

  p_label <- ifelse(
    is.na(d$perm_p),          "p = —",
    ifelse(d$perm_p < 0.001,  "p < 0.001",
                              sprintf("p = %.3f", d$perm_p))
  )

  d$tile_label <- ifelse(
    is_insufficient,
    sprintf("n = %d, %d", d$n_old, d$n_young),
    sprintf("%s\n(n = %d, %d)", p_label, d$n_old, d$n_young)
  )

  ggplot(d, aes(x = group, y = transition, fill = perm_p)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = tile_label),
              colour = ifelse(is_insufficient, "grey45", "grey10"),
              size = 16, size.unit = "pt", lineheight = 1.1) +
    scale_y_discrete(expand = expansion(add = c(Y_GAP_BOTTOM, Y_GAP_TOP))) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(colours = colorRampPalette(c("#155A94", "#EFF6FC"))(5),
                        values  = seq(0, 1, length.out = 5)^3,
                        limits = c(0, 1), na.value = "grey85",
                        name = "Permutation\np-value",
                        guide = guide_colourbar(barheight = LEGEND_BARHEIGHT)) +
    labs(x = NULL, y = NULL, title = subheading) +
    thesis_theme +
    theme(
      panel.grid          = element_blank(),
      panel.grid.major.y  = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      plot.title      = element_text(size = 12.5, face = "plain", colour = "grey20",
                                      hjust = 0.5, lineheight = 1.15, margin = margin(b = 10)),
      plot.subtitle   = element_blank(),
      axis.line       = element_blank(),
      axis.ticks      = element_blank(),
      axis.ticks.length = unit(0, "pt"),
      axis.text.x     = element_text(margin = margin(t = 4)),
      axis.text.y     = element_blank(),
      plot.margin     = margin(t = 6, r = 10, b = 6, l = 8)
    )
}

build_combined_boundary_panels <- function(rows, scenario) {
  watson_data <- lapply(rows, function(r) panel_build_watson_multi(scenario, r$bearing_type, r$collapse))
  perm_data   <- lapply(rows, function(r) panel_build_perm_multi(scenario, r$bearing_type, r$collapse))
  xmax <- max(vapply(watson_data, function(d) max(d$U2, na.rm = TRUE), numeric(1)),
              CRIT$U2, na.rm = TRUE) * 1.08

  for (i in seq_along(rows)) {
    r <- rows[[i]]
    is_top <- (i == 1)
    row_tag <- sprintf("key_boundaries_row%d_%s_%s_%s", i, scenario, r$bearing_type,
                       if (r$collapse) "eruption" else "deposit")

    p_watson <- plot_watson_panel(watson_data[[i]], xmax,
                                   show_crit_labels = is_top,
                                   subheading       = r$subheading)
    save_square_panel(p_watson, sprintf("panel_%s_watson.png", row_tag))

    p_perm <- plot_permutation_panel_from_data(perm_data[[i]], subheading = r$subheading)
    save_square_panel(p_perm, sprintf("panel_%s_perm.png", row_tag))
  }
}

KEY_BOUNDARY_ROWS <- make_scheme_rows(NULL)   # bundle field unused by the _multi builders above

build_combined_boundary_panels(KEY_BOUNDARY_ROWS, "max")
build_combined_boundary_panels(KEY_BOUNDARY_ROWS, "min")

