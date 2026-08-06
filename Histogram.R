
# ===============
# Load packages
# ===============

  library(dplyr)
  library(sf)
  library(ggplot2)
  library(rnaturalearth)
  library(readxl)
  library(patchwork)
  library(scico)
  library(ggrepel)
  library(terra)
  library(tidyterra)
  library(ggnewscale)
  library(elevatr)


# ==========
# Load data
# ==========

tephra <- read_excel("Mediterranean_tephra_data.xlsx", sheet = "Tephra data")
sitespans <- read_excel("Mediterranean_tephra_data.xlsx", sheet = "Site spans")


# ===========
# Site spans
# ===========

# Site span data frame
sitespans <- rename(sitespans,
  site                 = `Site`,
  min_coverage         = `Minimum site coverage (cal yr BP)`,
  max_coverage         = `Maximum site coverage (cal yr BP)`
) %>%
  mutate(
    min_coverage = suppressWarnings(as.numeric(min_coverage)),
    max_coverage = suppressWarnings(as.numeric(max_coverage))
  )

# Define age grid, spanning 0-25 ka in 1 yr increments
site_cov_yearly <- data.frame(bin_center = seq(0, 25000, by = 1)) %>%
  mutate(
    yr_start = bin_center,
    yr_end   = bin_center + 1
  )

# Function to count sites covering each calendar yr
count_sites <- function(start, end, sitespans) {
  count <-
    sum(
      sitespans$min_coverage <= end &
        sitespans$max_coverage >= start,
      na.rm = TRUE
)}

# Calculate counts for each yr
site_cov_yearly <- site_cov_yearly %>%
  rowwise() %>%
  mutate(
    site_count = count_sites(yr_start, yr_end, sitespans),
  ) %>%
  ungroup()

# Aggregate the per-yr counts into the same 500-yr bins used for the tephra
# deposits (0-500, 500-1000, etc.), so the site-coverage line is directly
# comparable bin-for-bin with the bars
site_cov <- data.frame(bin_start = seq(0, 24500, by = 500)) %>%
  mutate(
    bin_end    = bin_start + 500 - 0.1,
    bin_center = bin_start + 250
  ) %>%
  rowwise() %>%
  mutate(
    site_count = mean(site_cov_yearly$site_count[site_cov_yearly$bin_center >= bin_start &
                                                    site_cov_yearly$bin_center <= bin_end])
  ) %>%
  ungroup()

print(site_cov)


site_cov_area <- data.frame(
  bin_center = as.vector(rbind(site_cov$bin_start, site_cov$bin_end)),
  site_count = rep(site_cov$site_count, each = 2)
)

# ==========================================================
# Temporal distribution of tephra deposits by source region
# ==========================================================

# Groups every "Source volcano" entry into a broader source region and plots
# the number of tephra deposits per 500-yr bin by median age

# Grouping/attribution logic mirrors Dispersal_analysis.R
# "minimum" scenario = certain attributions, primary sources only
# "maximum" scenario = uncertain "(?)" attributions and secondary populations of
# deposits containing multiple populations also included

showtext::showtext_opts(dpi = 400)

# Columns needed to reuse the assign_group()/expand_sources() logic below
tephra_groups <- rename(tephra,
  source_volcano        = `Source volcano`,
  source_eruption       = `Source eruption`,
  age_median            = `Median age (cal yr BP)`,
  distance_km_secondary = `Great circle distance between source and deposit (km) - secondary`,
  bearing_secondary     = `Angle between source and deposit (bearing in degrees) - secondary`
)

# Volcano grouping scheme:
VOLCANO_GROUPS <- list(
  "Campanian" = list(
    exact   = c("Campi Flegrei", "Ischia", "Somma-Vesuvius", "Procida-Vivara"),
    pattern = "Campi Flegrei|Ischia|Somma.Vesuvius|Procida|Campanian|Phlegraean"
  ),
  "Aeolian" = list(
    exact   = c("Lipari", "Vulcano", "Palinuro Seamount"),
    pattern = "Lipari|Vulcano|Aeolian|Palinuro"
  ),
  "Etna" = list(
    exact   = c("Etna"),
    pattern = "Etna"
  ),
  "Aegean" = list(
    exact   = c("Santorini"),
    pattern = "Santorini"
  ),
  "Anatolian" = list(
    exact   = c("Erciyes Dagi", "Acigöl", "Süphan Dagi"),
    pattern = "Erciyes|Acigöl|Anatolian|Nemrut|Süphan"
  ),
  "Icelandic" = list(
    exact   = c("Katla", "Askja"),
    pattern = "Katla|Askja"
  ),
  "Eifel" = list(
    exact   = c("East Eifel volcanic field"),
    pattern = "Eifel"
  ),
  "Unknown" = list(
    exact   = c("Unknown"),
    pattern = "Unknown"
  )
)

# Same function as assign_group() in Dispersal_analysis.R:
assign_group <- function(volcano_names, scenario = c("min", "max")) {
  scenario <- match.arg(scenario)
  vapply(volcano_names, function(v) {
    for (gname in names(VOLCANO_GROUPS)) {
      g <- VOLCANO_GROUPS[[gname]]
      matched <- if (scenario == "min") v %in% g$exact else grepl(g$pattern, v, ignore.case = TRUE)
      if (matched) return(gname)
    }
    return(NA_character_)
  }, character(1))
}

expand_sources <- function(data) {
  parts <- strsplit(as.character(data$source_volcano), ",", fixed = TRUE)
  n1 <- trimws(vapply(parts, function(p) p[1], character(1)))
  n2 <- trimws(vapply(parts,
                      function(p) if (length(p) >= 2) p[2] else NA_character_,
                      character(1)))

  prim <- data
  prim$source_volcano <- n1
  prim$source_rank    <- "primary"

  sec_ok <- !is.na(n2) & !is.na(data$distance_km_secondary) & !is.na(data$bearing_secondary)

  if (any(sec_ok)) {
    sec <- data[sec_ok, ]
    sec$source_volcano <- n2[sec_ok]
    sec$source_rank    <- "secondary"
    out <- rbind(prim, sec[, names(prim)])
  } else {
    out <- prim
  }
  out
}

# Bin tephras by median age into 500-yr bins (0-500, 500-1000, etc.)
deposits <- tephra_groups %>%
  filter(!is.na(age_median)) %>%
  mutate(bin_start  = floor(age_median / 500) * 500,
         bin_center = bin_start + 250) %>%
  expand_sources()

deposits$group_min <- assign_group(deposits$source_volcano, "min")
deposits$group_max <- assign_group(deposits$source_volcano, "max")

# Minimum scenario data frame
df_min <- deposits %>% filter(!is.na(group_min), source_rank == "primary") %>% rename(group = group_min)
# Maximum scenario data frame
df_max <- deposits %>% filter(!is.na(group_max))                            %>% rename(group = group_max)

cat(sprintf("\nMinimum scenario: %d deposits\n", nrow(df_min)))
cat(sprintf("Maximum scenario: %d deposits\n\n", nrow(df_max)))

cat("Group counts (minimum scenario):\n")
print(table(df_min$group))

cat("\nGroup counts (maximum scenario):\n")
print(table(df_max$group))

# Legend
region_levels <- c("Campanian", "Aeolian", "Etna", "Aegean",
                   "Anatolian", "Icelandic", "Eifel", "Unknown")

# Colour scheme
region_colour_palette <- setNames(
  c("#B8DEC3", "#97552B", "#5F4C81", "#20374dff",
    "#C77DAA", "#4e7048ff", "#b23a3aff", "#4D7FB3"),
  region_levels
)

# Final figure size
FIG_WIDTH_IN  <- 11
FIG_HEIGHT_IN <- 8

# Builds one stacked histogram (per source region) for a given scenario's data
make_source_histogram <- function(data, scenario_label) {
  region_counts <- data %>%
    count(group, name = "n_region") %>%
    mutate(label = sprintf("%s (n = %d)", group, n_region))

  levels_present <- region_levels[region_levels %in% region_counts$group]

  data <- data %>%
    left_join(region_counts, by = "group") %>%
    mutate(label = factor(label, levels = region_counts$label[match(levels_present, region_counts$group)]))

  region_colours <- setNames(
    region_colour_palette[levels_present],
    region_counts$label[match(levels_present, region_counts$group)]
  )

  n_deposits <- nrow(data)

  # Eruptions with > 5 deposits are labelled in a single row above the bars
  bar_top_max <- max(data %>% count(bin_center, name = "bar_top") %>% pull(bar_top))

  y_axis_max  <- ceiling(max(bar_top_max, site_cov$site_count) / 10) * 10
  y_breaks    <- seq(0, y_axis_max, by = 10)
  y_minor     <- seq(0, y_axis_max, by = 5)

  major_breaks <- seq(0, 25000, 5000)
  minor_breaks <- seq(0, 25000, 1000)

  tick_major_len <- 0.018 * y_axis_max
  tick_minor_len <- 0.009 * y_axis_max
  tick_label_y   <- -0.05 * y_axis_max
  title_y        <- -0.10 * y_axis_max

  label_bottom_y <- -0.18 * y_axis_max
  label_top_y    <- 0.95 * y_axis_max

  erup_labels <- data %>%
    filter(!is.na(source_eruption), source_eruption != "", source_eruption != "Unknown") %>%
    group_by(source_eruption) %>%
    summarise(n_erup = n(), erup_age = mean(age_median),
              group = dplyr::first(group), .groups = "drop") %>%
    filter(n_erup > 5) %>%
    arrange(desc(erup_age))

  spread_close_labels <- function(ages, min_gap) {
    ord <- order(-ages)
    x   <- ages[ord]
    n   <- length(x)
    i <- 1
    while (i < n) {
      j <- i
      while (j < n && (x[j] - x[j + 1]) < min_gap) j <- j + 1
      if (j > i) {
        count   <- j - i + 1
        centre  <- mean(x[i:j])
        offsets <- seq((count - 1) / 2, -(count - 1) / 2, length.out = count) * min_gap
        x[i:j]  <- centre + offsets
      }
      i <- j + 1
    }
    out <- numeric(n)
    out[ord] <- x
    out
  }
  erup_labels$label_x <- spread_close_labels(erup_labels$erup_age, min_gap = 1100)

  # Edge tick labels
  tick_label_df <- data.frame(
    x     = major_breaks,
    label = format(major_breaks, big.mark = ","),
    hjust = ifelse(major_breaks == max(major_breaks), 0,
             ifelse(major_breaks == min(major_breaks), 1, 0.5))
  )

  manual_axis_layer <- list(
    annotate("segment", x = 25000, xend = 0, y = 0, yend = 0,
             linewidth = 0.6, colour = "grey20"),
    annotate("segment", x = 25000, xend = 25000, y = 0, yend = y_axis_max,
             linewidth = 0.6, colour = "grey20"),
    annotate("segment", x = major_breaks, xend = major_breaks, y = 0, yend = -tick_major_len,
             linewidth = 0.4, colour = "grey20"),
    annotate("segment", x = minor_breaks, xend = minor_breaks, y = 0, yend = -tick_minor_len,
             linewidth = 0.3, colour = "grey20"),
    geom_text(data = tick_label_df, aes(x = x, y = tick_label_y, label = label, hjust = hjust),
              inherit.aes = FALSE, size = 10.5, size.unit = "pt", colour = "grey20", vjust = 1),
    annotate("text", x = 0, y = title_y, label = "Age (cal yrs BP)", hjust = 1, vjust = 1,
             size = 11.5, size.unit = "pt", colour = "grey20")
  )

  eruption_label_layer <- geom_text(data = erup_labels,
              aes(x = label_x, y = label_top_y, label = source_eruption),
              inherit.aes = FALSE, angle = 90, hjust = 1, vjust = 0.5,
              size = 3.6, lineheight = 0.85)

  ggplot(data, aes(x = bin_center)) +
    geom_area(data = site_cov_area, aes(x = bin_center, y = pmin(site_count, y_axis_max), fill = "Total sites with\nactive records"),
              inherit.aes = FALSE) +   # drawn first, so the bars sit on top of it unobstructed; capped so it never bleeds above the axis's real-data range into the label zone
    scale_fill_manual(name = NULL, values = c("Total sites with\nactive records" = "#f5f1f1"),
                       guide = guide_legend(order = 2)) +
    ggnewscale::new_scale_fill() +
    geom_bar(aes(fill = label), width = 500, colour = "grey30", linewidth = 0.15) +   # Matches the 500-yr bin width, so bars sit edge-to-edge
    manual_axis_layer +
    eruption_label_layer +
    scale_fill_manual(name = "Volcanic source", values = region_colours,
                       guide = guide_legend(order = 1)) +
    scale_x_reverse(limits = c(25000, 0), expand = c(0, 0)) +
    scale_y_continuous(name = "Number of tephra deposits and sites \n with active records per bin",
                       breaks = y_breaks, minor_breaks = y_minor,
                       guide = guide_axis(minor.ticks = TRUE),
                       expand = c(0, 0)) +
    coord_cartesian(ylim = c(label_bottom_y, y_axis_max)) +
    labs(title = sprintf("Temporal distribution of tephra deposits by source region - %s scenario (n=%d)",
                         scenario_label, n_deposits)) +
    thesis_theme +
    theme(panel.grid.major.y      = element_blank(),
          axis.line.x             = element_blank(),
          axis.ticks.x            = element_blank(),
          axis.text.x             = element_blank(),
          axis.title.x            = element_blank(),
          axis.line.y             = element_blank(),
          axis.ticks.y            = element_line(linewidth = 0.4),
          plot.title              = element_text(margin = margin(b = 0.3, unit = "cm")),
          plot.margin             = margin(t = 0.3, r = 0.3, b = 0.3, l = 0.3, unit = "cm"),
          axis.text.y             = element_text(size = 10.5, colour = "grey20"),
          axis.title.y            = element_text(size = 11.5, colour = "grey20", margin = margin(r = 18)),
          axis.minor.ticks.length = unit(0.08, "cm"),
          axis.ticks.length       = unit(0.16, "cm"),
          legend.position         = "right",
          legend.title            = element_text(size = 10.5, face = "bold"),
          legend.text             = element_text(size = 10.5),
          legend.key.size         = unit(0.8, "lines"),
          legend.key.spacing.y    = unit(5, "pt"),
          legend.margin           = margin(3, 5, 3, 5),
          legend.spacing.y        = unit(0.5, "cm"))
}

p3_min <- make_source_histogram(df_min, "Minimum")
p3_max <- make_source_histogram(df_max, "Maximum")

print(p3_min)
print(p3_max)

ggsave("tephra_temporal_sources_min.png", p3_min,
       width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN, units = "in", dpi = 400)
ggsave("tephra_temporal_sources_max.png", p3_max,
       width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN, units = "in", dpi = 400)
