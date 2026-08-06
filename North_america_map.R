
# ==============
# Load packages
# ==============

  library(dplyr)
  library(sf)
  library(ggplot2)
  library(rnaturalearth)
  library(ggrepel)

source("map_theme.R")

# ======================
# Load and prepare data
# ======================

# Paths
TEPHRA_CSV  <- "North_am_sites.csv"
KMZ_PATH    <- "194_GVP_Holocene_Volcanoes.kmz"
PLEIST_CSV  <- "GVP_Volcano_List_Pleistocene.csv"
KMZ_EXTRACT <- "kmz_extracted_na"
OUTPUT_PNG  <- "tephra_map_northamerica.png"

# Style constants
LAND_COLOUR   <- "#d9d3c2"
OCEAN_COLOUR  <- "#c8dce8"
COAST_COLOUR  <- "#5a5a5a"
BORDER_COLOUR <- "#5a5a5a"

# Lambert Conformal Conic projection
PROJ <- "+proj=lcc +lat_1=33 +lat_2=55 +lat_0=45 +lon_0=-90 +datum=WGS84"
#

# Load tephra study sites 
df <- read.csv(TEPHRA_CSV, fileEncoding = "UTF-8-BOM", check.names = FALSE) %>%
  select(Site, Latitude, Longitude)
for (col in c("Latitude", "Longitude")) {
  df[[col]] <- gsub("−", "-", as.character(df[[col]]))
  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
}
depos   <- df %>% filter(!is.na(Latitude), !is.na(Longitude))
n_sites <- n_distinct(depos$Site)

cat(sprintf("Loaded %d tephra sites\n", n_sites))


# Load Holocene volcanoes from KMZ
dir.create(KMZ_EXTRACT, showWarnings = FALSE)
utils::unzip(KMZ_PATH, exdir = KMZ_EXTRACT)

kml_file <- list.files(KMZ_EXTRACT, pattern = "\\.kml$", recursive = TRUE, full.names = TRUE)[1]
layers   <- st_layers(kml_file)$name
gdf <- bind_rows(lapply(layers, function(l) st_read(kml_file, layer = l, quiet = TRUE))) %>%
  st_zm(drop = TRUE) %>%
  st_set_crs(4326)

# Load Pleistocene volcanoes
pleist <- read.csv(PLEIST_CSV, encoding = "UTF-8", skip = 1, check.names = FALSE)
pleist$Latitude  <- suppressWarnings(as.numeric(pleist$Latitude))
pleist$Longitude <- suppressWarnings(as.numeric(pleist$Longitude))
pleist <- pleist %>% filter(!is.na(Latitude), !is.na(Longitude))
pleist_gdf <- st_as_sf(pleist, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)


# Laurentide/Cordilleran ice-sheet extent (Dalton et al. 2020)
ICE_DIR   <- "Dalton_ice_sheet"
ICE_FILES <- c(
  "1ka_0.91cal_Dalton_et_al_2020_QSR.shp"  = 0.91,
  "2ka_2cal_Dalton_et_al_2020_QSR.shp"     = 2.0,
  "4ka_4.5cal_Dalton_et_al_2020_QSR.shp"   = 4.5,
  "6.5ka_7.3cal_Dalton_et_al_2020_QSR.shp" = 7.3,
  "9ka_10.3cal_Dalton_et_al_2020_QSR.shp"  = 10.3,
  "11ka_12.8cal_Dalton_et_al_2020_QSR.shp" = 12.8,
  "13ka_15.5cal_Dalton_et_al_2020_QSR.shp" = 15.5
)

load_ice_isochrone <- function(file, cal_ka) {
  st_read(file.path(ICE_DIR, file), quiet = TRUE) %>%
    st_transform(4326) %>%   
    st_make_valid() %>%      
    st_union() %>%           
    st_as_sf() %>%
    mutate(cal_ka = cal_ka)
}

ice_isochrones <- bind_rows(Map(load_ice_isochrone, names(ICE_FILES), ICE_FILES)) %>%
  mutate(age_label = sprintf("%.1f ka cal BP", cal_ka),
         age_label = factor(age_label,
                             levels = sprintf("%.1f ka cal BP", sort(unique(cal_ka), decreasing = TRUE))))

ICE_AGE_PALETTE <- setNames(
  colorRampPalette(c("#4a86ac", "#082a45"))(nlevels(ice_isochrones$age_label)),
  levels(ice_isochrones$age_label))


# =====
# Plot
# =====

# 50m Natural Earth features 
land      <- ne_countries(scale = "medium", returnclass = "sf")
coastline <- ne_coastline(scale = "medium", returnclass = "sf")
borders   <- ne_download(scale = 50, type = "admin_0_boundary_lines_land",
                          category = "cultural", returnclass = "sf")
lakes     <- ne_download(scale = 50, type = "lakes", category = "physical", returnclass = "sf")

# Filter volcano data to map extent
NA_LON_MIN <- -130; NA_LON_MAX <- -50
NA_LAT_MIN <-    5; NA_LAT_MAX <-  75

hol_coords <- st_coordinates(gdf)
hol_na <- gdf[hol_coords[, "X"] >= NA_LON_MIN & hol_coords[, "X"] <= NA_LON_MAX &
                hol_coords[, "Y"] >= NA_LAT_MIN & hol_coords[, "Y"] <= NA_LAT_MAX, ]

pleist_na <- pleist_gdf %>%
  filter(between(Longitude, NA_LON_MIN, NA_LON_MAX), between(Latitude, NA_LAT_MIN, NA_LAT_MAX))

# Legend labels
pleist_label <- sprintf("Active Pleistocene volcanoes (n = %d)", nrow(pleist_na))
hol_label    <- sprintf("Active Holocene volcanoes (n = %d)", nrow(hol_na))
site_label   <- sprintf("Tephra study sites (n = %d)", n_sites)

p <- ggplot() +
  geom_sf(data = land, fill = LAND_COLOUR, colour = NA) +
  geom_sf(data = lakes, fill = OCEAN_COLOUR, colour = COAST_COLOUR, linewidth = 0.2, alpha = 0.75) +
  geom_sf(data = coastline, colour = COAST_COLOUR, linewidth = 0.2, alpha = 0.75) +
  geom_sf(data = borders, colour = BORDER_COLOUR, linewidth = 0.2, linetype = "dotted", alpha = 0.9) +
  # Ice-sheet isochrones 
  geom_sf(data = ice_isochrones, aes(colour = age_label), fill = NA, linewidth = 0.45) +
  scale_colour_manual(name = "Ice sheet margin", values = ICE_AGE_PALETTE) +
  # Pleistocene volcanoes
  geom_sf(data = pleist_na, aes(fill = pleist_label),
          shape = 24, size = 2.4, colour = COAST_COLOUR, stroke = 0.4) +
  # Holocene volcanoes
  geom_sf(data = hol_na, aes(fill = hol_label),
          shape = 24, size = 2.4, colour = COAST_COLOUR, stroke = 0.4) +
  # Tephra study sites
  geom_point(data = depos, aes(Longitude, Latitude, fill = site_label),
             shape = 21, size = 2.4, colour = COAST_COLOUR, stroke = 0.4) +
  scale_fill_manual(name = NULL, breaks = c(pleist_label, hol_label, site_label),
                     values = setNames(c("#d6d893", "#e0acca", "#874037"),
                                       c(pleist_label, hol_label, site_label))) +
  # Both legends stacked outside on the right
  guides(fill = guide_legend(position = "right", order = 1),
         colour = guide_legend(position = "right", order = 2,
                                keywidth = unit(0.7, "lines"), keyheight = unit(0.7, "lines"),
                                theme = theme(legend.key.spacing.y = unit(6, "pt")))) +
  coord_sf(crs = PROJ, default_crs = st_crs(4326),
           xlim = c(NA_LON_MIN, NA_LON_MAX), ylim = c(NA_LAT_MIN, NA_LAT_MAX)) +
  ggtitle("Tephra study sites in Northeast North America") +
  map_theme +
  theme(panel.background = element_rect(fill = OCEAN_COLOUR, colour = NA),
        panel.grid = element_line(colour = "grey60", linewidth = 0.2, linetype = "dashed"),
        axis.text.x       = element_text(size = 13, colour = "grey40"),
        axis.text.y       = element_text(size = 13, colour = "grey40"),
        axis.ticks        = element_line(colour = "grey40"),
        legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = NA),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.text = element_text(size = 13.5),
        legend.title = element_text(size = 13.5),
        legend.key.size = unit(1.1, "lines"))

# Save 
ggsave(OUTPUT_PNG, p, width = 12, height = 9, dpi = 600, bg = "white")
print(p)
cat(sprintf("Saved -> %s\n", OUTPUT_PNG))
