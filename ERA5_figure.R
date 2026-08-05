
# ===============================
# ~~~ ERA5 reanalysis figure ~~~
# ===============================

# Modern (1991-2020) DJF vs JJA mean 250 hPa circulation over the Mediterranean from ERA5 monthly-mean reanalysis
# Requires the ERA5 250 hPa monthly-mean u/v wind component netCDF file 

# Contents:
# ------------
#   1. Load packages
#   2. Configuration
#   3. Load and prepare ERA5 data
#   4. Plot


# =================
# 1. Load packages
# =================

  library(ncdf4)
  library(data.table)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)


# =================
# 2. Configuration
# =================

# R theme setup from TraCE script 
src <- readLines("Trace21k_simulation.R")
main_start <- grep("^## -+ MAIN", src)[1]
eval(parse(text = src[seq_len(main_start - 1)]), envir = .GlobalEnv)

ERA5_FILE      <- "5a8e1fc88438bd1c89aca37ccaad82de.nc"
ITALY_LON      <- 12.5
ITALY_LAT      <- 42
MAP_LON        <- c(-6, 36)      
MAP_LAT        <- c(27, 50)   
ARROW_SCALE_M  <- 100000         
THIN_STRIDE    <- 6              

COAST_COLOUR  <- "#2b2b2b"
ARROW_COLOUR  <- "grey15"

SPEED_BREAKS  <- c(0, 10, 15, 20, 25, 30, 35, 40, 45, Inf)
SPEED_COLOURS <- c("#ffffff", "#f6e8c3", "#f0c987", "#f0a04b", "#e8791f",
                    "#d6432b", "#c2185b", "#9c27b0", "#5e0e8f")

MAP_CRS <- st_crs(sprintf(
  "+proj=eqc +lat_ts=%f +lon_0=%f +datum=WGS84 +units=m +no_defs",
  mean(MAP_LAT), ITALY_LON))


# ==============================
# 3. Load and prepare ERA5 data
# ==============================

nc  <- nc_open(ERA5_FILE)
lon <- ncvar_get(nc, "longitude")
lat <- ncvar_get(nc, "latitude")
vt  <- ncvar_get(nc, "valid_time")
u   <- ncvar_get(nc, "u")
v   <- ncvar_get(nc, "v")
nc_close(nc)

times  <- as.POSIXct(vt, origin = "1970-01-01", tz = "UTC")
months <- as.integer(format(times, "%m"))
cat(sprintf("loaded %d monthly steps, %s to %s\n",
            length(vt), format(min(times), "%Y-%m"), format(max(times), "%Y-%m")))

season_of <- function(m) {
  ifelse(m %in% c(12, 1, 2), "DJF", ifelse(m %in% c(6, 7, 8), "JJA", NA))
}
season <- season_of(months)

seasonal_mean <- function(arr, sel) apply(arr[, , sel, drop = FALSE], c(1, 2), mean, na.rm = TRUE)

lon_idx <- seq(1, length(lon), by = THIN_STRIDE)
lat_idx <- seq(1, length(lat), by = THIN_STRIDE)

make_season_table <- function(s) {
  sel <- which(season == s)
  cat(sprintf("  %s: averaging %d months\n", s, length(sel)))
  ubar <- seasonal_mean(u, sel)
  vbar <- seasonal_mean(v, sel)
  CJ(li = lon_idx, lj = lat_idx)[
    , `:=`(lon = lon[li], lat = lat[lj],
           u = ubar[cbind(li, lj)], v = vbar[cbind(li, lj)], season = s)]
}

cat("computing seasonal composites...\n")
clim <- rbindlist(lapply(c("DJF", "JJA"), make_season_table))
clim[, speed := sqrt(u^2 + v^2)]
clim[, season := factor(season, levels = c("DJF", "JJA"))]


xy0 <- st_coordinates(st_transform(
  st_as_sf(clim, coords = c("lon", "lat"), crs = 4326, remove = FALSE), MAP_CRS))
clim[, `:=`(x0 = xy0[, 1], y0 = xy0[, 2])]
clim[, `:=`(xend = x0 + (u / speed) * ARROW_SCALE_M,
            yend = y0 + (v / speed) * ARROW_SCALE_M)]

make_speed_grid <- function(s) {
  sel <- which(season == s)
  ubar <- seasonal_mean(u, sel)
  vbar <- seasonal_mean(v, sel)
  CJ(li = seq_along(lon), lj = seq_along(lat))[
    , `:=`(lon = lon[li], lat = lat[lj],
           speed = sqrt(ubar[cbind(li, lj)]^2 + vbar[cbind(li, lj)]^2), season = s)]
}
speed_grid <- rbindlist(lapply(c("DJF", "JJA"), make_speed_grid))
speed_grid[, season := factor(season, levels = c("DJF", "JJA"))]
speed_grid[, speed_band := cut(speed, SPEED_BREAKS, right = FALSE)]

xyg <- st_coordinates(st_transform(
  st_as_sf(speed_grid, coords = c("lon", "lat"), crs = 4326, remove = FALSE), MAP_CRS))
speed_grid[, `:=`(x = xyg[, 1], y = xyg[, 2])]


# ========
# 4. Plot
# ========

box_pts <- st_as_sf(data.frame(
  lon = c(MAP_LON[1], MAP_LON[2], MAP_LON[1], MAP_LON[2], ITALY_LON, ITALY_LON),
  lat = c(MAP_LAT[1], MAP_LAT[1], MAP_LAT[2], MAP_LAT[2], MAP_LAT[1], MAP_LAT[2])),
  coords = c("lon", "lat"), crs = 4326)
box_xy  <- st_coordinates(st_transform(box_pts, MAP_CRS))
PROJ_XLIM <- range(box_xy[, 1])
PROJ_YLIM <- range(box_xy[, 2])

land <- st_transform(st_crop(ne_countries(scale = "medium", returnclass = "sf"),
                              xmin = MAP_LON[1] - 2, xmax = MAP_LON[2] + 2,
                              ymin = MAP_LAT[1] - 2, ymax = MAP_LAT[2] + 2),
                      MAP_CRS)

p <- ggplot() +
  geom_raster(data = speed_grid, aes(x, y, fill = speed_band)) +
  geom_sf(data = land, fill = NA, colour = COAST_COLOUR, linewidth = 0.2) +
  geom_segment(data = clim, aes(x0, y0, xend = xend, yend = yend),
               colour = ARROW_COLOUR,
               arrow = arrow(length = unit(0.06, "cm"), type = "closed"),
               linewidth = 0.3) +
  scale_fill_manual(name = expression("Wind speed (m s"^-1*")"),
                     values = setNames(SPEED_COLOURS, levels(speed_grid$speed_band)),
                     labels = c("0-10", "10-15", "15-20", "20-25", "25-30",
                                "30-35", "35-40", "40-45", "45+"),
                     drop = FALSE) +
  scale_x_continuous(breaks = seq(10 * ceiling(MAP_LON[1] / 10), 10 * floor(MAP_LON[2] / 10), 10)) +
  scale_y_continuous(breaks = seq(10 * ceiling(MAP_LAT[1] / 10), 10 * floor(MAP_LAT[2] / 10), 10)) +
  coord_sf(crs = MAP_CRS, xlim = PROJ_XLIM, ylim = PROJ_YLIM, expand = FALSE) +
  facet_wrap(~season, ncol = 1) +
  labs(x = NULL, y = NULL, title = "ERA5 250 hPa mean circulation, 1991-2020") +
  FIG_THEME +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.3),
        panel.spacing.y = unit(1.5, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(size = 13),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

ggsave("era5_modern_circulation.png", p, width = 8.5, height = 12, dpi = 400)
cat("\nwrote era5_modern_circulation.png\n")
