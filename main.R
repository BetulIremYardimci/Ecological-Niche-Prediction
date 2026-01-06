# ==============================================================================
# Ecological Niche Modeling of Spermophilus xanthoprymnus
# Main Analysis Script (main.R)
# ==============================================================================
# This script models current and future habitat suitability using a GAM
# Climate data: WorldClim v2.1 (historical 1970–2000) & CMIP6 MIROC6 SSP370
# Outputs: results/figures/*.png and results/rasters/*.tif
# ==============================================================================

suppressPackageStartupMessages({
  library(here)      # robust relative paths
  library(mgcv)      # GAM
  library(raster)    # raster I/O (kept for compatibility)
  library(dismo)     # randomPoints
  library(ggplot2)   # plotting
  library(maps)      # map_data("world")
  library(viridis)   # scale_fill_viridis_c
  library(raster)

})

# ------------------------------------------------------------------------------
#CONFIG
# ------------------------------------------------------------------------------
CFG <- list(
  seed = 42,
  background_n = 1000,
  occurrence_file = here("data", "raw", "occurrence.txt"),
  hist_dir = here("data", "raw", "historical_wc2"),
  fut_dir  = here("data", "raw", "future_wc2"),
  out_fig  = here("results", "figures"),
  out_rst  = here("results", "rasters")
)

# Historical layers (WorldClim bioclim)
HIST_FILES <- list(
  bio1  = "wc2.1_5m_bio_1.tif",   # Mean Annual Temperature
  bio11 = "wc2.1_5m_bio_11.tif",  # Mean Temperature of Coldest Quarter
  bio12 = "wc2.1_5m_bio_12.tif",  # Annual Precipitation
  bio19 = "wc2.1_5m_bio_19.tif"   # Precipitation of Coldest Quarter
)

# Future layers (CMIP6 MIROC6 SSP370 2021–2040)
# NOTE: These are used to approximate bio1/bio11/bio12/bio19.
# If you have future bioclim layers directly (bio_1, bio_11, etc.), load those instead.
FUT_FILES <- list(
  tmin = "wc2.1_5m_tmin_MIROC6_ssp370_2021-2040.tif",
  tmax = "wc2.1_5m_tmax_MIROC6_ssp370_2021-2040.tif",
  prec = "wc2.1_5m_prec_MIROC6_ssp370_2021-2040.tif",
  bioc = "wc2.1_5m_bioc_MIROC6_ssp370_2021-2040.tif"  # make sure this corresponds to bio19 or replace accordingly
)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
stop_if_missing <- function(paths) {
  miss <- paths[!file.exists(paths)]
  if (length(miss)) stop("Missing required file(s):\n- ", paste(miss, collapse = "\n- "), call. = FALSE)
}

ensure_dir(CFG$out_fig)
ensure_dir(CFG$out_rst)
set.seed(CFG$seed)

# ------------------------------------------------------------------------------
# 1) DATA LOADING
# ------------------------------------------------------------------------------
stop_if_missing(CFG$occurrence_file)
occurrences <- read.csv(CFG$occurrence_file, sep = "\t", stringsAsFactors = FALSE)

req_cols <- c("decimalLongitude", "decimalLatitude")
if (!all(req_cols %in% names(occurrences))) {
  stop("occurrence.txt must include columns: decimalLongitude, decimalLatitude", call. = FALSE)
}

occurrences <- occurrences[!is.na(occurrences$decimalLatitude) & !is.na(occurrences$decimalLongitude), ]
cat("Loaded", nrow(occurrences), "occurrence records\n")

# Historical climate variables
hist_paths <- file.path(CFG$hist_dir, unlist(HIST_FILES, use.names = FALSE))
stop_if_missing(hist_paths)

bio1  <- raster(file.path(CFG$hist_dir, HIST_FILES$bio1))
bio11 <- raster(file.path(CFG$hist_dir, HIST_FILES$bio11))
bio12 <- raster(file.path(CFG$hist_dir, HIST_FILES$bio12))
bio19 <- raster(file.path(CFG$hist_dir, HIST_FILES$bio19))

climate_stack <- stack(bio1, bio11, bio12, bio19)
names(climate_stack) <- names(HIST_FILES)
cat("Historical climate data loaded successfully\n")

# Future climate variables
fut_paths <- file.path(CFG$fut_dir, unlist(FUT_FILES, use.names = FALSE))
stop_if_missing(fut_paths)

future_tmin <- raster(file.path(CFG$fut_dir, FUT_FILES$tmin))
future_tmax <- raster(file.path(CFG$fut_dir, FUT_FILES$tmax))
future_prec <- raster(file.path(CFG$fut_dir, FUT_FILES$prec))
future_bioc <- raster(file.path(CFG$fut_dir, FUT_FILES$bioc))

# Derive future versions of selected historical variables (approximation)
future_bio1  <- (future_tmin + future_tmax) / 2
future_bio11 <- future_tmin
future_bio12 <- future_prec
future_bio19 <- future_bioc

future_projection <- stack(future_bio1, future_bio11, future_bio12, future_bio19)
names(future_projection) <- names(HIST_FILES)

# Align future stack to historical grid
future_projection <- resample(future_projection, climate_stack, method = "bilinear")
cat("Future climate data loaded and aligned\n")

# ------------------------------------------------------------------------------
# 2) DATA VISUALIZATION
# ------------------------------------------------------------------------------
png(file.path(CFG$out_fig, "occurrence_map.png"), width = 900, height = 650, res = 120)
plot(climate_stack[[1]], main = "S. xanthoprymnus - Occurrence Points")
points(occurrences$decimalLongitude, occurrences$decimalLatitude, col = "red", pch = 20, cex = 0.7)
dev.off()
cat("Occurrence map saved to results/figures/occurrence_map.png\n")

# ------------------------------------------------------------------------------
# 3) MODEL PREPARATION
# ------------------------------------------------------------------------------
presence_data <- extract(climate_stack, occurrences[, c("decimalLongitude", "decimalLatitude")])
presence_data <- cbind(as.data.frame(presence_data), presence = 1)

background <- randomPoints(climate_stack, CFG$background_n)
background_data <- extract(climate_stack, background)
background_data <- cbind(as.data.frame(background_data), presence = 0)

model_data <- rbind(presence_data, background_data)
model_data <- na.omit(model_data)
model_data <- as.data.frame(model_data)

cat("Model data prepared:", nrow(model_data), "total samples\n")
cat("  - Presence:", sum(model_data$presence == 1), "\n")
cat("  - Background:", sum(model_data$presence == 0), "\n")

# ------------------------------------------------------------------------------
# 4) GAM MODEL FITTING
# ------------------------------------------------------------------------------
# GAM chosen for flexibility in capturing non-linear species–environment relationships

gam_model <- gam(
  presence ~ s(bio1) + s(bio11) + s(bio12) + s(bio19),
  data = model_data,
  family = binomial()
)

cat("\n=== GAM Model Summary ===\n")
print(summary(gam_model))

# ------------------------------------------------------------------------------
# 5) CURRENT SUITABILITY PREDICTION
# ------------------------------------------------------------------------------
current_prediction <- predict(climate_stack, gam_model, type = "response")

writeRaster(
  current_prediction,
  filename = file.path(CFG$out_rst, "current_suitability.tif"),
  overwrite = TRUE
)

png(file.path(CFG$out_fig, "current_suitability.png"), width = 1100, height = 650, res = 120)
plot(current_prediction, main = "Current Habitat Suitability (1970–2000)")
dev.off()
cat("Current suitability saved to results/rasters/current_suitability.tif and results/figures/current_suitability.png\n")

# ------------------------------------------------------------------------------
# 6) FUTURE SUITABILITY PREDICTION
# ------------------------------------------------------------------------------
future_prediction <- predict(future_projection, gam_model, type = "response")

writeRaster(
  future_prediction,
  filename = file.path(CFG$out_rst, "future_suitability.tif"),
  overwrite = TRUE
)

png(file.path(CFG$out_fig, "future_suitability.png"), width = 1100, height = 650, res = 120)
plot(future_prediction, main = "Future Habitat Suitability (2021–2040, MIROC6 SSP370)")
dev.off()
cat("Future suitability saved to results/rasters/future_suitability.tif and results/figures/future_suitability.png\n")

# ------------------------------------------------------------------------------
# 7) WORLD MAP VISUALIZATION (GGPLOT2)
# ------------------------------------------------------------------------------
# Convert predictions to data frames for ggplot
current_df <- as.data.frame(rasterToPoints(current_prediction))
future_df  <- as.data.frame(rasterToPoints(future_prediction))

colnames(current_df) <- c("x", "y", "suitability")
colnames(future_df)  <- c("x", "y", "suitability")

world <- map_data("world")

p_current <- ggplot() +
  geom_tile(data = current_df, aes(x = x, y = y, fill = suitability)) +
  geom_path(data = world, aes(x = long, y = lat, group = group), color = "black", linewidth = 0.2) +
  scale_fill_viridis_c(name = "Suitability") +
  ggtitle("Current Habitat Suitability (1970–2000)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(CFG$out_fig, "current_worldmap.png"), p_current, width = 12, height = 7, dpi = 300)
cat("Current world map saved to results/figures/current_worldmap.png\n")

p_future <- ggplot() +
  geom_tile(data = future_df, aes(x = x, y = y, fill = suitability)) +
  geom_path(data = world, aes(x = long, y = lat, group = group), color = "black", linewidth = 0.2) +
  scale_fill_viridis_c(name = "Suitability") +
  ggtitle("Future Habitat Suitability (2021–2040, MIROC6 SSP370)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(CFG$out_fig, "future_worldmap.png"), p_future, width = 12, height = 7, dpi = 300)
cat("Future world map saved to results/figures/future_worldmap.png\n")

# ------------------------------------------------------------------------------
# 8) SUITABILITY CHANGE ANALYSIS (Future - Current)
# ------------------------------------------------------------------------------
change_raster <- future_prediction - current_prediction

writeRaster(
  change_raster,
  filename = file.path(CFG$out_rst, "suitability_change_future_minus_current.tif"),
  overwrite = TRUE
)

png(file.path(CFG$out_fig, "suitability_change.png"), width = 1100, height = 650, res = 120)
plot(change_raster, main = "Change in Suitability (Future - Current)")
dev.off()

cat("Suitability change saved to results/rasters/suitability_change_future_minus_current.tif and results/figures/suitability_change.png\n")

# Summary statistics
vals <- values(change_raster)
cat("\n=== Suitability Change Summary ===\n")
cat("Mean change:", mean(vals, na.rm = TRUE), "\n")
cat("Max increase:", max(vals, na.rm = TRUE), "\n")
cat("Max decrease:", min(vals, na.rm = TRUE), "\n")

cat("\n✓ Analysis completed successfully!\n")
cat("Outputs saved under results/figures and results/rasters\n")
