### Drone Multispectral Image Co-registration in R
# Author: Sofiya Tumanova
# Data used: DJI multispectral TIFFs (Red, Green, NIR) provided by prof. Duccio Rocchini PhD
# Description: Interactive co-registration of multispectral drone images
# Requires: R >= 4.2, terra package

# -----------------------------
# 0. Setup
# -----------------------------
library(terra)

# -----------------------------
# 1. User parameters (adjust based on data)
# -----------------------------
DATA_DIR <- "data"          # Folder where TIFF bands are stored
ZOOM_WINDOW <- 150          # Pixels for zoomed GCP selection (adjust accordingly)
OUTPUT_DIR <- "output"      # Folder to save aligned stacks

# Ensure output folder exists
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# -----------------------------
# 2. Load bands
# -----------------------------
r <- rast(file.path(DATA_DIR, "DJI_20251115110614_0003_MS_R.TIF"))
g <- rast(file.path(DATA_DIR, "DJI_20251115110614_0003_MS_G.TIF"))
nir <- rast(file.path(DATA_DIR, "DJI_20251115110614_0003_MS_NIR.TIF"))

# Optional: flip vertically if needed
r <- flip(r, "vertical")
g <- flip(g, "vertical")
nir <- flip(nir, "vertical")

# -----------------------------
# 3. Functions
# -----------------------------

# 3a. Select GCP with zoomed-in plot
select_gcp_zoom <- function(band, approx_x, approx_y, zoom_window, band_name="Band") {
  # Ensure window stays within raster
  xmin <- max(approx_x - zoom_window, 1)
  xmax <- min(approx_x + zoom_window, ncol(band))
  ymin <- max(approx_y - zoom_window, 1)
  ymax <- min(approx_y + zoom_window, nrow(band))
  
  patch <- crop(band, ext(xmin, xmax, ymin, ymax))
  plot(patch, main=paste("Zoomed-in", band_name, ": click on GCP"))
  gcp <- locator(1)
  
  x <- round(gcp$x + xmin - 1)
  y <- round(gcp$y + ymin - 1)
  cat(sprintf("GCP in %s: x=%d, y=%d\n", band_name, x, y))
  
  return(c(x=x, y=y))
}

# 3b. Align bands based on pixel shifts
align_band <- function(base_band, target_band, gcp_base, gcp_target) {
  dx <- gcp_base["x"] - gcp_target["x"]
  dy <- gcp_base["y"] - gcp_target["y"]
  cat(sprintf("Shift for target band: dx=%d, dy=%d\n", dx, dy))
  return(shift(target_band, dx=dx, dy=dy))
}

# -----------------------------
# 4. Select GCPs interactively
# -----------------------------
approx_x <- ncol(r)/2
approx_y <- nrow(r)/2

gcp_red <- select_gcp_zoom(r, approx_x, approx_y, ZOOM_WINDOW, "Red")
gcp_green <- select_gcp_zoom(g, approx_x, approx_y, ZOOM_WINDOW, "Green")
gcp_nir <- select_gcp_zoom(nir, approx_x, approx_y, ZOOM_WINDOW, "NIR")

# -----------------------------
# 5. Align bands
# -----------------------------
g_aligned <- align_band(r, g, gcp_red, gcp_green)
nir_aligned <- align_band(r, nir, gcp_red, gcp_nir)

# -----------------------------
# 6. Crop to common overlap
# -----------------------------
ext_common <- intersect(ext(r), intersect(ext(g_aligned), ext(nir_aligned)))
r_c <- crop(r, ext_common)
g_c <- crop(g_aligned, ext_common)
nir_c <- crop(nir_aligned, ext_common)

# -----------------------------
# 7. Stack bands and save
# -----------------------------
stack_aligned <- c(r_c, g_c, nir_c)
names(stack_aligned) <- c("Red", "Green", "NIR")

# Save aligned stack
writeRaster(stack_aligned, file.path(OUTPUT_DIR, "aligned_stack.tif"), overwrite=TRUE)

# False-color visualization: NIR-Red-Green
plotRGB(stack_aligned, r=3, g=1, b=2, stretch="hist", main="False Color (Aligned)")

# -----------------------------
# 8. Overlap check
# -----------------------------
plot(r_c, col=gray.colors(256), main="Red band overlay")
plot(g_c, col=rgb(0,1,0,0.3), add=TRUE)
plot(nir_c, col=rgb(0,0,1,0.3), add=TRUE)

# -----------------------------
# 9. Side-by-side visualization
# -----------------------------
# Full 3-band comparison
par(mfrow=c(1,2))
stack_before <- c(r, g, nir)
names(stack_before) <- c("Red", "Green", "NIR")
plotRGB(stack_before, r=3, g=1, b=2, stretch="hist", main="Before Coregistration")
plotRGB(stack_aligned, r=3, g=1, b=2, stretch="hist", main="After Coregistration")
par(mfrow=c(1,1))

# R + G comparison
par(mfrow=c(1,2))
rg_before <- c(r, g)
rg_after <- c(r_c, g_c)
plotRGB(rg_before, r=1, g=2, b=1, stretch="hist", main="R + G BEFORE")
plotRGB(rg_after, r=1, g=2, b=1, stretch="hist", main="R + G AFTER")
par(mfrow=c(1,1))

# R + NIR comparison
par(mfrow=c(1,2))
rn_before <- c(r, nir)
rn_after <- c(r_c, nir_c)
plotRGB(rn_before, r=1, g=2, b=1, stretch="hist", main="R + NIR BEFORE")
plotRGB(rn_after, r=1, g=2, b=1, stretch="hist", main="R + NIR AFTER")
par(mfrow=c(1,1))

# -----------------------------
# 10. Done
# -----------------------------
cat("Co-registration complete. Aligned stack saved in", OUTPUT_DIR, "\n")
