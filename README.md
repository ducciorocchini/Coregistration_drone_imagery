# Coregistration of Multispectral Drone Images in R

> Script: coreg_correl.R

This document describes a simple, correlation-based workflow to **coregister multispectral drone imagery** (Green, Red, NIR) using **R** and the **terra** package.

The approach aligns bands by searching for the spatial shift that maximizes pixel-wise correlation with a reference band (here: **Red**).

---

## Requirements

* R (>= 4.0 recommended)
* Packages:

  * `terra`

```r
install.packages("terra")
```

---

## Input Data

The script assumes three single-band GeoTIFF files acquired by a multispectral drone sensor:

* Green band
* Red band (used as reference)
* Near-Infrared (NIR) band

Example filenames:

```text
DJI_20251115110614_0003_MS_G.tiff
DJI_20251115110614_0003_MS_R.tiff
DJI_20251115110614_0003_MS_NIR.tiff
```

All files are expected to be in the same directory.

---

## Load and Visualize the Original Data

```r
library(terra)
setwd("~/Downloads/")

gre <- rast("DJI_20251115110614_0003_MS_G.tiff")
red <- rast("DJI_20251115110614_0003_MS_R.tiff")
nir <- rast("DJI_20251115110614_0003_MS_NIR.tiff")

ms <- c(gre, red, nir)
plotRGB(ms, r = 3, g = 2, b = 1, stretch = "lin")
```

This RGB composite highlights any spatial misalignment between bands.

---

## Convert Raster Bands to Matrices

For pixel-wise correlation analysis, raster layers are converted to matrices.

```r
band_names <- c("Green", "Red", "NIR")

band_matrices <- lapply(1:nlyr(ms), function(i) {
  as.matrix(ms[[i]], wide = TRUE)
})

names(band_matrices) <- band_names
```

---

## Correlation-Based Alignment Function

This function searches for the x/y shift that maximizes correlation between a reference matrix and a target matrix.

```r
align_band_simple <- function(ref_mat, target_mat, max_shift = 20) {
  best_cor <- -Inf; best_dx <- 0; best_dy <- 0
  rows <- nrow(ref_mat); cols <- ncol(ref_mat)

  for (dx in -max_shift:max_shift)
    for (dy in -max_shift:max_shift) {
      shifted <- matrix(NA, rows, cols)
      r_dst <- (1:rows) + dy; c_dst <- (1:cols) + dx
      ok_r <- r_dst >= 1 & r_dst <= rows
      ok_c <- c_dst >= 1 & c_dst <= cols

      shifted[r_dst[ok_r], c_dst[ok_c]] <-
        target_mat[(1:rows)[ok_r], (1:cols)[ok_c]]

      overlap <- !is.na(ref_mat) & !is.na(shifted)

      if (sum(overlap) > 1000) {
        cc <- cor(ref_mat[overlap], shifted[overlap], use = "complete.obs")
        if (!is.na(cc) && cc > best_cor) {
          best_cor <- cc; best_dx <- dx; best_dy <- dy
        }
      }
    }

  list(dx = best_dx, dy = best_dy, cor = best_cor)
}
```

---

## Apply the Optimal Shift

```r
apply_shift <- function(mat, dx, dy) {
  rows <- nrow(mat); cols <- ncol(mat)
  shifted <- matrix(NA, rows, cols)

  r_dst <- (1:rows) + dy; c_dst <- (1:cols) + dx
  ok_r <- r_dst >= 1 & r_dst <= rows
  ok_c <- c_dst >= 1 & c_dst <= cols

  shifted[r_dst[ok_r], c_dst[ok_c]] <-
    mat[(1:rows)[ok_r], (1:cols)[ok_c]]

  shifted
}
```

---

## Coregister Bands

The **Red** band is used as the reference. Green and NIR are aligned to it.

```r
ref_red <- band_matrices$Red
aligned_matrices <- list(Red = ref_red)

for (b in c("Green", "NIR")) {
  res <- align_band_simple(ref_red, band_matrices[[b]])
  aligned_matrices[[b]] <- apply_shift(band_matrices[[b]], res$dx, res$dy)
}
```

---

## Rebuild Raster Stack and Visualize Results

```r
grea <- rast(aligned_matrices$Green)
reda <- rast(aligned_matrices$Red)
nira <- rast(aligned_matrices$NIR)

after <- c(grea, reda, nira)
plotRGB(after, r = 3, g = 2, b = 1, stretch = "lin")
```

### Before vs After Comparison

```r
par(mfrow = c(1, 2))
plotRGB(ms, r = 3, g = 2, b = 1, stretch = "lin")
plotRGB(after, r = 3, g = 2, b = 1, stretch = "lin")
```

---

## Notes and Limitations

* This is a **rigid translation-only** alignment (no rotation or scaling).
* Best suited for **small misalignments** (≤ ~20 pixels).
* Correlation assumes similar radiometric structure between bands.
* For larger distortions or higher accuracy, consider feature-based or external tools (e.g. OpenCV, GDAL, or photogrammetry software).

---

## Use Cases

* Multispectral drone imagery
* NDVI preprocessing
* Band alignment prior to classification or texture analysis

---

## License

Use freely, modify boldly, and cite politely 😄
