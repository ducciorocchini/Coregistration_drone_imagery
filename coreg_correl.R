library(terra)
setwd("~/Desktop/")
# gre = rast("DJI_20251115110614_0003_MS_G.tiff")
# red = rast("DJI_20251115110614_0003_MS_R.tiff")
# nir = rast("DJI_20251115110614_0003_MS_NIR.tiff")

gre = rast("ms_crop_gre.tif")
red = rast("ms_crop_red.tif")
nir = rast("ms_crop_nir.tif")

ms = c(gre, red, nir)
plotRGB(ms, r=3, g=2, b=1, stretch="lin")

# other sets
#library(terra)
#setwd("/Users/ducciorocchini/Desktop/DJI_202511151151_002")
#gre = rast("DJI_20251128093951_0003_MS_G.TIF")
#red = rast("DJI_20251128093951_0003_MS_R.TIF")
#nir = rast("DJI_20251128093951_0003_MS_NIR.TIF")



band_names <- c("Green", "Red", "NIR")

# Convert raster layers to matrices
band_matrices <- lapply(1:nlyr(ms), function(i) {
  as.matrix(ms[[i]], wide = TRUE)
  })
names(band_matrices) <- band_names

align_band_simple <- function(ref_mat, target_mat, max_shift=20) {
  best_cor <- -Inf; best_dx <- 0; best_dy <- 0
  rows <- nrow(ref_mat); cols <- ncol(ref_mat)
  for(dx in -max_shift:max_shift) for(dy in -max_shift:max_shift) {
    shifted <- matrix(NA, rows, cols)
    r_dst <- (1:rows)+dy; c_dst <- (1:cols)+dx
    ok_r <- r_dst >= 1 & r_dst <= rows; ok_c <- c_dst >= 1 & c_dst <= cols
    shifted[r_dst[ok_r], c_dst[ok_c]] <- target_mat[(1:rows)[ok_r], (1:cols)[ok_c]]
    overlap <- !is.na(ref_mat) & !is.na(shifted)
    if(sum(overlap)>1000) {
      cc <- cor(ref_mat[overlap], shifted[overlap], use="complete.obs")
      if(!is.na(cc) && cc>best_cor) { best_cor <- cc; best_dx <- dx; best_dy <- dy }
    }
  }
  list(dx=best_dx, dy=best_dy, cor=best_cor)
}

apply_shift <- function(mat, dx, dy) {
  rows <- nrow(mat); cols <- ncol(mat); shifted <- matrix(NA, rows, cols)
  r_dst <- (1:rows)+dy; c_dst <- (1:cols)+dx
  ok_r <- r_dst >=1 & r_dst<=rows; ok_c <- c_dst>=1 & c_dst<=cols
  shifted[r_dst[ok_r], c_dst[ok_c]] <- mat[(1:rows)[ok_r], (1:cols)[ok_c]]
  shifted
}

ref_red <- band_matrices$Red
aligned_matrices <- list(Red=ref_red)
for(b in c("Green","NIR")) {
  res <- align_band_simple(ref_red, band_matrices[[b]])
  aligned_matrices[[b]] <- apply_shift(band_matrices[[b]], res$dx, res$dy)
}
ref_red
reda <- rast(aligned_matrices$Red)
grea <- rast(aligned_matrices$Green)
nira <- rast(aligned_matrices$NIR)
after <- c(grea, reda, nira)
plotRGB(after, r=3, g=2, b=1, stretch="lin")  

par(mfrow=c(1,2))
plotRGB(ms, r=3, g=2, b=1, stretch="lin")
plotRGB(after, r=3, g=2, b=1, stretch="lin")  
