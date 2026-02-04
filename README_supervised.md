# Drone Multispectral Image Co-registration in R

**Author:** Sofiya Tumanova  
**Data:** Multispectral TIFFs (Red, Green, NIR)  
**Requires:** R >= 4.2, `terra` package  

---

## Overview

Interactive co-registration of drone multispectral bands:

- Select the same Ground Control Point (GCP) in Red, Green, and NIR bands.
- Align Green and NIR to Red using pixel shifts.
- Crop to common overlap and save stacked raster.
- Compare before vs after co-registration in false-color and R+G / R+NIR plots.

---

## Usage

1. Download the contents of the `data/` folder.  
2. Open `clean_coregistration_script_github.R` in R or RStudio.  
3. Run the script and click on the GCP in each band.  
4. Aligned raster is saved in `output/aligned_stack.tif`.  

---

## Example Plots

### 3-band False Color
![Before vs After](result_images/before%20vs.%20after%20coregistration.jpeg)

### Red + Green
![R+G Before vs After](result_images/red+green_before_after.jpeg)

### Red + NIR
![R+NIR Before vs After](result_images/red+nir_before_after.jpeg)


