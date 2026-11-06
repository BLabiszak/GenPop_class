
plot_genomic_vulnerability <- function(coords_sf, gv_values, title = "Genomic Vulnerability", n_breaks = 8) {
  library(sf)
  library(dplyr)
  library(gstat)
  library(raster)
  library(rnaturalearth)
  library(ggplot2)
  library(RColorBrewer)
  library(ggpubr)
  
  # ---- 1. Add vulnerability values ----
  coords_sf$gv <- round(gv_values, 3)
  world <- ne_countries(scale = "medium", returnclass = "sf")
  transect <- world[world$admin %in% c("Poland", "Finland", "Estonia", "Latvia", "Lithuania"),]
  
  # ---- 2. Define bounding box (PL + Baltic + FI) ----
  bbox_vals <- c(xmin = 14.13467, xmax = 31.51066, ymin = 48.01889, ymax = 70.07551)
  bbox <- st_bbox(bbox_vals, crs = 4326)
  transect_sf <- st_as_sfc(bbox)
  bbox_sp <- as(transect_sf, "Spatial")
  
  # ---- 3. Filter coords to bbox ----
  coords_sf <- st_filter(coords_sf, transect_sf)
  
  # ---- 4. IDW interpolation over bbox grid ----
  coords_sp <- as(coords_sf, "Spatial")
  
  grd <- as.data.frame(spsample(bbox_sp, type = "regular", n = 50000))
  names(grd) <- c("x", "y")
  coordinates(grd) <- ~x + y
  gridded(grd) <- TRUE
  proj4string(grd) <- CRS("+proj=longlat +datum=WGS84")
  
  idw_result <- idw(gv ~ 1, coords_sp, newdata = grd, idp = 2.0)
  r <- raster(idw_result)
  crs(r) <- CRS("+proj=longlat +datum=WGS84")
  
  r_cropped <- crop(r, transect)
  r_masked <- mask(r_cropped, transect)
  
  # ---- 5. Convert raster to data frame ----
  r_df <- as.data.frame(r_masked, xy = TRUE)
  colnames(r_df)[3] <- "gv"
  
  # ---- 6. Dynamic breaks based on gv range ----
  gv_range <- range(r_df$gv, na.rm = TRUE)
  breaks <- seq(gv_range[1], gv_range[2], length.out = n_breaks + 1)
  
  # ---- 7. Crop world basemap to bbox ----
  world_sf <- ne_countries(scale = "medium", returnclass = "sf")
  world_crop <- st_crop(world_sf, bbox)
  
  # ---- 8. Extract lon/lat for plotting points ----
  coords_sf <- coords_sf %>%
    mutate(lon = st_coordinates(.)[, 1],
           lat = st_coordinates(.)[, 2])
  
  # ---- 9. Final ggplot ----
  ggplot() + 
    geom_raster(data = r_df, aes(x = x, y = y, fill = cut(gv, breaks = breaks, include.lowest = TRUE))) +
    geom_sf(data = world_crop, fill = NA, color = "black", linewidth = 0.3) +
    geom_point(data = coords_sf, aes(x = lon, y = lat), size = 2, color = "black") +
    scale_fill_brewer(palette = "YlOrRd", name = "Genomic\nVulnerability", direction = 1) +
    coord_sf(xlim = bbox_vals[c("xmin", "xmax")], ylim = bbox_vals[c("ymin", "ymax")], expand = FALSE) +
    labs(x = "Longitude", y = "Latitude", title = title) +
    theme_pubr(legend = "right")
}
#plot_genomic_vulnerability(coords_sf = coords_sf,gv_values = gv_df_RCP26$gv,n_breaks = 8)



  