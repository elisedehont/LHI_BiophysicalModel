install.packages("raster")
install.packages("sf")
install.packages("sp")
install.packages("ggplot2")
install.packages("ncdf4")
install.packages("dplyr")
install.packages("transformr")
install.packages("gifski")
install.packages(c("terra", "spatstat"))
install.packages("akima")
install.packages(c("rnaturalearth", "rnaturalearthdata"))
install.packages("geosphere")
install.packages("ggforce")
install.packages("ggspatial")
install.packages("geosphere")
install.packages("gdistance")
install.packages("reticulate")
install.packages("reshape2")
install.packages("pheatmap")
install.packages("patchwork")
install.packages("gstat")
install.packages("FNN")
install.packages("writexl")
install.packages("ggalluvial")
install.packages("readxl")
library(dplyr)
library(ggplot2)
library(gganimate)
library(transformr)
library(gifski)
library(sf)
library(raster)
library(ggplot2)
##csv to tif#
library(spatstat) 
library(akima)
###transform .nc data into velocity rasters#######
library(ncdf4)   # To read NetCDF files
library(raster)  # To manipulate and visualize raster data
library(terra)   # Alternative for raster operations
library(lubridate)
library(sp)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(geosphere)
library(ggforce)
library(geosphere)
library(gdistance)
library(reticulate)
library(reshape2)
library(pheatmap)
library(patchwork)
library(gstat)
library(FNN)
library(ggalluvial)
library(writexl)
library(readxl)
#open file
nc_file <- "C:/Users/33778/Desktop/StageM2/currents_149_168true.nc"
nc <- nc_open(nc_file)

# Print variable names and dimensions to check
print(nc$var)
print(nc$dim)
###U (eastward velocity, uo) and V (northward velocity, vo)##
uo <- ncvar_get(nc, "uo")
vo <- ncvar_get(nc, "vo")

#extract time values
time_values <- ncvar_get(nc, "time")  # Check if 'time' is the correct name
time_units <- ncatt_get(nc, "time", "units")$value  # Get unit metadata
print(time_units)
#convert time values to dates
time_origin <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
time_dates <- time_origin + time_values
time_dates <- time_dates[format(time_dates, "%m") %in% c("10", "11", "12", "01", "02")]
unique(format(time_dates, "%Y-%m"))
table(format(time_dates, "%m"))#oct to february  OK


#print first and last dates
print(range(time_dates)) ##10-2011 to 02 2021 OK
diff_time <- diff(time_dates)  # Compute time differences
print(unique(diff_time))  # Print unique time intervals
print(head(time_dates, 20))  #first 20 time points OK
print(tail(time_dates, 20))  #last 20 time points OK

#Convert dates to year-month format
time_months <- format(time_dates, "%Y-%m")

#Create an empty list to store monthly averages
monthly_avg_list <- list()
#longitude and latitude
longitude <- ncvar_get(nc, "longitude")
latitude <- ncvar_get(nc, "latitude")
#empty list to store annual rasters
annual_raster_list <- list()
#Loop through each month
for (month in unique(time_months[grep("-(10|11|12|01|02)$", time_months)])) {
  print(paste("Processing month:", month))  
  
#data indices for the month
  indices <- which(time_months == month)
  
#extraction of monthly velocity data (U and V) 
  print(dim(uo))  #should be lon, lat, time or lon, lat, depth, time
  print(dim(vo))
  if (length(dim(uo)) == 4) {
    uo <- uo[,,1,]  
    vo <- vo[,,1,]
  }
  uo_month <- apply(uo[,,indices], c(1,2), mean, na.rm = TRUE)
  vo_month <- apply(vo[,,indices], c(1,2), mean, na.rm = TRUE)
  
  
  
#convert to rasters
  uo_month_raster <- raster(t(uo_month))
  vo_month_raster <- raster(t(vo_month))
  uo_month_raster <- flip(uo_month_raster, direction = "y")
  vo_month_raster <- flip(vo_month_raster, direction = "y")

  
#Defining the raster extent
  extent(uo_month_raster) <- extent(min(longitude), max(longitude), min(latitude), max(latitude))
  extent(vo_month_raster) <- extent(min(longitude), max(longitude), min(latitude), max(latitude))
  #compute velocity and direction
  velocity_month <- sqrt(uo_month_raster^2 + vo_month_raster^2)
  direction_month <- atan2(vo_month_raster, uo_month_raster) * (180 / pi)
  
#results storage
  monthly_avg_list[[month]] <- list(velocity = velocity_month, direction = direction_month)
  
#saving rasters for each month
  writeRaster(velocity_month, paste0("velocity_", month, ".tif"), format = "GTiff", overwrite = TRUE)
  writeRaster(direction_month, paste0("direction_", month, ".tif"), format = "GTiff", overwrite = TRUE)
}

#inspection OK
nc_close(nc)
###rasters OK to work with


#load currents direction and velocity rasters for a month
#define the month and year
month_year <- "2014-12"
#rasters download
SCD <- rast("C:/Users/33778/Desktop/StageM2/velocity/velocity_2014-12.tif")
SCV <- rast("C:/Users/33778/Desktop/StageM2/direction/direction_2014-12.tif")

#Define the area of interest (AOI)
AOI <- ext(149, 168, -32, -17)
world <- ne_countries(scale = "medium", returnclass = "sf")
land <- st_crop(world, c(xmin = 149, xmax = 168, ymin = -32, ymax = -17))

#Cut the rasters
SCD_crop <- crop(SCD, AOI)
SCV_crop <- crop(SCV, AOI)

#convert into data.frame
vel_df <- as.data.frame(SCD_crop, xy = TRUE)
dir_df <- as.data.frame(SCV_crop, xy = TRUE)

#Naming colums
colnames(vel_df) <- c("x", "y", "velocity")
colnames(dir_df) <- c("x", "y", "direction")

# Merge the two data.frames and delete the NAs
current_df <- na.omit(merge(vel_df, dir_df, by = c("x", "y")))

#function for calculating arrow coordinates
deg2rad <- function(deg) (pi * deg) / 180
arrowCOORD <- function(x, y, r, s = 0.08) {
  ENDx <- x + sin(deg2rad(r)) * s
  ENDy <- y + cos(deg2rad(r)) * s
  STAx <- x - sin(deg2rad(r)) * s
  STAy <- y - cos(deg2rad(r)) * s
  ar1X <- ENDx + sin(deg2rad(r + 135)) * (s / 2)
  ar1Y <- ENDy + cos(deg2rad(r + 135)) * (s / 2)
  ar2X <- ENDx + sin(deg2rad(r - 135)) * (s / 2)
  ar2Y <- ENDy + cos(deg2rad(r - 135)) * (s / 2)
  return(list('ENDx' = ENDx, 'ENDy' = ENDy, 'STAx' = STAx, 'STAy' = STAy,
              'ar1X' = ar1X, 'ar1Y' = ar1Y, 'ar2X' = ar2X, 'ar2Y' = ar2Y))
}

#reducing the arrows density
step <- 5
current_df_sampled <- current_df[seq(1, nrow(current_df), by = step), ]
##add Lord Howe island-can be done with rnaturalearthdata later
lord_howe <- data.frame(
  lon = 159.0833,  # Longitude Lord Howe Island
  lat = -31.5500   # Latitude Lord Howe Island
)

#plot
ggplot() +
  geom_sf(data = land, fill = "darkgreen", color = "black") +
  geom_segment(data = current_df_sampled, 
               aes(x = x, y = y,
                   xend = x + sin(deg2rad(direction)) * velocity * 0.1,
                   yend = y + cos(deg2rad(direction)) * velocity * 0.1,
                   color = velocity),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_point(data = lord_howe, aes(x = lon, y = lat), color = "red", size = 3) +
  geom_text(data = lord_howe, aes(x = lon, y = lat, label = "Lord Howe Island"), 
            hjust = 0, vjust = -1, color = "black", size = 3) +
  scale_color_gradient(low = "lightblue", high = "purple", 
                       name = "Velocity") +
  labs(title = paste("Current Vectors -", month_year),
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  coord_sf(xlim = c(149, 168), ylim = c(-32, -17))

######currents characterisation####


#Lord Howe Island
lord_howe <- data.frame(
  lon = 159.0833,
  lat = -31.5500
)

#buffer around LHI (~8°)
buffer_radius <- 8  #degrees
buffer_extent <- ext(
  lord_howe$lon - buffer_radius, lord_howe$lon + buffer_radius,
  lord_howe$lat - buffer_radius, lord_howe$lat + buffer_radius
)

#Data extraction in the area of interest (buffer)
SCD_buffer <- crop(SCD_crop, buffer_extent)  # Vitesse
SCV_buffer <- crop(SCV_crop, buffer_extent)  # Direction


current_df <- current_df[current_df$month %in% c(10, 11, 12, 1, 2), ]
current_df <- current_df %>% filter(!is.na(direction) & !is.na(velocity))
time_values <- as.Date(paste0(current_df$year, "-", current_df$month, "-01"), format = "%Y-%m-%d")


#dataframe
vel_df <- as.data.frame(SCD_buffer, xy = TRUE)
dir_df <- as.data.frame(SCV_buffer, xy = TRUE)
colnames(vel_df) <- c("x", "y", "velocity")
colnames(dir_df) <- c("x", "y", "direction")

#merge data and delete NAs
current_df <- na.omit(merge(vel_df, dir_df, by = c("x", "y")))

#calculating the direction of points towards Lord Howe
bearing_to_lord_howe <- function(lon, lat) {
  bearing(c(lon, lat), c(lord_howe$lon, lord_howe$lat))  
}

current_df$bearing_to_lh <- mapply(bearing_to_lord_howe, current_df$x, current_df$y)

#determining if the current is flowing into or out of Lord Howe
angle_diff <- abs(current_df$direction - current_df$bearing_to_lh)  #difference between current direction and direction towards LHI
current_df$type <- ifelse(angle_diff < 90, "Arriving", "Leaving")  #<90°= incoming current, >90° = outgoing

#statistics
current_stats <- current_df %>%
  group_by(type) %>%
  summarise(
    mean_velocity = mean(velocity, na.rm = TRUE),
    sd_velocity = sd(velocity, na.rm = TRUE),
    mean_direction = mean(direction, na.rm = TRUE)
  )

#results
print(current_stats)
#1 Arriving         mean velo :0.252       sd dir : 0.106         mean dir: -46.6
#2 Leaving          mean velo: 0.305       sd vel : 0.132         mean dir : -37.9

#current vectors map around LHI
ggplot() +
  geom_segment(data = current_df, 
               aes(x = x, y = y,
                   xend = x + sin(deg2rad(direction)) * velocity * 0.3,
                   yend = y + cos(deg2rad(direction)) * velocity * 0.3,
                   color = velocity),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_point(data = lord_howe, aes(x = lon, y = lat), color = "red", size = 3) +
  geom_text(data = lord_howe, aes(x = lon, y = lat, label = "Lord Howe Island"), 
            hjust = 0, vjust = -1, color = "black", size = 4) +
  scale_color_gradient(low = "lightblue", high = "purple", name = "Velocity") +
  labs(title = "Current Vectors (Sep-Jan) Around Lord Howe", x = "Longitude", y = "Latitude") +
  theme_minimal()

#incoming/outgoing currents map
ggplot() +
  geom_point(data = current_df, aes(x = x, y = y, color = type), alpha = 0.7) +
  geom_point(data = lord_howe, aes(x = lon, y = lat), color = "red", size = 3) +
  scale_color_manual(values = c("Arriving" = "blue", "Leaving" = "orange")) +
  labs(title = "Arriving vs Leaving Currents (Sep-Jan)", x = "Longitude", y = "Latitude", color = "Current Type") +
  theme_minimal()

#histogram of currents speed
ggplot(current_df, aes(x = velocity, fill = type)) +
  geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
  scale_fill_manual(values = c("Arriving" = "blue", "Leaving" = "orange")) +
  labs(title = "Velocity Distribution (Sep-Jan)", 
       x = "Velocity (m/s)", y = "Count", fill = "Current Type") +
  theme_minimal()

####zoom

#world data
world <- ne_countries(scale = "large", returnclass = "sf")

#coastlines data
coastline <- ne_coastline(scale = "large", returnclass = "sf")

#coordinates for Lord Howe Island
zoom_lon <- c(158.0, 160.0)  # Longitudes
zoom_lat <- c(-32.0, -31.0)  # Latitudes

#currents map and LHI design
ggplot() +
  # Fond de carte des côtes
  geom_sf(data = world, fill = "gray80", color = "black") +
  geom_sf(data = coastline, color = "black") +  
  
#current vectors
  geom_segment(data = current_df, 
               aes(x = x, y = y,
                   xend = x + sin(deg2rad(direction)) * velocity * 0.3,
                   yend = y + cos(deg2rad(direction)) * velocity * 0.3,
                   color = velocity),
               arrow = arrow(length = unit(0.3, "cm")), alpha = 0.7) +
  
#LHI
  geom_point(data = lord_howe, aes(x = lon, y = lat), size = 3) +
  geom_text(data = lord_howe, aes(x = lon, y = lat, label = "Lord Howe Island"), 
            hjust = 0, vjust = -1, color = "black", size = 4) +
  
#velocity colors
  scale_color_gradient(low = "lightblue", high = "purple", name = "Velocity") +
  
#zoom limits
  coord_sf(xlim = zoom_lon, ylim = zoom_lat, expand = FALSE) +
  
#labs
  labs(title = "Current Vectors (Sep-Jan) Around Lord Howe", 
       x = "Longitude", y = "Latitude") +
  theme_minimal()

#####unzoom

#coordinates
buffer_extent <- ext(
  149, 168,  # Longitude
  -32, -17   # Latitude
)

#data extraction
SCD_buffer <- crop(SCD_crop, buffer_extent) 
SCV_buffer <- crop(SCV_crop, buffer_extent)  

#dataframe
vel_df <- as.data.frame(SCD_buffer, xy = TRUE)
dir_df <- as.data.frame(SCV_buffer, xy = TRUE)
colnames(vel_df) <- c("x", "y", "velocity")
colnames(dir_df) <- c("x", "y", "direction")

#merge and delete NAs
current_df <- na.omit(merge(vel_df, dir_df, by = c("x", "y")))

#coastlines
coastline <- ne_countries(scale = "medium", returnclass = "sf")
australia_coast <- coastline %>% filter(admin == "Australia")

#current vectors
ggplot() +
  geom_sf(data = australia_coast, fill = "gray80", color = "black") + 
  geom_segment(data = current_df, 
               aes(x = x, y = y,
                   xend = x + sin(deg2rad(direction)) * velocity * 0.3,
                   yend = y + cos(deg2rad(direction)) * velocity * 0.3,
                   color = velocity),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_point(data = lord_howe, aes(x = lon, y = lat), color = "red", size = 1) +
  geom_text(data = lord_howe, aes(x = lon, y = lat, label = "Lord Howe Island"), 
            hjust = 0, vjust = -1, color = "black", size = 4) +
  scale_color_gradient(low = "lightblue", high = "purple", name = "Velocity") +
  coord_sf(xlim = c(149, 168), ylim = c(-32, -17), expand = FALSE) +  
  labs(title = "Currents between East Australia & Lord Howe", x = "Longitude", y = "Latitude") +
  theme_minimal()

#incoming/outgoing current calcul
bearing_to_lord_howe <- function(lon, lat) {
  bearing(c(lon, lat), c(lord_howe$lon, lord_howe$lat))  
}

current_df$bearing_to_lh <- mapply(bearing_to_lord_howe, current_df$x, current_df$y)


angle_diff <- abs(current_df$direction - current_df$bearing_to_lh)

#incoming/outgoing
current_df$type <- ifelse(angle_diff < 90, "Arriving", "Leaving")
ggplot() +
  geom_point(data = current_df, aes(x = x, y = y, color = type), alpha = 0.7) +
  geom_point(data = lord_howe, aes(x = lon, y = lat), color = "red", size = 3) +
  scale_color_manual(values = c("Arriving" = "blue", "Leaving" = "orange")) +
  labs(title = "Arriving vs Leaving Currents (Sep-Jan)", x = "Longitude", y = "Latitude", color = "Current Type") +
  theme_minimal()

#histogram of currents speed
ggplot(current_df, aes(x = velocity, fill = type)) +
  geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
  scale_fill_manual(values = c("Arriving" = "blue", "Leaving" = "orange")) +
  labs(title = "Velocity Distribution (Sep-Jan)", 
       x = "Velocity (m/s)", y = "Count", fill = "Current Type") +
  theme_minimal()



#####connectivity matrice : geographical distance
#file
reef_sites <- read.csv("C:/Users/33778/Desktop/StageM2/Elise_Connectivity_ReefSites.csv", sep = ";")
#check
head(reef_sites)


#convert data to spatial object
reef_sf <- st_as_sf(reef_sites, coords = c("Lon", "Lat"), crs = 4326)

#world
world <- ne_countries(scale = "medium", returnclass = "sf")

#map with reefs
ggplot() +
  geom_sf(data = world, fill = "antiquewhite") +
  geom_sf(data = reef_sf, color = "red", size = 2) +
  coord_sf(xlim = c(149, 168), ylim = c(-32, -17), expand = FALSE) +  
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", which_north = "true") +
  theme_minimal()

#calculate the geodesic distances between all the points
#extract data
coords <- reef_sites[, c("Lon", "Lat")]

#calculate connectivity matrix (meters)
distance_matrix <- distm(coords, fun = distHaversine)

#convert to kilometers
distance_matrix_km <- distance_matrix / 1000

#add reef names
rownames(distance_matrix_km) <- reef_sites$Reef_Name
colnames(distance_matrix_km) <- reef_sites$Reef_Name

#see
print(round(distance_matrix_km, 2))
#see and save matrix
round(distance_matrix_km, 2)
write.csv(distance_matrix_km, "C:/Users/33778/Desktop/Matrice_Connectivite_Recifs.csv")

#####realistic connectivity matrix based on current speed (each cell represents the cost of travelling between two reefs via currents) : 

#defining the velocity file
file_velocity <- "C:/Users/33778/Desktop/StageM2/velocity/"

#listing velocity .tif files
file_list <- list.files(path = file_velocity, pattern = "\\.tif$", full.names = TRUE)

#load the rasters in a stack
stack_velocity <- rast(file_list)
#calculate mean velocity pixel by pixel
mean_velocity <- mean(stack_velocity, na.rm = TRUE)

#plot and save
plot(mean_velocity, main = "Mean velocity (sept-dec, 2011-2021)")
writeRaster(mean_velocity, "C:/Users/33778/Desktop/StageM2/MeanVelocity.tif", overwrite = TRUE)

#defining the direction file
file_direction <- "C:/Users/33778/Desktop/StageM2/direction/"

#listing direction .tif files
file_list2 <- list.files(path = file_direction, pattern = "\\.tif$", full.names = TRUE)

#load the rasters in a stack
stack_direction <- rast(file_list2)
#calculate mean velocity pixel by pixel
mean_direction <- mean(stack_direction, na.rm = TRUE)

#plot and save
plot(mean_direction, main = "Mean direction (sept-dec, 2011-2021)")
writeRaster(mean_direction, "C:/Users/33778/Desktop/StageM2/MeanDirection.tif", overwrite = TRUE)



###anisotropic friction
#load rasters
raster_velocity <- rast("C:/Users/33778/Desktop/StageM2/MeanVelocity.tif")
raster_direction <- rast("C:/Users/33778/Desktop/StageM2/MeanDirection.tif")
#create U and V rasters from directional components
raster_U <- raster_velocity * sin(raster_direction)
raster_V <- raster_velocity * cos(raster_direction)
#calculate anisotropic friction
raster_friction_directional <- 1 / sqrt(raster_U^2 + raster_V^2)

# Replace infinite values (0 speed) with NA
raster_friction_directional[raster_friction_directional == Inf] <- NA
#see and download
plot(raster_friction_directional, main="Friction Directionnelle")
writeRaster(raster_friction_directional, "C:/Users/33778/Desktop/StageM2/raster_friction_directional.tif", overwrite = TRUE)

##we can see that the resolution is too low so interpolation in SAGA
raster_friction_smoothed <- raster("C:/Users/33778/Desktop/StageM2/friction_smoothed.tif")
# Visualiser le résultat
plot(raster_friction_smoothed, main="smoothed friction raster")
#replacing NaN values because Circuitscape cant read it
raster_friction_smoothed[is.na(raster_friction_smoothed)] <- -9999

#export to .asc format because Circuitscape cant read .tif
writeRaster(raster_friction_smoothed, "C:/Users/33778/Desktop/StageM2/friction_smoothed.asc", NAflag=-9999, overwrite=TRUE)
#asc format for circuitscape

###circuitscape raster pairwise

resistance_matrix <- read.table("C:/Users/33778/Desktop/StageM2/connectivity_matrix_csc_resistances", quote="\"", comment.char="")
reef_sites <- read.csv("C:/Users/33778/Desktop/StageM2/reefs_format.csv", sep = ";")
reef_names <- reef_sites$id
resistance_matrix_clean <- resistance_matrix[-1, -1]
rownames(resistance_matrix_clean) <- reef_names
colnames(resistance_matrix_clean) <- reef_names
max_resistance <- max(resistance_matrix_clean, na.rm = TRUE)
#standardisation to balance differences in resistance
#and avoid extreme effects while maintaining proportionality
resistance_matrix_normalized <- resistance_matrix_clean / max_resistance
connectivity_matrix <- exp(-resistance_matrix_normalized)
write.csv(connectivity_matrix, "C:/Users/33778/Desktop/StageM2/Connectivity_matrix.csv")

library(pheatmap)
pheatmap(connectivity_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "blue", "darkblue"))(100),
         border_color = "grey60")
#illegible beacuse 439 lines

###pca
pca_res <- prcomp(connectivity_matrix, scale. = TRUE)

#results from the 2 first principal components
pca_df <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  Reefs = rownames(connectivity_matrix)  # si tu as nommé tes lignes avec les récifs
)
ggplot(pca_df, aes(x = PC1, y = PC2, label = Reefs)) +
  geom_point(color = "darkgreen", size = 2, alpha = 0.7) +
  geom_text(size = 3, vjust = -1, check_overlap = TRUE) +
  theme_minimal() +
  labs(title = "PCA sur la matrice de connectivité",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "%)"))

#complex heat map
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
Heatmap(as.matrix(connectivity_matrix),
        name = "Connectivité",
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = colorRampPalette(c("white", "purple", "orange", "yellow"))(100))

###matrix analysis

importance_import <- colSums(connectivity_matrix)
importance_export <- rowSums(connectivity_matrix)
top_sources <- sort(importance_export, decreasing = TRUE)[1:30]
print("Top 30 reef sources :")
print(top_sources)

top_sinks <- sort(importance_import, decreasing = TRUE)[1:30]
print("Top 30 reef sinks :")
print(top_sinks)
#values very similar, similarity test :
#distribution of connectivity weights
#vectorise the matrix without the diagonal
connectivity_values <- connectivity_matrix[lower.tri(connectivity_matrix, diag = FALSE)]
#calculate the variation coefficient
cv_connectivity <- sd(connectivity_values) / mean(connectivity_values)
cv_connectivity
cv_percent <- cv_connectivity * 100
paste("variation coefficient (%) :", round(cv_percent, 2))


###simulation
set.seed(42)  #reproductibility
#Nafis, N. (2023, 24 janvier). The story behind ‘random.seed(42)’ in machine learning. Medium. https://medium.com/geekculture/the-story-behind-random-seed-42-in-machine-learning-b838c4ac290a
larval_production <- runif(nrow(connectivity_matrix), min = 0.5, max = 5)
recruitment_success <- runif(ncol(connectivity_matrix), min = 0.0001, max = 0.01)
connectivity_matrix_bio <- connectivity_matrix

#apply production factor to each line (emission)
connectivity_matrix_bio <- sweep(connectivity_matrix_bio, 1, larval_production, FUN = "*")

#apply installation success to each column (reception)
connectivity_matrix_bio <- sweep(connectivity_matrix_bio, 2, recruitment_success, FUN = "*")
importance_export_bio <- rowSums(connectivity_matrix_bio)
importance_import_bio <- colSums(connectivity_matrix_bio)
top_sources <- names(sort(importance_export_bio, decreasing = TRUE))[1:30]
top_sinks   <- names(sort(importance_import_bio, decreasing = TRUE))[1:30]
#quick visualization
barplot(sort(importance_export_bio, decreasing = TRUE)[1:30], main="Top 10 reef sources", col="orange")
barplot(sort(importance_import_bio, decreasing = TRUE)[1:30], main="Top 10 sink reefs", col="blue")

###quadrant plot
#dataframe
reef_df <- data.frame(
  reef = colnames(connectivity_matrix_bio),
  export = rowSums(connectivity_matrix_bio),
  import = colSums(connectivity_matrix_bio)
)

#calculate metrics
reef_df$net_flow <- reef_df$export - reef_df$import
reef_df$inward_degree <- reef_df$import  

#transformation into percentiles
reef_df$net_flow_percentile <- ecdf(reef_df$net_flow)(reef_df$net_flow) * 100
reef_df$inward_degree_percentile <- ecdf(reef_df$inward_degree)(reef_df$inward_degree) * 100

#ecological classification
reef_df$role <- "Other"
reef_df$role[reef_df$net_flow_percentile >= 90] <- "Source"
reef_df$role[reef_df$net_flow_percentile <= 10] <- "Sink"
reef_df$role[reef_df$inward_degree_percentile >= 90] <- "Connectivity Hub"

#add colours
colors <- c("Source" = "#003399", "Sink" = "#FF3399", "Connectivity Hub" = "#66CC99", "Other" = "grey80")

#plot
p1 <- ggplot(reef_df, aes(x = inward_degree_percentile, y = net_flow_percentile)) +
  geom_point(aes(color = role), shape = 21, size = 2, stroke = 0.6, fill = "white") +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = 90, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 90, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "black") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Reef Connectivity Roles",
    subtitle = "Based on Net Larval Flow and Inward Connectivity Percentiles",
    x = "Inward Connections (percentile)",
    y = "Net Larval Flow (percentile)",
    color = "Reef Role"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  )

#barplot on the right
role_counts <- as.data.frame(table(reef_df$role))
colnames(role_counts) <- c("Role", "Count")

p2 <- ggplot(role_counts, aes(x = reorder(Role, -Count), y = Count, fill = Role)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = colors) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4) +
  labs(
    title = "Number of Reefs per Role",
    x = NULL, y = "Number of Reefs"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

#gather the 2 plots with patwork package
p1 + p2 + plot_layout(widths = c(2.5, 1))

###map
#merge roles with their reef names
reef_map_df <- merge(reef_sites, reef_df[, c("reef", "role")], by.x = "id", by.y = "reef", all.x = TRUE)

#NA data becomes "other"
reef_map_df$role[is.na(reef_map_df$role)] <- "Other"

#define colors
colors <- c("Source" = "#003399", "Sink" = "#FF3399", "Connectivity Hub" = "#66CC99", "Other" = "grey80")

#australia map
australia_LHI <- ne_countries(scale = "medium", country = "Australia", returnclass = "sf")

#filter key reefs (no "other")
reef_map_df_filtered <- subset(reef_map_df, role != "Other")

#map with key reefs (no "other" if applying data=reef_map_df_filtered)
#sGBR
ggplot() +
  geom_sf(data = australia_LHI, fill = "grey90", color = "black") +
  geom_point(data = reef_map_df, aes(x = lon, y = lat, color = role), size = 3, alpha = 0.9) +
  scale_color_manual(values = colors) +
  coord_sf(xlim = c(149, 154), ylim = c(-24, -19), expand = FALSE) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Reef Connectivity Roles - Southern Great Barrier Reef",
    subtitle = "Sources, Sinks, Connectivity Hubs and other reefs",
    x = "Longitude", y = "Latitude", color = "Reef Role"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  )

#LHI
land <- ne_download(scale = "large", type = "land", category = "physical", returnclass = "sf")
lord_howe <- land[st_intersects(land, st_sfc(st_point(c(159.08, -31.55)), crs = 4326), sparse = FALSE), ]

ggplot() +
  geom_sf(data = lord_howe, fill = "grey90", color = "black") +  # Affichage de Lord Howe Island
  geom_point(data = reef_map_df, aes(x = lon, y = lat, color = role), size = 3, alpha = 0.9) +
  scale_color_manual(values = colors) +
  coord_sf(xlim = c(159, 159.2), ylim = c(-31.6, -31.5), expand = FALSE) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Reef Connectivity Roles - Lord Howe Island",
    subtitle = "Sources, Sinks, Connectivity Hubs and other reefs",
    x = "Longitude", y = "Latitude", color = "Reef Role"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  )
#all additional reefs
world <- ne_countries(scale = "large", returnclass = "sf")

australia_nc <- world[world$name %in% c("Australia", "New Caledonia"), ]

ggplot() +
  geom_sf(data = australia_nc, fill = "grey90", color = "black") +
  geom_point(data = reef_map_df, aes(x = lon, y = lat, color = role), size = 3, alpha = 0.9) +
  scale_color_manual(values = colors) +
  coord_sf(xlim = c(149, 168), ylim = c(-32, -17), expand = FALSE) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Reef Connectivity Roles - All the reefs of interest",
    subtitle = "Sources, Sinks, Connectivity Hubs and other reefs",
    x = "Longitude", y = "Latitude", color = "Reef Role"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  )



####classify by reef class
write_xlsx(reef_map_df, "reef_map.xlsx")
#add reef class to the dataframe
reef_map_df <- read_excel("C:/Users/33778/Desktop/StageM2/reef_map.xlsx")
#calculate %
reef_percentages <- reef_map_df %>%
  group_by(class, role) %>%
  summarise(count = n()) %>%
  group_by(class) %>%
  mutate(percentage = count / sum(count) * 100)

#create plot
ggplot(reef_percentages, aes(x = class, y = percentage, fill = role)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_stack(vjust = 0.5), 
            size = 3) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = colors) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Distribution of connectivity roles by reef class",
    x = "Reef class",
    y = "Percentage",
    fill = "Reef role"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



##add temperatures to the model
#temperature data for the AOI
temp_df <- read.csv("C:/Users/33778/Desktop/StageM2/temperatures_AOI.csv", sep=";")

#mean temperature metric per location
temp_mean <- temp_data %>%
  group_by(LON, LAT) %>%
  summarise(mean_SST = mean(SST, na.rm = TRUE),
            mean_DHW = mean(DHW, na.rm = TRUE))

#search for the nearest temperature point for each reef
nn <- get.knnx(temp_mean[, c("LON", "LAT")], reef_sites[, c("lon", "lat")], k = 1)

#associate average SST/DHW values with the reef
reef_sites$mean_SST <- temp_mean$mean_SST[nn$nn.index]
reef_sites$mean_DHW <- temp_mean$mean_DHW[nn$nn.index]

#thermal distance function
dist_matrix <- as.matrix(dist(reef_sites[, c("mean_SST", "mean_DHW")], method = "euclidean"))

#convert to similarity (close values = more similar)
sigma <- mean(dist_matrix)  # ou autre valeur selon la sensibilité voulue
thermal_similarity <- exp(- (dist_matrix^2) / (2 * sigma^2))
#weighted connectivity matrix
connectivity_weighted <- connectivity_matrix_bio * thermal_similarity
#new quadrant plot
reef_df_temp <- data.frame(
  reef = colnames(connectivity_weighted),
  export = rowSums(connectivity_weighted),
  import = colSums(connectivity_weighted)
)

#metrics calcul
reef_df_temp$net_flow <- reef_df_temp$export - reef_df_temp$import
reef_df_temp$inward_degree <- reef_df_temp$import

#switching to percentile
reef_df_temp$net_flow_percentile <- ecdf(reef_df_temp$net_flow)(reef_df_temp$net_flow) * 100
reef_df_temp$inward_degree_percentile <- ecdf(reef_df_temp$inward_degree)(reef_df_temp$inward_degree) * 100

#ecological classification
reef_df_temp$role <- "Other"
reef_df_temp$role[reef_df_temp$net_flow_percentile >= 90] <- "Source"
reef_df_temp$role[reef_df_temp$net_flow_percentile <= 10] <- "Sink"
reef_df_temp$role[reef_df_temp$inward_degree_percentile >= 90] <- "Connectivity Hub"

#colors
colors <- c("Source" = "#003399", "Sink" = "#FF3399", "Connectivity Hub" = "#66CC99", "Other" = "grey80")


plot_temp <- ggplot(reef_df_temp, aes(x = inward_degree_percentile, y = net_flow_percentile)) +
  geom_point(aes(color = role), shape = 21, size = 2.2, stroke = 0.7, fill = "white") +
  scale_color_manual(values = colors) +
  geom_vline(xintercept = 90, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 90, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "black") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Quadrant Plot with Thermal Similarity",
    subtitle = "Reef roles weighted by SST/DHW similarity",
    x = "Inward Connectivity (percentile)",
    y = "Net Larval Flow (percentile)",
    color = "Reef Role"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  )

p1 + plot_temp + plot_layout(ncol = 2) +
  plot_annotation(
    title = "Comparison of Reef Roles: Without vs. With Thermal Weighting",
    theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
  )
#merge the 2 dataframes to compare
role_comparison <- data.frame(
  reef = reef_df$reef,
  role_initial = reef_df$role,
  role_thermal = reef_df_temp$role
)

#adding a column to check witch reef changed role
role_comparison$change <- ifelse(role_comparison$role_initial == role_comparison$role_thermal,
                                 "No change", "Changed")

#alluvial plot
ggplot(role_comparison,
       aes(axis1 = role_initial, axis2 = role_thermal)) +
  geom_alluvium(aes(fill = role_initial), width = 1/12, alpha=0.8) +
  geom_stratum(width = 1/12, fill = "grey30", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4, color = "white") +
  scale_x_discrete(limits = c("Initial Role", "Thermal Role"), expand = c(0.1, 0.1)) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(title = "Transitions of Reef Roles with Thermal Weighting",
       y = "Number of Reefs",
       fill = "Initial Role")

table(role_comparison$change)

#barplot
library(ggplot2)
ggplot(role_comparison, aes(x = change, fill = change)) +
  geom_bar(width = 0.6) +
  scale_fill_manual(values = c("Changed" = "#FF6666", "No change" = "#66CC99")) +
  theme_minimal() +
  labs(title = "Number of Reefs That Changed Role After Thermal Weighting",
       x = "Change of Role",
       y = "Number of Reefs")


###un raster c'est sûrement mieux pour la suite