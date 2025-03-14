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
#open file
nc_file <- "C:/Users/33778/Desktop/StageM2/DirectionVelocity_data.nc"
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
time_dates <- time_dates[format(time_dates, "%m") %in% c("09", "10", "11", "12", "01")]
unique(format(time_dates, "%Y-%m"))
table(format(time_dates, "%m"))#sept to january OK


#print first and last dates
print(range(time_dates)) ##09-2011 to 01 2021 OK
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
for (month in unique(time_months[grep("-(09|10|11|12|01)$", time_months)])) {
  print(paste("Processing month:", month))  
  
#data indices for the month
  indices <- which(time_months == month)
  
#extraction of monthly velocity data (U and V) 
  uo_month <- uo[,,indices]  
  vo_month <- vo[,,indices]
  
  
#convert to rasters
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
AOI <- ext(150, 160, -32, -20)
world <- ne_countries(scale = "medium", returnclass = "sf")
land <- st_crop(world, c(xmin = 150, xmax = 160, ymin = -32, ymax = -20))

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
  coord_sf(xlim = c(150, 160), ylim = c(-32, -20))

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


current_df <- current_df[current_df$month %in% c(9, 10, 11, 12, 1), ]
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
  150, 160,  # Longitude
  -32, -20   # Latitude
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
  geom_point(data = lord_howe, aes(x = lon, y = lat), color = "red", size = 3) +
  geom_text(data = lord_howe, aes(x = lon, y = lat, label = "Lord Howe Island"), 
            hjust = 0, vjust = -1, color = "black", size = 4) +
  scale_color_gradient(low = "lightblue", high = "purple", name = "Velocity") +
  coord_sf(xlim = c(153, 160), ylim = c(-36, -24), expand = FALSE) +  
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
reef_sf <- st_as_sf(reef_sites, coords = c("LON", "LAT"), crs = 4326)

#world
world <- ne_countries(scale = "medium", returnclass = "sf")

#map with reefs
ggplot() +
  geom_sf(data = world, fill = "antiquewhite") +
  geom_sf(data = reef_sf, color = "red", size = 2) +
  coord_sf(xlim = c(150, 160), ylim = c(-32, -20), expand = FALSE) +  
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", which_north = "true") +
  theme_minimal()

#calculate the geodesic distances between all the points
#extract data
coords <- reef_sites[, c("LON", "LAT")]

#calculate connectivity matrix (meters)
distance_matrix <- distm(coords, fun = distHaversine)

#convert to kilometers
distance_matrix_km <- distance_matrix / 1000

#add reef names
rownames(distance_matrix_km) <- reef_sites$Reef_Name
colnames(distance_matrix_km) <- reef_sites$Reef_Name

#see
print(round(connectivity_matrix, 2))
#see and save matrix
round(distance_matrix_km, 2)
write.csv(distance_matrix_km, "C:/Users/33778/Desktop/Matrice_Connectivite_Recifs.csv")

#####realistic connectivity matrix based on current speed (each cell represents the cost of travelling between two reefs via currents) : 

#defining the velocity file
file_velocity <- "C:/Users/33778/Desktop/StageM2/velocity/"

#listing velocity .tif files
file_list <- list.files(path = file_velocity, pattern = "\\.tif$", full.names = TRUE)

#load the rasters in a stack
stack_vitesse <- rast(file_list)
#calculate mean velocity pixel by pixel
mean_velocity <- mean(stack_vitesse, na.rm = TRUE)

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
writeRaster(raster_friction_directional, "C:/Users/33778/Desktop/StageM2/friction_directional.tif", overwrite = TRUE)

###circuitscape raster pairwise

resistance_matrix <- read.table("C:/Users/33778/Desktop/StageM2/connectivity_matrix_csc_resistances", quote="\"", comment.char="")
reef_sites <- read.csv("C:/Users/33778/Desktop/StageM2/reefs_format.csv", sep = ";")
reef_names <- reef_sites$id
resistance_matrix_clean <- resistance_matrix[-1, -1]
rownames(resistance_matrix_clean) <- reef_names
colnames(resistance_matrix_clean) <- reef_names
connectivity_matrix <- exp(-resistance_matrix_clean)
write.csv(connectivity_matrix, "C:/Users/33778/Desktop/StageM2/Connectivity_matrix.csv")
library(pheatmap)

pheatmap(connectivity_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "blue", "darkblue"))(100),
         border_color = "grey60")
