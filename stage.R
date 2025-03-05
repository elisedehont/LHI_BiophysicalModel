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
library(dplyr)
library(ggplot2)
library(gganimate)
library(transformr)
library(gifski)
library(sf)
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
library(rnaturalearth)
library(rnaturalearthdata)
library(geosphere)
library(ggforce)
#open file
nc_file <- "C:/Users/33778/Desktop/StageM2/data/DirectionVelocity_data.nc"
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

#print first and last dates
print(range(time_dates)) ##06-2011 to 05 2021 OK
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
  uo_month_raster <- raster(t(apply(uo_month, c(1,2), mean, na.rm = TRUE)))
  vo_month_raster <- raster(t(apply(vo_month, c(1,2), mean, na.rm = TRUE)))
  
  #Defining the raster extent
  extent(uo_month_raster) <- extent(min(longitude), max(longitude), min(latitude), max(latitude))
  extent(vo_month_raster) <- extent(min(longitude), max(longitude), min(latitude), max(latitude))
  #compute velocity and direction
  velocity_month <- sqrt(uo_month_raster^2 + vo_month_raster^2)
  direction_month <- atan2(vo_month_raster, uo_month_raster) * (180 / pi)
  
  #results storage
  monthly_avg_list[[month]] <- list(velocity = velocity_month, direction = direction_month)
  
  # Saving rasters for each month
  writeRaster(velocity_month, paste0("velocity_", month, ".tif"), format = "GTiff", overwrite = TRUE)
  writeRaster(direction_month, paste0("direction_", month, ".tif"), format = "GTiff", overwrite = TRUE)
}

#inspection OK
nc_close(nc)
###rasters OK to work with


#load currents direction and velocity rasters for a month
#define the month and year
month_year <- "2013-05"
#rasters download
SCD <- rast("C:/Users/33778/Desktop/StageM2/velocity_2013-05.tif")
SCV <- rast("C:/Users/33778/Desktop/StageM2/direction_2013-05.tif")

#Define the area of interest (AOI)
AOI <- ext(153.9583, 162.9584, -35.04167, -25.95833)
world <- ne_countries(scale = "medium", returnclass = "sf")
land <- st_crop(world, c(xmin = 153.9583, xmax = 162.9584, ymin = -35.04167, ymax = -25.95833))

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
  coord_sf(xlim = c(153.9583, 162.9584), ylim = c(-35.04167, -25.95833))

######currents characterisation####


#Lord Howe Island
lord_howe <- data.frame(
  lon = 159.0833,
  lat = -31.5500
)

#buffer around LHI (~1°)
buffer_radius <- 1  #degrees
buffer_extent <- ext(
  lord_howe$lon - buffer_radius, lord_howe$lon + buffer_radius,
  lord_howe$lat - buffer_radius, lord_howe$lat + buffer_radius
)

#Data extraction in the area of interest (buffer)
SCD_buffer <- crop(SCD_crop, buffer_extent)  # Vitesse
SCV_buffer <- crop(SCV_crop, buffer_extent)  # Direction


current_df <- current_df[current_df$month %in% c(9, 10, 11, 12, 1), ]
current_df$month <- sprintf("%02d", as.numeric(current_df$month))
time_values <- as.Date(paste0(current_df$year, "-", current_df$month, "-01"), format = "%Y-%m-%d")


# Dataframe
vel_df <- as.data.frame(SCD_buffer, xy = TRUE)
dir_df <- as.data.frame(SCV_buffer, xy = TRUE)
colnames(vel_df) <- c("x", "y", "velocity")
colnames(dir_df) <- c("x", "y", "direction")

# Merge data and delete NAs
current_df <- na.omit(merge(vel_df, dir_df, by = c("x", "y")))

# Calculating the direction of points towards Lord Howe
bearing_to_lord_howe <- function(lon, lat) {
  bearing(c(lon, lat), c(lord_howe$lon, lord_howe$lat))  
}

current_df$bearing_to_lh <- mapply(bearing_to_lord_howe, current_df$x, current_df$y)

# Determining whether a current is flowing into or out of Lord Howe
angle_diff <- abs(current_df$direction - current_df$bearing_to_lh)  # Difference between current direction and direction towards LHI
current_df$type <- ifelse(angle_diff < 90, "Arriving", "Leaving")  # <90° = incoming current, >90° = outgoing

# Statistics
current_stats <- current_df %>%
  group_by(type) %>%
  summarise(
    mean_velocity = mean(velocity, na.rm = TRUE),
    sd_velocity = sd(velocity, na.rm = TRUE),
    mean_direction = mean(direction, na.rm = TRUE)
  )

# Results
print(current_stats)
#1 Arriving         mean velo :0.191       sd dir : 0.132         mean dir: -53.7
#2 Leaving          mean velo: 0.361       sd vel : 0.128         mean dir : -58.5

#current vectors map around LHI
ggplot() +
  geom_segment(data = current_df_filtered, 
               aes(x = x, y = y,
                   xend = x + sin(deg2rad(direction)) * velocity * 0.1,
                   yend = y + cos(deg2rad(direction)) * velocity * 0.1,
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
  geom_point(data = current_df_filtered, aes(x = x, y = y, color = type), alpha = 0.7) +
  geom_point(data = lord_howe, aes(x = lon, y = lat), color = "red", size = 3) +
  scale_color_manual(values = c("Arriving" = "blue", "Leaving" = "orange")) +
  labs(title = "Arriving vs Leaving Currents (Sep-Jan)", x = "Longitude", y = "Latitude", color = "Current Type") +
  theme_minimal()

#histogram of currents speed
ggplot(current_df_filtered, aes(x = velocity, fill = type)) +
  geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
  scale_fill_manual(values = c("Arriving" = "blue", "Leaving" = "orange")) +
  labs(title = "Velocity Distribution (Sep-Jan)", 
       x = "Velocity (m/s)", y = "Count", fill = "Current Type") +
  theme_minimal()































###temperature###
#open file
nc_file2 <- "C:/Users/33778/Desktop/StageM2/data/seawaterTemperatures.nc"
nc2 <- nc_open(nc_file2)
# Print variable names and dimensions to check
print(nc2$var)
print(nc2$dim)
#extract dim and variables values
time_values2 <- ncvar_get(nc2, "time")  # Check if 'time' is the correct name
time_units2 <- ncatt_get(nc2, "time", "units")$value  # Get unit metadata
print(time_units2)
depth_vals <- ncvar_get(nc2, "depth")
temp <- ncvar_get(nc2, "thetao")  # Température
#convert time values to dates
time_origin2 <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
time_dates2 <- time_origin2 + time_values2
#print first and last dates
print(range(time_dates2)) ##06-2011 to 05 2021 OK
diff_time <- diff(time_dates2)  # Compute time differences
print(unique(diff_time))  # Print unique time intervals
print(head(time_dates2, 20))  #first 20 time points OK
print(tail(time_dates2, 20))  #last 20 time points OK

#Convert dates to year-month format
time_months2 <- format(time_dates2, "%Y-%m")

#Create an empty list to store monthly averages
monthly_avg_list2 <- list()
#longitude and latitude
longitude2 <- ncvar_get(nc2, "longitude")
latitude2 <- ncvar_get(nc2, "latitude")
thetao <- ncvar_get(nc2, "thetao")
#loop for each month
for (month in unique(time_months2)) {
  print(paste("Processing month:", month))
  
  #data indices for the month
  indices2 <- which(time_months2 == month)
  print(paste("Month:", month, "- Indices:", length(indices2)))
  
  if (length(indices2) > 1) {
    #Extraction of surface temperatures 
    temp_month <- apply(thetao[,,indices2, drop = FALSE], c(1,2), mean, na.rm = TRUE)
  } else {
    temp_month <- thetao[,,1,indices2]  # Si un seul point, pas besoin de moyenne
  }
  
  #convert to raster
  temp_month_raster <- raster(t(temp_month))
  extent(temp_month_raster) <- extent(min(longitude2), max(longitude2), min(latitude2), max(latitude2))
  
  #store raster
  monthly_avg_list2[[month]] <- temp_month_raster
  
  #save
  writeRaster(temp_month_raster, paste0("temperature_", month, ".tif"), format = "GTiff", overwrite = TRUE)
}

# Fermeture du fichier NetCDF
nc_close(nc2)
##temperature rasters ready to use

#load the temperature for the month
temp_raster <- rast(monthly_avg_list2[[month_year]])

#checking+harmonising resolutions
if (!inherits(monthly_avg_list2[[month_year]], "SpatRaster")) {
  temp_raster <- rast(monthly_avg_list2[[month_year]])  #convert to raster
} else {
  temp_raster <- monthly_avg_list2[[month_year]]
}

print(temp_raster)
SCD_crop <- rast("C:/Users/33778/Desktop/StageM2/velocity_2015-06.tif")
SCV_crop <- rast("C:/Users/33778/Desktop/StageM2/direction_2015-06.tif")
temp_raster <- resample(temp_raster, SCD_crop, method = "bilinear")

temp_raster <- crop(temp_raster, AOI)

#Stacking aligned rasters (Temperature, velocity, Direction)
stacked_rasters <- c(temp_raster, SCD_crop, SCV_crop)
names(stacked_rasters) <- c("Temperature", "Velocity", "Direction")
#convert to data.frame
stacked_df <- as.data.frame(stacked_rasters, xy = TRUE)

#delete NAs
stacked_df <- na.omit(stacked_df)

#show data
head(stacked_df)
library(ggplot2)

ggplot() +
  geom_raster(data = stacked_df, aes(x = x, y = y, fill = Temperature)) +
  geom_segment(data = stacked_df, 
               aes(x = x, y = y,
                   xend = x + sin(deg2rad(Direction)) * Velocity * 0.1,
                   yend = y + cos(deg2rad(Direction)) * Velocity * 0.1,
                   color = Velocity),
               arrow = arrow(length = unit(0.2, "cm"))) +
  scale_fill_gradient(low = "blue", high = "red", name = "Temperature") +
  scale_color_gradient(low = "lightblue", high = "purple", name = "Velocity") +
  geom_point(data = lord_howe, aes(x = lon, y = lat), color = "red", size = 3) +
  geom_text(data = lord_howe, aes(x = lon, y = lat, label = "Lord Howe Island"), 
            hjust = 0, vjust = -1, color = "black", size = 3) +
  labs(title = paste("Temperature & Currents -", month_year),
       x = "Longitude", y = "Latitude") +
  theme_minimal()
##DHW
#open file
nc_file3 <- "C:/Users/33778/Desktop/StageM2/DHW1.nc"
nc3 <- nc_open(nc_file3)

# Print variable names and dimensions to check
print(nc3$var)
print(nc3$dim)
# Extraction des dimensions et des valeurs
time_values3 <- ncvar_get(nc3, "time")  # Vérifier le nom correct de la variable temps
time_units3 <- ncatt_get(nc3, "time", "units")$value  
print(time_units3)

#convert time values to dates
time_origin3 <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
time_dates3 <- time_origin + time_values

#print first and last dates
print(range(time_dates3)) ##05-2011 to 05 2021 OK
diff_time3 <- diff(time_dates3)  # Compute time differences
print(unique(diff_time3))  # Print unique time intervals
print(head(time_dates3, 20))  #first 20 time points OK
print(tail(time_dates3, 20))  #last 20 time points OK
#Convert dates to year-month format
time_months3 <- format(time_dates3, "%Y-%m")

#Create an empty list to store monthly averages
dhw_rasters <- list()
#longitude and latitude
dhw <- ncvar_get(nc3, "degree_heating_week")  # Variable DHW
lon <- ncvar_get(nc3, "longitude")
lat <- ncvar_get(nc3, "latitude")
dim(dhw) #3dim ok


#loop for each month
for (month in unique(time_months3)) {
  idx <- which(time_months3 == month)  
  
  if (length(dim(dhw)) == 3) {
    dhw_month <- apply(dhw[,,idx], c(1,2), mean, na.rm = TRUE)  
  } else {
    dhw_month <- dhw[,,idx]  
  }
  
  dhw_month <- as.matrix(dhw_month)
  r3 <- raster(t(dhw_month), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
  projection(r3) <- CRS("+proj=longlat +datum=WGS84")
  
  dhw_rasters[[month]] <- r3  
}

#save
for (month in names(dhw_rasters)) {
  writeRaster(dhw_rasters[[month]], filename=paste0("DHW_", month, ".tif"), format="GTiff", overwrite=TRUE)
}
##DHW rasters ready to use
#graph
# Load necessary libraries
library(terra)
library(sp)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggplot2)

# Define the month and year
month_year <- "2013-05"

# Rasters download
SCD <- rast("C:/Users/33778/Desktop/StageM2/velocity_2013-05.tif")
SCV <- rast("C:/Users/33778/Desktop/StageM2/direction_2013-05.tif")
DHW <- rast("C:/Users/33778/Desktop/StageM2/DHW_2013-05.tif")

# Define the area of interest (AOI)
AOI <- ext(153.9583, 162.9584, -35.04167, -25.95833)

# Load world map data
world <- ne_countries(scale = "medium", returnclass = "sf")
land <- st_crop(world, c(xmin = 153.9583, xmax = 162.9584, ymin = -35.04167, ymax = -25.95833))

# Cut the rasters
SCD_crop <- crop(SCD, AOI)
SCV_crop <- crop(SCV, AOI)
DHW_crop <- crop(DHW, AOI) #cropp DHW

# Convert into data.frame
vel_df <- as.data.frame(SCD_crop, xy = TRUE)
dir_df <- as.data.frame(SCV_crop, xy = TRUE)
dhw_df <- as.data.frame(DHW_crop, xy = TRUE) #DHw data in df

# Naming columns
colnames(vel_df) <- c("x", "y", "velocity")
colnames(dir_df) <- c("x", "y", "direction")
colnames(dhw_df) <- c("x", "y", "dhw") #DHw colnames

# Merge the two data.frames and delete the NAs
current_df <- na.omit(merge(vel_df, dir_df, by = c("x", "y")))

# Merge current data with DHW data
current_df <- merge(current_df, dhw_df, by = c("x", "y"), all.x = TRUE)

# Function for calculating arrow coordinates
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

# Reducing the arrows density
step <- 5
current_df_sampled <- current_df[seq(1, nrow(current_df), by = step), ]

# Add Lord Howe island-can be done with rnaturalearthdata later
lord_howe <- data.frame(
  lon = 159.0833, # Longitude Lord Howe Island
  lat = -31.5500 # Latitude Lord Howe Island
)

# Plot
ggplot() +
  geom_raster(data = current_df, aes(x = x, y = y, fill = dhw)) + # Add DHW as raster
  scale_fill_gradientn(colors = c("blue", "yellow", "red"), name = "DHW (°C-weeks)") + # color scale for DHW
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
  labs(title = paste("Current Vectors and Degree Heating Week -", month_year),
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  coord_sf(xlim = c(153.9583, 162.9584), ylim = c(-35.04167, -25.95833))


 ####example for larvae tracking
#larvae data simulation
set.seed(123)  #reproductibility
n_larves <- 150
n_steps <- 40
larves_data <- data.frame(
  particle_id = rep(1:n_larves, each = n_steps),
  time = rep(1:n_steps, times = n_larves),
  x = numeric(n_larves * n_steps),
  y = numeric(n_larves * n_steps)
)

#initial position of the larvae
larves_data$x[larves_data$time == 1] <- rnorm(n_larves, mean = 159.0833, sd = 0.1)
larves_data$y[larves_data$time == 1] <- rnorm(n_larves, mean = -31.5500, sd = 0.1)

#movement simulation
for (i in 2:n_steps) {
  prev_positions <- larves_data[larves_data$time == (i-1), ]
  new_positions <- larves_data[larves_data$time == i, ]
  
  #random motion influenced by currents
  new_positions$x <- prev_positions$x + rnorm(n_larves, mean = 0.05, sd = 0.02)
  new_positions$y <- prev_positions$y + rnorm(n_larves, mean = 0.05, sd = 0.02)
  
  larves_data[larves_data$time == i, c("x", "y")] <- new_positions[, c("x", "y")]
}

#basic graph
p <- ggplot() +
  geom_segment(data = current_df_sampled, 
               aes(x = x, y = y, 
                   xend = x + sin(deg2rad(direction)) * velocity * 0.1,
                   yend = y + cos(deg2rad(direction)) * velocity * 0.1,
                   color = velocity),
               arrow = arrow(length = unit(0.2, "cm"))) +
  geom_point(data = lord_howe, aes(x = lon, y = lat), color = "red", size = 3) +
  geom_text(data = lord_howe, aes(x = lon, y = lat, label = "Lord Howe Island"), 
            hjust = 0, vjust = -1, color = "black", size = 3) +
  scale_color_gradient(low = "lightblue", high = "purple", name = "Velocity") +
  labs(title = paste("Current Vectors and Larval Dispersion -", month_year),
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  coord_fixed(ratio = 1)

#add larvae and animation
p_animated <- p +
  geom_point(data = larves_data, aes(x = x, y = y, group = particle_id), 
             color = "green", size = 1, alpha = 0.7) +
  transition_time(time) +
  labs(subtitle = "Time step: {frame_time}")

#generate
anim <- animate(p_animated, nframes = n_steps, fps = 2)
anim
#save
anim_save("larval_dispersion.gif", animation = anim)

