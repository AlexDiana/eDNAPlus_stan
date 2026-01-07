# Collecting climate variables and plotting field map

library(sf)
library(leaflet)
library(terra)
library(dplyr)
library(tidyr)
library(reticulate)
library(rgee)
library(geojsonio)
library(ggplot2)
library(patchwork)
library(maptiles)
library(ggspatial)

options(scipen = 999)  # turn off scientific notation

# for reticulate to work, I have had to manually set my python location in the Renviron file.
# i did this by usethis::edit_r_environ() and then adding RETICULATE_PYTHON=/Users/peggybevan/miniconda3/bin/python to the file, then restarting rstudio.

# Check the configuration (should point to /Users/peggybevan/miniconda3/bin/python)
# py_config()

# ee_install_set_pyenv(py_path = "/Users/peggybevan/miniconda3/bin/python",
# py_env = "base")

site_info = read.csv("Data/NatureAir_SpatialManifest.csv", fileEncoding = "latin1")
# just get long and lat for each site

# correct typos
site_info$Point_ID[site_info$Point_ID == "A" & site_info$Notes == 'North Facing'] = "AN"
site_info$Point_ID[site_info$Point_ID == "A" & site_info$Notes == 'South Facing'] = "AS"
site_info$Point_ID[site_info$Point_ID == "A" & site_info$Notes == 'Other side of roost building'] = "AO"
site_info$Point_ID[site_info$Point_ID == "A"] = "AN"

site_info$Point_ID[site_info$Point_ID == "AS1"] = "AS"

# only keep active sampler results
site_info = site_info %>%
  filter(SampleType=='Active')

# correct typo
site_info$Longitude[site_info$Point_ID=='AO' & site_info$Location=="Richmond Park"] = "-0.274390"
site_info$Longitude= as.numeric(site_info$Longitude)

site_locs = site_info %>%
  group_by(Location, Point_ID) %>%
  summarise(Latitude = mean(Latitude),
            Longitude = mean(Longitude))

cfpoints = st_as_sf(site_locs, coords = c("Longitude", "Latitude"), crs = 4326)


bat_roost = read_sf('Data/prelim_data/NatureAir.kml')
bat_roost = st_zm(bat_roost, drop = TRUE)
bat_roost <- bat_roost[bat_roost$Name == "Bat roost", ]
names(bat_roost) = c("Location", "Point_ID", "geometry")
bat_roost$Location = "Canada Farm"
bat_roost$Point_ID = "Roost"

cfpoints = rbind(cfpoints, bat_roost)

# richmond park
roostRP = read_sf('Data/prelim_data/NatureAir_RichmondPark.kml')
roostRP = st_zm(roostRP, drop = TRUE)
roostRP <- roostRP[roostRP$Name == "Bat roost", ]
names(roostRP) = c("Location", "Point_ID", "geometry")
roostRP$Location = "Richmond Park"
roostRP$Point_ID = "Roost"

cfpoints = rbind(cfpoints, roostRP)

# transform to UK national grid
cfpoints_proj <- st_transform(cfpoints, crs = 27700)

# Extracting remote sensed variables --------------------------------------
# UK CEH Land Cover Map 2024 for canada farm only
LCM = rast('Data/CEH_LCM_2024/data/LCM.tif')
plot(LCM)
unique_LCs = unique(LCM$LCM_1)

# using the CEH land cover data, converting each number into a real land-use label
lut = lut <- data.frame(
  code = unique_LCs,
  label = c(NA, "Broadleaf woodland", "Coniferous woodland", "Arable", "Improve grassland",
            "Neutral grassland", "Acid grassland", "Fen", "Heather", "Inland rock",
            "Freshwater", "Urban", "Suburban")
)

levels(LCM) = lut
plot(LCM$label, col = terrain.colors(length(lut$label)))
plot(cfpoints_proj, add = TRUE)

# crop LCM
ext = ext(cfpoints_proj[cfpoints_proj$Location=='Canada Farm',]) + 10
LCM_cf = crop(LCM, ext)
plot(LCM_cf$label, col = terrain.colors(length(lut$label)))
plot(cfpoints_proj["Point_ID"], add = TRUE, col = 'blue', pch = 20)


# land-use summary in area.
# %age of woodloand in the area:
lcm_values = values(LCM_cf$label)
lcm_values = lcm_values[!is.na(lcm_values)]
lcm_summary = as.data.frame(table(lcm_values))
lcm_summary$Percentage = (lcm_summary$Freq / sum(lcm_summary$Freq)) * 100
# 48% Fen
# 30% woodland
# 15$ Improved Grassland
# 5% Neutral Grassland
# 0.5% Urban
# 0.1% Suburban


# Extract dominant land-cover within 10m and 20m of each point ----------------------
# cfpoints_proj <- st_transform(cfpoints, crs(LCM))
# add 10m buffer
cfpoints_buffer <- st_buffer(cfpoints_proj[cfpoints_proj$Location=="Canada Farm",], dist = 10)
# add 20m buffer
cfpoints_buffer20 <- st_buffer(cfpoints_proj[cfpoints_proj$Location=="Canada Farm",] , dist = 20)

# Convert to terra vector
cfpoints_vect10 <- vect(cfpoints_buffer)
cfpoints_vect20 <- vect(cfpoints_buffer20)

# Extract all LCM values within each buffer (each point = a polygon)
extracted_vals10 <- terra::extract(LCM$label, cfpoints_vect10)
extracted_vals20 <- terra::extract(LCM$label, cfpoints_vect20)

dominant_cover10 <- extracted_vals10 %>%
  group_by(ID) %>%
  count(label, sort = TRUE) %>%
  slice(1) %>%
  rename(DominantLC_10m = label, Count = n)

dominant_cover20 <- extracted_vals20 %>%
  group_by(ID) %>%
  count(label, sort = TRUE) %>%
  slice(1) %>%
  rename(DominantLC_20m = label, Count = n)

cfpoints_final <- cfpoints_proj %>%
  filter(Location=="Canada Farm") %>%
  mutate(ID = row_number()) %>%
  left_join(dominant_cover10, by = "ID") %>%
  left_join(dominant_cover20, by = "ID") %>%
  select(-c("ID", "Count.x", "Count.y"))

bat_roost = cfpoints_final[cfpoints_final$Point_ID=="Roost",]

cfpoints_final$dist_m <- as.numeric(st_distance(cfpoints_final, bat_roost)[,1])


# Richmond park landcover
# UK CEH Land Cover Map 2024 for canada farm only
LCMRP = rast('Data/CEH_LCM_2024_RP/data/LCM.tif')
plot(LCMRP)
unique_LCs = unique(LCMRP$LCM_1)

levels(LCMRP) = lut
plot(LCMRP$label, col = terrain.colors(length(lut$label)))

# crop LCM
ext = ext(cfpoints_proj[cfpoints_proj$Location=='Richmond Park',]) + 10
LCM_rp = crop(LCMRP, ext)
plot(LCM_rp$label, col = terrain.colors(length(lut$label)))
plot(cfpoints_proj["Point_ID"], add = TRUE, col = 'blue', pch = 20)


# land-use summary in area.
# %age of woodloand in the area:
lcm_values = values(LCM_rp$label)
lcm_values = lcm_values[!is.na(lcm_values)]
lcm_summary = as.data.frame(table(lcm_values))
lcm_summary$Percentage = (lcm_summary$Freq / sum(lcm_summary$Freq)) * 100
# 71% woodland
# 0.3% coniferous woodland
# 24% Improved grassland
# 1.3% Neutral Grassland
# 1.1% Acid grassland
#  1.3% Suburban

# Extract dominant land-cover within 10m and 20m of each point ----------------------
# cfpoints_proj <- st_transform(cfpoints, crs(LCM))
# add 10m buffer
rppoints_buffer <- st_buffer(cfpoints_proj[cfpoints_proj$Location=="Richmond Park",], dist = 10)
# add 20m buffer
rppoints_buffer20 <- st_buffer(cfpoints_proj[cfpoints_proj$Location=="Richmond Park",] , dist = 20)

# Convert to terra vector
rppoints_vect10 <- vect(rppoints_buffer)
rppoints_vect20 <- vect(rppoints_buffer20)

# Extract all LCM values within each buffer (each point = a polygon)
extracted_vals10 <- terra::extract(LCMRP$label, rppoints_vect10)
extracted_vals20 <- terra::extract(LCMRP$label, rppoints_vect20)

dominant_cover10 <- extracted_vals10 %>%
  group_by(ID) %>%
  count(label, sort = TRUE) %>%
  slice(1) %>%
  rename(DominantLC_10m = label, Count = n)

dominant_cover20 <- extracted_vals20 %>%
  group_by(ID) %>%
  count(label, sort = TRUE) %>%
  slice(1) %>%
  rename(DominantLC_20m = label, Count = n)

rppoints_final <- cfpoints_proj %>%
  filter(Location=="Richmond Park") %>%
  mutate(ID = row_number()) %>%
  left_join(dominant_cover10, by = "ID") %>%
  left_join(dominant_cover20, by = "ID") %>%
  select(-c("ID", "Count.x", "Count.y"))

bat_roost = rppoints_final[rppoints_final$Point_ID=="Roost",]

rppoints_final$dist_m <- as.numeric(st_distance(rppoints_final, bat_roost)[,1])

points_final = rbind(cfpoints_final, rppoints_final)

# Using google earth engine to download stuff
np <- reticulate::import("numpy", convert = FALSE)
pd = reticulate::import('pandas')
ee = reticulate::import('ee')
# getting this working took a lot of time - i had to do a lot of authentication and setting correct paths etc.
ee$Initialize()

# To get seasonal changes in NDVI, we need a dates of each field visit

# March, May, June, July, October
cfdates = data.frame(Location = "Canada Farm", Visit = 1:5, Date = c("04/03/2025", "06/05/2025", "10/06/2025", "29/07/2025", "08/10/2025"))
rpdates = data.frame(Location = "Richmond Park", Visit = 1:4, Date = c("01/07/2025", "06/08/2025", "03/09/2025", "01/10/2025"))
dates = rbind(cfdates, rpdates)
dates$Date = as.Date(dates$Date, format = "%d/%m/%Y")
dates

# this selects images that have less than 40% cloud cover - this is quite high, but neccessary for Uk winters!
s2 = ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED")$filter(ee$Filter$lt("CLOUDY_PIXEL_PERCENTAGE", 40))$
  map(function(image) {
    ndvi = image$normalizedDifference(c("B8", "B4"))$rename("NDVI")
    image$addBands(ndvi)
  })

get_ndvi_date = function(point, date, buffer_m = 20) {
  # revist time for sentinal is every 5 days, but the data available is much less than this
  # need to add a window of at least 14 days to make sure all visits are included
  # Sampling always began on tuesday
  # I've set it so that it will look for data within three weeks of the start date.
  start = ee$Date(date)
  begin = start$advance(-14, "day")
  end = start$advance(21, "day")

  # Filter collection to date range and point location
  img <- s2$
    filterDate(begin, end)$
    filterBounds(point$geometry())$
    sort("system:time_start")$
    first()

  # getInfo() converts EE object to R list or NULL
  img_info <- tryCatch(img$getInfo(), error = function(e) NULL)
  if (is.null(img_info)) {
    message(sprintf("No image found for point %s on date %s",
                    point$get("Name")$getInfo(), as.character(date)))
    return(NA_real_)
  }
  buffer_geom <- point$geometry()$buffer(buffer_m)

  # Extract mean NDVI value within buffer around point
  ndvi_value <- img$select("NDVI")$reduceRegion(
    reducer = ee$Reducer$mean(),
    geometry = buffer_geom,
    scale = 10
  )$get("NDVI")

  # Return numeric value or NA if no image found
  ee$Number(ndvi_value)$getInfo()
}

cfpoints_d = crossing(cfpoints[cfpoints$Location=="Canada Farm",-1], dates[dates$Location=="Canada Farm",-1])
rppoints_d = crossing(cfpoints[cfpoints$Location=="Richmond Park",-1], dates[dates$Location=="Richmond Park",-1])
cfpoints_d$Location = "Canada Farm"
rppoints_d$Location = "Richmond Park"
points_d = rbind(cfpoints_d, rppoints_d)

points_d$Date <- as.character(points_d$Date)
points_d <- st_as_sf(points_d, crs = 4326)
points_ee <- sf_as_ee(points_d)

# Loop over each point and date combo
results <- lapply(1:nrow(points_d), function(i) {
  print(paste0("i = ", i))
  point <- points_ee$filter(ee$Filter$eq("Point_ID", points_d$Point_ID[i]))$first() #get point
  ndvi <- get_ndvi_date(point, as.character(points_d$Date[i])) #extract ndvi
  data.frame(Location = points_d$Location[i], site = points_d$Point_ID[i], date = points_d$Date[i], ndvi = ndvi) #return data frame
})

ndvi_df <- do.call(rbind, results)
head(ndvi_df)
# add visit
ndvi_df$Visit = c(rep(1:5, 20), rep(1:4, 18))

write.csv(ndvi_df, "output/ndvi20m.csv", row.names = FALSE)


# Climate variables from ERA5 map -----------------------------------------

# Climate variables from the ERA5 map
# The ERA5 map has 30km resolution, so i will just get one value for the whole site.
# because our samples were only taken at night (6pm-9am), I am going to take hourly metrics
# and summarise over the time frame.

# Example: ERA5-Land hourly
era5l <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY")

# Helper to add local hour (Europe/London) as a numeric property
# The data comes in UTC time zone so this is adjusting it
add_local_hour <- function(img) {
  hr <- ee$Number$parse(ee$Date(img$get("system:time_start"))
                        $format("H", "Europe/London"))
  img$set("hour_local", hr)
}

# we want Temp, windspeed, precipitation, relative humidity
with_bands <- era5l$
  select(c("temperature_2m","dewpoint_temperature_2m","u_component_of_wind_10m",
           "v_component_of_wind_10m","total_precipitation"))$
  map(function(img){
    t  <- img$select("temperature_2m") # Kelvin
    td <- img$select("dewpoint_temperature_2m")# Kelvin
    # RH is calculated from the Magnus formula (OK for reanalysis scale)
    rh <- td$subtract(273.15)$rename("td_c") #convert to celcius
    tC <- t$subtract(273.15)$rename("t2m_c") #convert to celcius
    es <- tC$multiply(17.62)$divide(tC$add(243.12))$exp()$multiply(6.112) # hPa
    e  <- rh$multiply(17.62)$divide(rh$add(243.12))$exp()$multiply(6.112)
    RH <- e$divide(es)$multiply(100)$rename("rh_pct")
    # windspeed is calculatted by adding these two wind components
    wspd <- img$expression("sqrt(u*u + v*v)", list(u = img$select("u_component_of_wind_10m"),
                                                   v = img$select("v_component_of_wind_10m")))$
      rename("wind10m_ms")

    # Add the new bands
    img_with_bands <- img$addBands(list(tC, RH, wspd))
    add_local_hour(img_with_bands)
  })


# Function to summarise one night (start date inclusive)
summarise_night <- function(geom, date_str, start_hour = 18, end_hour = 9) {
  start <- ee$Date(date_str)
  end   <- start$advance(1, "day")  # next day boundary
  # Keep hours 18–23 and 0–9 local
  night_coll <- with_bands$
    filterDate(start, end)$
    filterBounds(geom)$
    filter(ee$Filter$Or(
      ee$Filter$rangeContains("hour_local", start_hour, 23),
      ee$Filter$rangeContains("hour_local", 0, end_hour)
    ))
  img_mean <- night_coll$mean()
  img_sum  <- night_coll$sum()
  # Reduce to your buffer/point
  reducers <- ee$Reducer$mean()$combine(reducer2 = ee$Reducer$min(), sharedInputs = TRUE)$
    combine(reducer2 = ee$Reducer$max(), sharedInputs = TRUE)
  stats <- img_mean$select(c("t2m_c","rh_pct","wind10m_ms"))$
    addBands(img_sum$select("total_precipitation"))$
    reduceRegion(
      reducer  = reducers,
      geometry = geom,
      scale    = 9000,   # ~9 km native; you can set 10000 for speed
      bestEffort = TRUE
    )
  stats
}

# get everysingle date for each site

dates_all = site_info[,c('Location', 'DateCollected')] %>%
  distinct()

dates_all$Date = as.Date(dates_all$DateCollected, format = "%d/%m/%Y")

# Because the map is 30km resolution, we only need one geolocation for each site.
# I'll use the bat roost lat/long as this is most central
points_d <- st_as_sf(points_d, crs = 4326)
cfbatroost = points_d[points_d$Point_ID=='Roost',][1,]
rpbatroost = points_d[points_d$Point_ID=='Roost',][6,]
cfbatroost_ee <- sf_as_ee(cfbatroost)
rpbatroost_ee <- sf_as_ee(rpbatroost)


# Get first feature from the FeatureCollection
feat1cf <- cfbatroost_ee$first()
# Extract the geometry from that feature
geom1 <- feat1$geometry()

# Get first feature from the FeatureCollection
feat1rp <- rpbatroost_ee$first()
# Extract the geometry from that feature
geom2 <- feat1rp$geometry()

era_results_list <- lapply(1:nrow(dates_all), function(i) {
  if (dates_all$Location[i]=="Canada Farm") {
    geom = geom1
  } else {
    geom = geom2
  }
  era_list = summarise_night(
    geom = geom,
    date_str = as.character(dates_all$Date[i])
  )$getInfo()
  era_df = as.data.frame(era_list)
  era_df$Date = dates_all$Date[i]
  era_df$Location = dates_all$Location[i]
  era_df
})

era_all <- do.call(rbind, era_results_list)
era_all

# remove extra precipitation columns
era_all1 = era_all %>%
  mutate(total_precipitation_mm = total_precipitation_mean*1000) %>%
  select(Location, Date, rh_pct_min, t2m_c_mean, total_precipitation_mm, wind10m_ms_mean)

# Combine all metadata into one dataset -----------------------------------

# Now combine the whole dataset
# We need one row for each site, visit, date combination
# site, distance, landcover, climate vars

# site_info currently has one row for every sample.

# add distance and landcover from cfpoints_final
# remove geometry
points_ngeom <- st_drop_geometry(points_final)
site_info_vars1 = left_join(site_info, points_ngeom, by = c("Location", "Point_ID"))

# add climate vars
site_info_vars1$Date = as.Date(site_info_vars1$DateCollected, format = "%d/%m/%Y")
site_info_vars2 = left_join(site_info_vars1, era_all1, by = c('Location', "Date"))

# add ndvi
names(ndvi_df)[2] = 'Point_ID'
site_info_vars3 = left_join(site_info_vars2, ndvi_df[,-3], by = c("Location", "Point_ID", "Visit"))

# save!
write.csv(site_info_vars3, 'output/allsitevars.csv', row.names = FALSE)
