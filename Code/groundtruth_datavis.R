# raw data vis
library(ggplot2)
library(dplyr)
library(lubridate)


# data
sd = read.csv("output/spatial_data_longform_withquant_speciesID.csv")
siteinfo = read.csv("output/all_site_vars2.csv")

# total library size
total_libsize = sum(sd$read_count[sd$Location == 'Canada Farm'])
# 0.005% of total library size
thresh = total_libsize * 0.00005

# sum lib size of each sample
libsizes = sd %>%
  group_by(Sample) %>%
  summarise(total_libsize = sum(read_count))
libsizes$thresh = libsizes$total_libsize * 0.0005

mean(libsizes$thresh)

sd = sd %>%
  mutate(Month = case_when(
    Location == 'Canada Farm' & Visit == 1 ~ 'March',
    Location == 'Canada Farm' & Visit == 2 ~ 'May',
    Location == 'Canada Farm' & Visit == 3 ~ 'June',
    Location == 'Canada Farm' & Visit == 4 ~ 'July',
    Location == 'Canada Farm' & Visit == 5 ~ 'October',
    Location == 'Richmond Park' & Visit == 1 ~ 'July',
    Location == 'Richmond Park' & Visit == 2 ~ 'August',
    Location == 'Richmond Park' & Visit == 3 ~ 'September',
    Location == 'Richmond Park' & Visit == 4 ~ 'October',
  ))

sd$Month = factor(sd$Month, levels = c('March', 'May', 'June', 'July', 'August', 'September', 'October'))

bat_detections = sd %>%
  filter(Order == 'Chiroptera')
# need to add one row for march where lib size = 0 so that it appears in the plot

# Remove all bat detections that are less than 0.05% of sample library size
# Notes - this removes non-detections too/
# Applying these thresholds filters out 49
# bat detections (out of 162) and the new minimum read count for
# a bat would now be 22. The 49 lost detection are almost all
# pipistrelles, except for one Nyctalus (4 reads) and one
# Barbastelle (2 reads).
bat_detections = bat_detections %>%
  left_join(libsizes, by = c("Sample")) %>%
  filter(read_count >= thresh)

marchrow = bat_detections[1,] %>%
  mutate(Month = 'March', read_count = 0, total_libsize = 0, thresh = 0)

bat_detections = rbind(bat_detections, marchrow)

# calculate read numbers as a proportion of library size
bat_detections = bat_detections %>%
  mutate(read_prop = read_count / total_libsize)

# bin distance into groups
bat_detections = bat_detections %>%
  mutate(Distance_bin = case_when(
    dist_m < 20 ~ '0',
    dist_m >= 20 & dist_m < 60 ~ '50',
    dist_m >= 60 & dist_m < 110 ~ '100',
    dist_m >= 110  ~ '150',
  ))

bat_detections$Month = factor(bat_detections$Month, levels = c('March', 'May', 'June', 'July', 'August', 'September', 'October'))


# summary of all bat detections -------------------------------------------

# table of detections: how many samples had bats in?
nrow(bat_detections[bat_detections$read_count > 0,])

# how many detections of each species group?
# group by sample ID and species ID to find total detections

sumtab = sd %>%
  group_by(Location, Sample, SpeciesID) %>%
  summarise(total_reads = sum(read_count))

table(sumtab$SpeciesID[sumtab$total_reads > 0], sumtab$Location[sumtab$total_reads > 0])

month_cols = c(
  "March" = "#440154",
  "May" = "#31688e",
  "June" = "#21908CFF",
  "July" = "#35b779",
  "August" = "#fde725",
  "September" = "#fbb03b",
  "October" = "#ff7f0e"
)
# x = distance, y = normalised read count, faceted by species coloured by visit

ggplot(bat_detections[bat_detections$Location=="Canada Farm",], aes(x = dist_m, y = read_prop, colour = Month)) +
  geom_point() +
  facet_wrap(~ SpeciesID, scales = 'free_x')

ggplot(bat_detections[bat_detections$Location=="Richmond Park",], aes(x = dist_m, y = read_prop, colour = Month)) +
  geom_point() +
  facet_wrap(~ SpeciesID, scales = 'free_x')

# now lets calculate with point range
# averaging read proportions at each site and visit. so each data point contains 3 days each with 3 replicates (9)
bat_detections_sum = bat_detections %>%
  group_by(Location, Distance_bin, Month, SpeciesID) %>%
  summarise(avg_reads = mean(read_prop, na.rm = TRUE), total_positives = sum(read_count > 0),
            sd_reads = sd(read_prop, na.rm = TRUE)) %>%
  ungroup()
bat_detections_sum$lci = bat_detections_sum$avg_reads - bat_detections_sum$sd_reads
bat_detections_sum$uci = bat_detections_sum$avg_reads + bat_detections_sum$sd_reads
bat_detections_sum$lci[bat_detections_sum$lci < 0] = 0.01
bat_detections_sum$uci[bat_detections_sum$uci > 1] = 0.99

bat_detections_sum$Distance_bin = factor(bat_detections_sum$Distance_bin, levels = c('0', '50', '100', '150'))

ggplot(bat_detections_sum[!is.na(bat_detections_sum$lci),], aes(x = Distance_bin, y = avg_reads, ymin = lci, ymax = uci, colour = Month)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  facet_wrap(Location ~ SpeciesID)

species_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# lets try raw data, with point range on top.
ggplot(bat_detections, aes(x = Distance_bin, y = read_prop, colour = SpeciesID)) +
  geom_point(alpha = 0.3, position = position_dodge(width = 0.7)) +
  geom_pointrange(data = bat_detections_sum,
                  aes(x = Distance_bin, y = avg_reads,
                      ymin = lci, ymax = uci, colour = SpeciesID),
                  position = position_dodge(width = 0.7), size = 0.3) +
  facet_grid(Location ~ Month) +
  scale_colour_manual(values = species_colors) +
  theme_bw() +
  labs(x = 'Distance (m)', y = 'Normalised read count', colour = 'Bat taxa') +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        major.grid = element_blank(),
        minor.grid = element_line(size = 0.5),
        legend.position = "bottom")


ggplot(bat_detections_sum, aes(x = Distance_bin, y = avg_reads, colour = SpeciesID, fill = SpeciesID)) +
  geom_col(alpha = 0.7, position = position_dodge2(preserve = 'single')) +
  geom_errorbar(aes(ymin = lci, ymax = uci),
                position = position_dodge2(preserve = 'single'), size = 0.3) +
  facet_grid(Location ~ Month) +
  scale_colour_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  theme_bw() +
  labs(x = 'Distance (m)', y = 'Normalised read count', colour = 'Bat taxa', fill = 'Bat taxa') +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        major.grid = element_blank(),
        minor.grid = element_line(size = 0.5),
        legend.position = "bottom")


# final try - can i plot it in the same way i do for the acoustic data?
# using geom_smooth
ggplot(data = bat_detections, aes(x = dist_m, y = read_prop, colour = SpeciesID)) +
  geom_smooth(method = "loess", se = TRUE) +
  facet_grid(Location ~ Month) +
  scale_colour_manual(values = species_colors) +
  theme_bw() +
  labs(x = 'Distance (m)', y = 'Normalised read count', colour = 'Bat taxa') +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        major.grid = element_blank(),
        minor.grid = element_line(size = 0.5),
        legend.position = "bottom")



# find total detection number for each species
bat_dets_total = bat_detections_sum %>%
  group_by(Location, SpeciesID) %>%
  summarise(total_reads = sum(total_reads), total_positives = sum(total_positives)) %>%
  filter(total_positives == 0)

# i want to remove species that were never detected at each location
bat_detections_sum = bat_detections_sum %>%
  anti_join(bat_dets_total, by = c("Location", "SpeciesID"))

bat_detections_sum$Distance_bin = factor(bat_detections_sum$Distance_bin, levels = c('0', '50', '100', '150'))
# plot read count in each distance bin


convert to month
visit_cols <- c(
  "1" = "#440154",
  "2" = "#31688e",
  "3" = "#21908CFF",
  "4" = "#35b779",
  "5" = "#fde725"
)

month_cols = c(
  "March" = "#440154",
  "May" = "#31688e",
  "June" = "#21908CFF",
  "July" = "#35b779",
  "August" = "#fde725",
  "September" = "#fbb03b",
  "October" = "#ff7f0e"
)


bat_detections %>%
  filter(Location == 'Canada Farm' & read_count > 0) %>%
  nrow()
bat_detections %>%
  filter(Location == 'Richmond Park' & read_count > 0) %>%
  nrow()

#total species
bat_detections$SpeciesID2 = ifelse(bat_detections$Species == "", bat_detections$Genus, bat_detections$Species)

bat_detections_sum2 = bat_detections %>%
  group_by(Location, Distance_bin, Month, SpeciesID2) %>%
  summarise(total_reads = sum(read_count), total_positives = sum(read_count > 0)) %>%
  ungroup()

# frequency table of samples per species per visit per location
bat_detections_sum2 %>%
  filter(Location == 'Canada Farm') %>%
  group_by(SpeciesID2) %>%
  summarise(total_positives = sum(total_positives)) %>%
  arrange(desc(total_positives))

bat_detections_sum2 %>%
  filter(Location == 'Richmond Park') %>%
  group_by(SpeciesID2) %>%
  summarise(total_positives = sum(total_positives)) %>%
  arrange(desc(total_positives))

bat_detections_sum2 %>%
  group_by(SpeciesID2) %>%
  summarise(total_positives = sum(total_positives)) %>%
  arrange(desc(total_positives))


march_insert = bat_detections_sum[1,]
march_insert$Month = 'March'
bat_detections_sum = rbind(bat_detections_sum, march_insert)

ggplot(bat_detections_sum[bat_detections_sum$Location=='Canada Farm',],
       aes(x = Distance_bin, y = log(total_reads), fill = Month)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = 'single')) +
  # geom_line(aes(group = interaction(SpeciesID, Visit)), alpha = 0.3) +
  facet_wrap(~SpeciesID, scales = 'free') +
  scale_y_continuous(limits = c(0,12), expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  labs(title = "Canada Farm", x = "Distance (m)", y = "Log(Total Reads)", fill = "Visit") +
  scale_fill_manual(values = month_cols, breaks = levels(bat_detections_sum$Month)) +
  # only include labels where total_positives > 0
  theme(text = element_text(size = 18,),
        strip.background = element_blank(),
        strip.text = element_text(size = 18))


# plot read count in each distance bin
ggplot(bat_detections_sum[bat_detections_sum$Location=='Richmond Park',],
       aes(x = Distance_bin, y = log(total_reads), fill = Month)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = 'single')) +
  # geom_line(aes(group = interaction(SpeciesID, Visit)), alpha = 0.3) +
  facet_wrap(~SpeciesID, scales = 'free') +
  scale_y_continuous(limits = c(0,12), expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  labs(title = "Richmond Park", x = "Distance (m)", y = "Log(Total Reads)", fill = "Visit") +
  scale_fill_manual(values = month_cols, breaks = levels(bat_detections_sum$Month)) +
  # only include labels where total_positives > 0
  theme(text = element_text(size = 18,),
        strip.background = element_blank(),
        strip.text = element_text(size = 18))

# plot bats as one

# group by distance bin, visit, species,
bat_detections_grpd = bat_detections %>%
  group_by(Location, Distance_bin, Month) %>%
  summarise(total_reads = sum(read_count), total_positives = sum(read_count > 0)) %>%
  ungroup()

bat_detections_grpd$Distance_bin = factor(bat_detections_grpd$Distance_bin, levels = c('0', '50', '100', '150'))

# plot read count in each distance bin
ggplot(bat_detections_grpd[bat_detections_grpd$Location=='Canada Farm',],
       aes(x = Distance_bin, y = log(total_reads), fill = Month)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = 'single')) +
  # geom_line(aes(group = interaction(SpeciesID, Visit)), alpha = 0.3) +
  scale_y_continuous(limits = c(0,15), expand = c(0,0)) +
  theme_classic() +
  scale_fill_manual(values = month_cols, breaks = levels(bat_detections_grpd$Month)) +
  labs(title = "All Bats - Canada Farm", x = "Distance (m)", y = "Log(Total Reads)", fill = "Visit") +
  # only include labels where total_positives > 0
  theme(text = element_text(size = 18,),
        strip.background = element_blank(),
        strip.text = element_text(size = 18))



# try out plotting using point ID instead of binned distances
bat_dets_point = bat_detections %>%
  group_by(Location, Point_ID, Visit, SpeciesID) %>%
  summarise(total_reads = sum(read_count), total_positives = sum(read_count > 0)) %>%
  ungroup()

# plot read count in each distance bin
ggplot(bat_dets_point[bat_dets_point$Location=='Canada Farm',], aes(x = Point_ID, y = total_reads, fill = as.factor(Visit))) +
  geom_col(position = position_dodge2(width = 0.9, preserve = 'single')) +
  # geom_line(aes(group = interaction(SpeciesID, Visit)), alpha = 0.3) +
  facet_wrap(~SpeciesID, scales = 'free') +
  ylim(0,220000) +
  theme_classic() +
  labs(title = "Bat Detections - Richmond Park", x = "Distance (m)", y = "Total Reads", fill = "Visit") +
  # only include labels where total_positives > 0
  theme(text = element_text(size = 18,),
        strip.background = element_blank(),
        strip.text = element_text(size = 18))

# fix point name - places called A need to be sorted
ggplot(bat_dets_point[bat_dets_point$Location=='Canada Farm',], aes(x = Point_ID, y = SpeciesID, fill = total_reads)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  facet_wrap(~Visit, nrow = 1) +
  labs(title = "Bat Detections at Canada Farm", x = "Site ID", y = "Species ID", fill = "Total Reads")


ggplot(bat_detections[bat_detections$Location=='Richmond Park',], aes(x = Point_ID, y = SpeciesID, fill = read_count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  facet_wrap(~Visit) +
  labs(title = "Bat Detections at Richmond Park", x = "Site ID", y = "Species ID", fill = "Total Reads")


#

# Maps with species detections --------------------------------------------
library(sf)
library(leaflet)
library(terra)
library(ggplot2)
library(tidyterra)
# POA: make basemap which has land cover map and all ppints in grey
# for each species, overlay points where detected, sized by number of detections (total possible = 9, 3 nights * 3 reps)
# keep unique combos of pointID and location
siteinfo_un = siteinfo %>%
  select(Point_ID, Location, Latitude, Longitude) %>%
  distinct()
# make into sf object
cfpoints = siteinfo_un %>%
  filter(Location == 'Canada Farm') %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
cfpoints_proj = st_transform(cfpoints, crs = 27700)


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
# crop LCM
ext = ext(cfpoints_proj) + 10
LCM_cf = crop(LCM, ext)
plot(LCM_cf$label, col = terrain.colors(length(lut$label)))
plot(cfpoints_proj["Point_ID"], add = TRUE, col = 'grey', pch = 19)

# colours:


mapcols = c('grey70',"#1b7837", "#1b7837","#e7d4e8", "#af8dc3", "#762a83", "#f7f7f7", "#7fbf7b","#d9f0d3", "#543005", "#74add1", "#d73027", "#fc8d59")
names(mapcols) = lut$label

# save this map as basemap
basemap = ggplot() +
  geom_spatraster(data = LCM_cf, aes(fill = label)) +
  scale_fill_manual(values = mapcols, na.value = 'transparent', name = 'Land Cover')

basemap

# now we want detections for each bat species
# dataframe with locattion, point id, visit, species and total positive samples

# first need to summarise results by species ID and sum read count
sample_sum = bat_detections %>%
  group_by(Location, Point_ID, Month, Night, pcr_replicate, SpeciesID) %>%
  summarise(read_count = sum(read_count))


sample_sum2 = sample_sum %>%
  group_by(Location, Point_ID, Month, SpeciesID) %>%
  summarise(total_positives = sum(read_count > 0)) %>%
  ungroup()

# add point coordinates
sample_sum3 = sample_sum2 %>%
  left_join(siteinfo_un, by = c("Location", "Point_ID"))
sample_sum3 = st_as_sf(sample_sum3, coords = c("Longitude", "Latitude"), crs = 4326)
# convert to uk map
sample_sum3 = st_transform(sample_sum3, crs = 27700)
# plot
basemap +
  geom_sf(data = sample_sum3[sample_sum3$Location=='Canada Farm' & sample_sum3$SpeciesID=='Pipistrellus spp.',],
          aes(size = total_positives), color = 'black', shape = 21, fill = 'yellow', alpha = 0.7) +
  facet_wrap(~Month) +
  ggtitle("Pipistrellus spp. Detections at Canada Farm") +
  coord_sf(expand = F) +
  theme(strip.text = element_text(size = 16),
        text = element_text(size = 14))

basemap +
  geom_sf(data = sample_sum3[sample_sum3$Location=='Canada Farm' & sample_sum3$SpeciesID=='Barbastella barbastellus',],
          aes(size = total_positives), color = 'black', shape = 21, fill = 'yellow', alpha = 0.7) +
  facet_wrap(~Month) +
  ggtitle("Barabastella Detections at Canada Farm") +
  coord_sf(expand = F) +
  theme(strip.text = element_text(size = 16),
        text = element_text(size = 14))


basemap +
  geom_sf(data = sample_sum3[sample_sum3$Location=='Canada Farm' & sample_sum3$SpeciesID=='Eptesicus serotinus',],
          aes(size = total_positives), color = 'black', shape = 21, fill = 'yellow', alpha = 0.7) +
  facet_wrap(~Month) +
  ggtitle("Eptesicus serotinus Detections at Canada Farm") +
  coord_sf(expand = F) +
  theme(strip.text = element_text(size = 16),
        text = element_text(size = 14))

basemap +
  geom_sf(data = sample_sum3[sample_sum3$Location=='Canada Farm' & sample_sum3$SpeciesID=='Myotis spp.',],
          aes(size = total_positives), color = 'black', shape = 21, fill = 'yellow', alpha = 0.7) +
  facet_wrap(~Month) +
  ggtitle("Myotis spp. Detections at Canada Farm") +
  coord_sf(expand = F) +
  theme(strip.text = element_text(size = 16),
        text = element_text(size = 14))

basemap +
  geom_sf(data = sample_sum3[sample_sum3$Location=='Canada Farm' & sample_sum3$SpeciesID=='Rhinolophus',],
          aes(size = total_positives), color = 'black', shape = 21, fill = 'yellow', alpha = 0.7) +
  facet_wrap(~Month) +
  ggtitle("Rhinolophus Detections at Canada Farm") +
  coord_sf(expand = F) +
  theme(strip.text = element_text(size = 16),
        text = element_text(size = 14))


# Richmond park
# Richmond park landcover
# UK CEH Land Cover Map 2024 for canada farm only
LCMRP = rast('Data/CEH_LCM_2024_RP/data/LCM.tif')
plot(LCMRP)
unique_LCs = unique(LCMRP$LCM_1)

levels(LCMRP) = lut
plot(LCMRP$label, col = terrain.colors(length(lut$label)))

# make into sf object
rppoints = siteinfo_un %>%
  filter(Location == 'Richmond Park') %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
rppoints_proj = st_transform(rppoints, crs = 27700)

# crop LCM
ext = ext(rppoints_proj[rppoints_proj$Location=='Richmond Park',]) + 10
LCM_rp = crop(LCMRP, ext)
plot(LCM_rp$label, col = terrain.colors(length(lut$label)))
plot(rppoints_proj["Point_ID"], add = TRUE, col = 'blue', pch = 20)

# save this map as basemap
basemaprp = ggplot() +
  geom_spatraster(data = LCM_rp, aes(fill = label)) +
  scale_fill_manual(values = mapcols, na.value = 'transparent', name = 'Land Cover')

basemaprp

# now we want detections for each bat species
# dataframe with locattion, point id, visit, species and total positive samples

# plot
basemaprp +
  geom_sf(data = sample_sum3[sample_sum3$Location=='Richmond Park' & sample_sum3$SpeciesID=='Myotis spp.',],
          aes(size = total_positives), color = 'black', shape = 21, fill = 'yellow', alpha = 0.7) +
  facet_wrap(~Month) +
  ggtitle("Myotis Detections at Richmond Park") +
  coord_sf(expand = F) +
  theme(strip.text = element_text(size = 16),
        text = element_text(size = 14))

basemaprp +
  geom_sf(data = sample_sum3[sample_sum3$Location=='Richmond Park' & sample_sum3$SpeciesID=='Nyctalus',],
          aes(size = total_positives), color = 'black', shape = 21, fill = 'yellow', alpha = 0.7) +
  facet_wrap(~Month) +
  ggtitle("Nyctalus Detections at Richmond Park") +
  coord_sf(expand = F) +
  theme(strip.text = element_text(size = 16),
        text = element_text(size = 14))


basemaprp +
  geom_sf(data = sample_sum3[sample_sum3$Location=='Richmond Park' & sample_sum3$SpeciesID=='Pipistrellus spp.',],
          aes(size = total_positives), color = 'black', shape = 21, fill = 'yellow', alpha = 0.7) +
  facet_wrap(~Month) +
  ggtitle("Pipistrellus spp. Detections at Richmond Park") +
  coord_sf(expand = F) +
  theme(strip.text = element_text(size = 16),
        text = element_text(size = 14))


basemaprp +
  geom_sf(data = sample_sum3[sample_sum3$Location=='Richmond Park' & sample_sum3$SpeciesID=='Plecotus',],
          aes(size = total_positives), color = 'black', shape = 21, fill = 'yellow', alpha = 0.7) +
  facet_wrap(~Month) +
  ggtitle("Plecotus Detections at Richmond Park") +
  coord_sf(expand = F) +
  theme(strip.text = element_text(size = 16),
        text = element_text(size = 14))


# Acoustic data -----------------------------------------------------------
ad = read.csv("Data/NatureAir_acoustics_ARD.csv")
head(ad)


# the measure i want to use is total number of active ninutes divided by
# total number of minutes 15mins before sunset and 15 mins after sunrise

# this data has rows for each location, and a record for every minute of recording & every species.
# from 5pm - 8am every day
# add visit and day
ad$minutes = ifelse(nchar(ad$minutes) ==10,
                     paste0(ad$minutes, " 00:00:00"),
                            ad$minutes)

ad$minutes = ymd_hms(ad$minutes, tz = "GMT")
ad$Date = date(ad$minutes)
ad$Time = format(ad$minutes, format = "%H:%M:%S")
ad$Location = ifelse(grepl('CanadaFarm', ad$sample_location_ID), 'Canada Farm', 'Richmond Park')
ad$Point_ID = sapply(strsplit(as.character(ad$sample_location_ID), "_"), `[`, 2)
# Canada farm - 1, March, 2, May, 3 June, 4 August, 5 October
# Richmond park - 1, July, 2 August, 3 September, 4 October.
unique(ad$Date[grep('RichmondPark', ad$sample_location_ID)])

ad$Month = format(ad$Date, "%B")
# dates 30th sept should be relabelled as october.
# 1st august should be relabeled as july
ad$Month = ifelse(ad$Location == 'Richmond Park' & ad$Date == '2025-09-30', 'October', ad$Month)
ad$Month = ifelse(ad$Location == 'Canada Farm' & ad$Date == '2025-08-01', 'July', ad$Month)


# need to add a night and hour column
ad = ad %>%
  mutate(Hour = as.numeric(format(minutes, "%H")),
         Night = case_when(
           Hour >= 16 ~ as.Date(Date),
           Hour < 9 ~ as.Date(Date) - 1,
           TRUE ~ as.Date(NA)
         ))

# number of minute intervals in which species was detected per 60 minutes of survey time.
# I have gone on timeanddate.com and taken sunrise and sunset times for middle day of each visit.
# March 4th-6th CF 06:46 - 17:59
# May 6-8 CF 05:33-20:43
# June 10-12 CF 04:56 - 21:25
# July 1-3 RP 04:49 - 21:20
# July 29-30 CF 05:33 - 21:00
# August 5-7 RP 05:33 - 20:39
# September 2-3 RP 06:17 - 19:42
# September 30 - October 2 RP 07:02 - 18:38
# October 7-9 CF 07:23 - 18:33


# remove detections that are more than 15 minutes before sunset or after sunrise
win_tbl <- tribble(
  ~Location,         ~Month,     ~start_str, ~end_str,
  "Canada Farm",     "March",    "17:44:00", "07:01:00",
  "Canada Farm",     "May",      "20:28:00", "05:48:00",
  "Canada Farm",     "June",     "21:10:00", "04:41:00",
  "Canada Farm",     "July",     "20:45:00", "05:18:00",
  "Canada Farm",     "October",  "18:18:00", "07:38:00",
  "Richmond Park",   "July",     "21:05:00", "05:04:00",
  "Richmond Park",   "August",   "20:24:00", "05:48:00",
  "Richmond Park",   "September","19:27:00", "06:32:00",
  "Richmond Park",   "October",  "18:23:00", "07:17:00"
) %>%
  mutate(
    start_t = hms::as_hms(start_str),
    end_t   = hms::as_hms(end_str)
  ) %>%
  select(-start_str, -end_str)

# set time as same format as start and end times
ad$Time = hms::as_hms(ad$Time)

# filter out detetions that our not within our nighttime window
ad2 <- ad %>%
  inner_join(win_tbl, by = c("Location", "Month")) %>%
  filter(Time >= start_t | Time <= end_t)  # <-- across-midnight window
  select(-start_t, -end_t)

ad_dets = ad2 %>%
  group_by(Location, Point_ID, Month, species, Night) %>%
  summarise(total_pulses = sum(no_pulses, na.rm = TRUE),
            total_active_minutes = sum(no_pulses > 0, na.rm = TRUE),
            total_mintutes = n()) %>%
  ungroup()

ad_dets$activity = ad_dets$total_active_minutes / ad_dets$total_mintutes

# save
write.csv(ad_dets, "output/acoustic_activity_by_night.csv", row.names = FALSE)

ad_dets_month = ad_dets %>%
  group_by(Location, Point_ID, Month, species) %>%
  summarise(avg_activity = mean(activity),
            sd_activity = sd(activity)) %>%
  ungroup()

unique(ad_dets_month$species)
# canada farm detected species - myotis, pipistrellus, plecotus, rhinolophus, eptesicus, barbastella
# richmond park - myotis, pipistrellus, nyctalus, plecotus

# for each species, i want a plot that shoes activity at each site, and overlayed the airDNA detections
# add land use type to each point

siteinfo_un = siteinfo %>%
  select(Point_ID, Location, Latitude, Longitude, DominantLC_10m, dist_m) %>%
  distinct()

# add binned distance
siteinfo_un = siteinfo_un %>%
  mutate(Distance_bin = case_when(
  dist_m < 20 ~ '0',
  dist_m >= 20 & dist_m < 60 ~ '50',
  dist_m >= 60 & dist_m < 110 ~ '100',
  dist_m >= 110  ~ '150',
))
siteinfo_un$Distance_bin = factor(siteinfo_un$Distance_bin, levels = c('0', '50', '100', '150'))

# acoustic A sites were labelled differently to eDNA - fix here.
ad_dets_month = ad_dets_month %>%
  mutate(Point_ID2 = case_when(
    Location == 'Canada Farm' & Point_ID == 'AA' ~ "AS",
    Location == 'Canada Farm' & Point_ID == 'AB' ~ "AN",
    Location == 'Canada Farm' & Point_ID == 'AC' ~ "AO",
    Location == 'Richmond Park' & Point_ID == 'AA' ~ "AN",
    Location == 'Richmond Park' & Point_ID == 'AB' ~ "AS",
    Location == 'Richmond Park' & Point_ID == 'AC' ~ "AO",
    .default = Point_ID))

adm2 = ad_dets_month %>%
  left_join(siteinfo_un, by = c("Location", "Point_ID2" = "Point_ID"))

adm2 = st_as_sf(adm2, coords = c("Longitude", "Latitude"), crs = 4326)
# convert to uk map
adm2 = st_transform(adm2, crs = 27700)

sample_sum4 = sample_sum3[sample_sum3$total_positives > 0, ]

adm2$Month = factor(adm2$Month, levels = c('March', 'May', 'June', 'July', 'August', 'September', 'October'))

unique(adm2$species[adm2$Location=='Canada Farm' & adm2$avg_activity > 0])
unique(adm2$species[adm2$Location=='Richmond Park' & adm2$avg_activity > 0]))
# we only want to keep in spcies we are interested in
# canada farm detected species - myotis, pipistrellus, plecotus, rhinolophus, eptesicus, barbastella
# richmond park - myotis, pipistrellus, nyctalus, plecotus
adm3 = adm2 %>%
  filter((Location == 'Canada Farm' &
            species %in% c("Barbastella barbastellus",
                           "Myotis", "Nyctalus/Cnephaeus agg.", "Pipistrellus",
                           "Plecotus", "Rhinolophus")) |
           (Location == 'Richmond Park' &
              species %in% c("Myotis", "Nyctaus", "Nyctalus/Cnephaeus agg.", "Pipistrellus", "Plecotus")))

adm3 = adm3 %>%
  mutate(Habitat_type = case_when(DominantLC_10m %in% c('Broadleaf woodland') ~ 'Closed',
                                  DominantLC_10m %in% c('Improve grassland', 'Neutral grassland', 'Fen') ~ 'Open',
                                  TRUE ~ 'Other'))

# get average activity per distance bin
adm4 = adm3 %>%
  group_by(Location, Distance_bin, Month, species) %>%
  summarise(avg_activity = mean(avg_activity),
            sd_activity = sd(avg_activity)) %>%
  ungroup()

# get average activity across land-use types
adm5 = adm3 %>%
  group_by(Location, DominantLC_10m, Month, species) %>%
  summarise(avg_activity = mean(avg_activity),
            sd_activity = sd(avg_activity)) %>%
  ungroup()

# remove barbs for now as there are less than 3 eDNA detections
adm3 = adm3 %>%
  filter(species != 'Barbastella barbastellus')
# rename nyctalus.cnephaues to be serotine
adm3$species = ifelse(adm3$species == 'Nyctalus/Cnephaeus agg.', "Eptesicus sp.", adm3$species)

# Plot acoustic activity over distance ------------------------------------

# use smooth line
# Canada Farm
cf1 = ggplot(adm3[adm3$Location=='Canada Farm',],
       aes(x = dist_m, y = avg_activity, color = Month)) +
  geom_smooth() +
  # geom_line(aes(group = interaction(SpeciesID, Visit)), alpha = 0.3) +
  facet_wrap(~species, scales = 'free', nrow = 1) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme_classic() +
  labs(title = "Canada Farm", x = "Distance (m)", y = "Acoustic activity", color = "Visit") +
  scale_color_manual(values = month_cols, breaks = levels(bat_detections_sum$Month)) +
  # only include labels where total_positives > 0
  theme(text = element_text(size = 18,),
        strip.background = element_blank(),
        strip.text = element_text(size = 18))

# plot average activity over land use type
cf2 = ggplot(adm3[adm3$Location=='Canada Farm',],
       aes(x = Habitat_type, y = avg_activity, color = Month)) +
  geom_boxplot() +
  facet_wrap(~species, scales = 'free', nrow = 1) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme_classic() +
  labs(x = "Habitat type", y = "Acoustic activity", color = "Visit") +
  scale_color_manual(values = month_cols, breaks = levels(bat_detections_sum$Month)) +
  # only include labels where total_positives > 0
  theme(text = element_text(size = 18,),
        strip.background = element_blank(),
        strip.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1))

cf2
cf1+cf2 + plot_layout(nrow = 2)


# Richmond Park
rp1 = ggplot(adm3[adm3$Location=='Richmond Park' & adm3$species!='Myotis',],
             aes(x = dist_m, y = avg_activity, color = Month)) +
  geom_smooth() +
  # geom_line(aes(group = interaction(SpeciesID, Visit)), alpha = 0.3) +
  facet_wrap(~species, scales = 'free', nrow = 1) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme_classic() +
  labs(title = "Richmond Park", x = "Distance (m)", y = "Acoustic activity", color = "Visit") +
  scale_color_manual(values = month_cols, breaks = levels(bat_detections_sum$Month)) +
  # only include labels where total_positives > 0
  theme(text = element_text(size = 18,),
        strip.background = element_blank(),
        strip.text = element_text(size = 18))
rp1
# plot average activity over land use type
rp2 = ggplot(adm3[adm3$Location=='Richmond Park' & adm3$species!='Myotis',],
             aes(x = Habitat_type, y = avg_activity, color = Month)) +
  geom_boxplot() +
  facet_wrap(~species, scales = 'free', nrow = 1) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme_classic() +
  labs(x = "Habitat type", y = "Acoustic activity", color = "Visit") +
  scale_color_manual(values = month_cols, breaks = levels(bat_detections_sum$Month)) +
  # only include labels where total_positives > 0
  theme(text = element_text(size = 18,),
        strip.background = element_blank(),
        strip.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1))

rp2
rp1+rp2 + plot_layout(nrow = 2)

ggplot(adm4[adm4$Location=='Canada Farm',],
       aes(x = Distance_bin, y = avg_activity, fill = Month)) +
  geom_smooth() +
  # geom_col(position = position_dodge2(width = 0.9, preserve = 'single')) +
  # geom_line(aes(group = interaction(SpeciesID, Visit)), alpha = 0.3) +
  facet_wrap(~species, scales = 'free') +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  labs(title = "Canada Farm", x = "Distance (m)", y = "Acoustic activity", fill = "Visit") +
  scale_fill_manual(values = month_cols, breaks = levels(bat_detections_sum$Month)) +
  # only include labels where total_positives > 0
  theme(text = element_text(size = 18,),
        strip.background = element_blank(),
        strip.text = element_text(size = 18))


ggplot(adm5[adm5$Location=='Canada Farm',],
       aes(x = DominantLC_10m, y = avg_activity, fill = Month)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = 'single')) +
  # geom_line(aes(group = interaction(SpeciesID, Visit)), alpha = 0.3) +
  facet_wrap(~species, scales = 'free') +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  labs(title = "Canada Farm", x = "Land Cover", y = "Acoustic activity", fill = "Visit") +
  scale_fill_manual(values = month_cols, breaks = levels(bat_detections_sum$Month)) +
  # only include labels where total_positives > 0
  theme(text = element_text(size = 18,),
        strip.background = element_blank(),
        strip.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1))


# Richmond Park

ggplot(adm4[adm4$Location=='Richmond Park',],
       aes(x = Distance_bin, y = avg_activity, fill = Month)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = 'single')) +
  # geom_line(aes(group = interaction(SpeciesID, Visit)), alpha = 0.3) +
  facet_wrap(~species, scales = 'free') +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  labs(title = "Richmond Park", x = "Distance (m)", y = "Acoustic activity", fill = "Visit") +
  scale_fill_manual(values = month_cols, breaks = levels(bat_detections_sum$Month)) +
  # only include labels where total_positives > 0
  theme(text = element_text(size = 18,),
        strip.background = element_blank(),
        strip.text = element_text(size = 18))


ggplot(adm5[adm5$Location=='Richmond Park',],
       aes(x = DominantLC_10m, y = avg_activity, fill = Month)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = 'single')) +
  # geom_line(aes(group = interaction(SpeciesID, Visit)), alpha = 0.3) +
  facet_wrap(~species, scales = 'free') +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  scale_fill_manual(values = month_cols, breaks = levels(bat_detections_sum$Month)) +
  # only include labels where total_positives > 0
  theme(text = element_text(size = 18,),
        strip.background = element_blank(),
        strip.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, hjust = 1))


# Plot acoustic acitivity on maps -----------------------------------------


library(ggnewscale)

basemap +
  geom_sf(data = adm2[adm2$Location=='Canada Farm' & adm2$species=='Barbastella barbastellus',],
          aes(size = log(avg_activity),
              colour = avg_activity == 0),
          shape = 22,
          fill = 'lightblue', alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = 'grey', "FALSE" = 'blue'), guide = 'none') +
  theme(strip.text = element_text(size = 16),
        text = element_text(size = 14)) +
  geom_sf(data = sample_sum4[sample_sum4$Location=='Canada Farm' &
                               sample_sum4$SpeciesID=='Barbastella barbastellus',],
          aes(size = total_positives),
          color = 'black', shape = 21, fill = 'yellow', alpha = 0.7) +
  coord_sf(expand = F) +
  facet_wrap(~Month) +
  ggtitle("Canada Farm, Barbastella barbastellus")


# ella's thoughts on bat data
# pulses can be overinflated by one bat
# has to be more than 3 pulses per minute to be counted as detection
# if you only have one detection, could be uncertain. its quite unlikely that you would just hear one pulse.
# number of minutes per night / total minutes per species.
# called occupancy activity measure

