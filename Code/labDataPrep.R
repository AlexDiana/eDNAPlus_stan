# SPATIAL DATA PREP
library(dplyr)
library(tidyr)
library(ggplot2)

Sd = read.csv("Data/NatureAir_SpatialDATA.csv")
siteinfo = read.csv("output/allsitevars.csv")

# fix typo
# when
siteinfo$SampleID[siteinfo$SampleID == "CFANN3V1" & siteinfo$Visit == 2] = "CFANN3V2"
siteinfo$SampleID[siteinfo$SampleID == "CFAON3V1" & siteinfo$Visit == 2] = "CFAON3V2"
siteinfo$SampleID[siteinfo$SampleID == "CFASN3V1" & siteinfo$Visit == 2] = "CFASN3V2"
# there are two rows where sampleID = RPAS1N1V1
# one needs to be RPASN1V1, the other RPANN1V1
siteinfo[163,]$SampleID = "RPASN1V1"
siteinfo[164,]$SampleID = "RPANN1V1"
siteinfo[164,]$Point_ID = "AN"


# convert Sd to longform
sd_long = pivot_longer(Sd, cols = 15:ncol(Sd), names_to = 'Sample', values_to = 'read_count')


# remvoe string from Sample column containing 'rep' followed by a number. Add the number to a new column called 'replicate'
sd_long = sd_long %>%
  mutate(pcr_replicate = as.numeric(gsub(".*rep(\\d+).*", "\\1", Sample))) %>%
  mutate(SampleID = gsub("rep\\d+", "", Sample))

# to check - there are some mock community results that to not have PCR rep numbers. ignoring for now, these are NA

# correct typo
sd_long$SampleID[grep("CFAAN1V2", sd_long$SampleID)] = "CFANN1V2"
sd_long$SampleID[grep("RPAN1N1V1", sd_long$SampleID)] = "RPANN1V1"
sd_long$SampleID[grep("RPAN1N1V2", sd_long$SampleID)] = "RPANN1V2"
sd_long$SampleID[grep("RPAN1N1V3", sd_long$SampleID)] = "RPANN1V3"
sd_long$SampleID[grep("RPAS1N1V1", sd_long$SampleID)] = "RPASN1V1"
sd_long$SampleID[grep("RPAS1N1V2", sd_long$SampleID)] = "RPASN1V2"
sd_long$SampleID[grep("RPAS1N1V3", sd_long$SampleID)] = "RPASN1V3"

# when Sample = "RPIN1V3rep1.1", change PCR rep to 3 and sample id to "RPIN1V3"
sd_long$SampleID[sd_long$Sample == "RPIN1V3rep1.1"] = "RPIN1V3"
sd_long$pcr_replicate[sd_long$Sample == "RPIN1V3rep1.1"] = 3

# join with metadata from siteinfo
sd_long1 = left_join(sd_long, siteinfo, by = c("SampleID" = "SampleID"), relationship = "many-to-many")


# check which samples did not line up.
mismatch = sd_long1 %>%
  filter(is.na(Point_ID)) %>%
  filter(!grepl("Mock", SampleID)) %>% #not interested in mock community
  filter(!grepl("EBA|EBT", SampleID)) # not interested in extraction blanks

unique(mismatch$SampleID)
# CFAON3V2

# something is missing V5.
# do all the rows in siteinfo match up with something?
siteinfo_ids = unique(siteinfo$SampleID)
sd_ids = unique(sd_long$SampleID)
missing_ids = siteinfo_ids[!siteinfo_ids %in% sd_ids]
# there are 53 samples that aren't in the spatial data. assuming this is because the samples failed somehow.
# i can't find anything that matches with 'EB' - no idea what this is.
# removing for now later on will see if something is missing V5


# remove any rows that didn't match up to site info - these are mainly the EB sample

sd_long1 = filter(sd_long1, !is.na(Point_ID))

# summarise overall detections across PCR replicates and days
sd_detections = sd_long1 %>%
  filter(SampleType == 'Active') %>%
  group_by(Family,Genus,Species, Location, Point_ID, Visit) %>%
  summarise(total_reads = sum(read_count)) %>%
  ungroup() %>%
  filter(Family!="")

# fix empty species/genus names
sd_detections$Species[sd_detections$Species == ""] = NA
sd_detections$Genus[sd_detections$Genus == ""] = NA

sd_detections$SpeciesID = ifelse(is.na(sd_detections$Species), sd_detections$Genus, sd_detections$Species)
sd_detections$SpeciesID = ifelse(is.na(sd_detections$Genus), sd_detections$Family, sd_detections$SpeciesID)

ggplot(sd_detections[sd_detections$Location=='Canada Farm',], aes(x = Point_ID, y = SpeciesID, fill = total_reads)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  facet_wrap(~Visit, nrow = 1) +
  labs(title = "All Detections at Canada Farm", x = "Site ID", y = "Species ID", fill = "Total Reads")


# what is the mean read count of humans and cows?
mean(sd_detections$total_reads[sd_detections$Species == "Homo sapiens"], na.rm = TRUE)
mean(sd_detections$total_reads[sd_detections$Species == "Bos taurus"], na.rm = TRUE)
mean(sd_detections$total_reads[sd_detections$SpeciesID == "Plecotus"], na.rm = TRUE)
mean(sd_detections$total_reads[sd_detections$Species == "Pipistrellus pipistrellus"], na.rm = TRUE)
mean(sd_detections$total_reads[sd_detections$Species == "Pipistrellus pipistrellus"], na.rm = TRUE)

mean(sd_detections$total_reads[sd_detections$SpeciesID %in% c('Homo sapiens', "Bos taurus")], na.rm = TRUE)
mean(sd_detections$total_reads[!sd_detections$SpeciesID %in% c('Homo sapiens', "Bos taurus")], na.rm = TRUE)

# just bats
sd_bats = sd_long1 %>%
  filter(Order == "Chiroptera")
unique(sd_bats$Species)

# get detection frequencies for each species. sum read_count over visit
bat_detections = sd_bats %>%
  filter(SampleType == 'Active') %>%
  group_by(Genus,Species, Location, Point_ID, Visit) %>%
  summarise(total_reads = sum(read_count)) %>%
  ungroup()

bat_detections$Species[bat_detections$Species == ""] = NA
bat_detections$SpeciesID = ifelse(is.na(bat_detections$Species), bat_detections$Genus, bat_detections$Species)

# remove mock communities and controls
bat_detections = bat_detections %>%
  filter(!is.na(Point_ID))

# fix point name - places called A need to be sorted


ggplot(bat_detections[bat_detections$Location=='Canada Farm',], aes(x = Point_ID, y = SpeciesID, fill = total_reads)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  facet_wrap(~Visit, nrow = 1) +
  labs(title = "Bat Detections at Canada Farm", x = "Site ID", y = "Species ID", fill = "Total Reads")


ggplot(bat_detections[bat_detections$Location=='Richmond Park',], aes(x = Point_ID, y = SpeciesID, fill = total_reads)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  facet_wrap(~Visit) +
  labs(title = "Bat Detections at Richmond Park", x = "Site ID", y = "Species ID", fill = "Total Reads")

# Combine detections with quant data ----------------------------------------


sq = read.csv('Data/NatureAir_SpatialQuant.csv')
head(sq)
# at the moment ignoring mock, EBA, and EBT
sq1 = sq %>%
  filter(!grepl("Mock", MiseqRunID)) %>%
  filter(!grepl("EBA|EBT", MiseqRunID))

# fixing a typo - again not sure if this Site I is missing a rep
sd_long1$Sample[sd_long1$Sample == "RPIN1V3rep1.1"] = "RPIN1V3rep3"
sq1[991,]$MiseqRunID = "RPIN1V3rep3"

# we need the pre-index concentration
# start by seeing what matches up
sq_id = unique(sq1$MiseqRunID)
sd_id = unique(sd_long1$Sample)
# which sq ids are not in sd_long?
missing_sq = sd_id[!sd_id %in% sq_id] #all the ids in sd_long, match with something in sq
missing_sd = sq_id[!sq_id %in% sd_id] #there are 31 samples that were quanted but not metabarcoding results
missing_ids #there are 53 samples in siteinfo that are not in sd_long

# there are lots of missing samples - assuming this is because nothing is detected.
# will need to include these as blanks at some point.

# join sq1 with sd_long1 by Sample and MiseqRunID
sd_long2 = left_join(sd_long1, sq1, by = c("Sample" = "MiseqRunID"), relationship = "many-to-many")

sd_long2$PreIndexConcentration.nM. = as.numeric(sd_long2$PreIndexConcentration.nM.)
sd_long2$PreIndexConcentration.ng.ul. = as.numeric(sd_long2$PreIndexConcentration.ng.ul.)
sd_long2$PostIndexConcentration.nM. = as.numeric(sd_long2$PostIndexConcentration.nM.)
sd_long2$PostIndexConcentration.ng.ul. = as.numeric(sd_long2$PostIndexConcentration.ng.ul.)

# save this
write.csv(sd_long2, "output/spatial_data_longform_withquant.csv", row.names = FALSE)


# the other thing to do is combine species.
# there are 219 species level assignments, but many are likely to be the same species.
species_list = unique(sd_long2[,c('ASV_ID',"Class", "Order",'Family','Genus','Species', "Status")])

# lets group some of these together.
species_list$SpeciesID = ifelse(species_list$Species != "", species_list$Species, species_list$Genus)
species_list$SpeciesID = ifelse(species_list$SpeciesID != "", species_list$SpeciesID, species_list$Family)

species_list$SpeciesID[species_list$Status == 'Contaminant'] = "Contaminant"
species_list$SpeciesID[species_list$Status== "Unassigned"] = "Unassigned"

species_list$SpeciesID[species_list$Genus== "Meles"] = "Meles meles" #badger
species_list$SpeciesID[species_list$Order== "Primates" & species_list$Species==""] = "Contaminant" #humans
species_list$SpeciesID[species_list$Genus== "Microtus"] = "Microtus agrestis" #vole
species_list$SpeciesID[species_list$Genus== "Myodes"] = "Myodes glareolus" #vole
species_list$SpeciesID[species_list$Family== "Sciuridae" & species_list$Species==""] = "Sciurus"
species_list$SpeciesID[species_list$Class== "Actinopterygii"] = "Fish"
species_list$SpeciesID[species_list$Genus== "Pica"] = "Pica pica" #magpie
species_list$SpeciesID[species_list$Class== "Aves" & species_list$Family==""] = "Unassigned bird"

# all birds that could not be assigned below family can be unassigned bird
species_list$SpeciesID[species_list$Class== "Aves" & species_list$Genus==""] = "Unassigned bird"

# now lets group some species together to genus level
# Apodemus
species_list$SpeciesID[species_list$Genus== "Apodemus"] = "Apodemus spp."
# Sciurus
species_list$SpeciesID[species_list$Family== "Sciuridae"] = "Sciurus spp."
# bufo
species_list$SpeciesID[species_list$Genus== "Bufo"] = "Bufo spp."

# assign NUMT as NUMT
species_list$SpeciesID[species_list$Status== "NUMT"] = "NUMT"
species_list$SpeciesID[grepl("Synthetic", species_list$Status)] = "Synthetic"

unique(species_list$SpeciesID)

# there are 74 species. another thing we could do is group all birds together




# for now, join species_ID back up with sd_long2
sd_long3 = left_join(sd_long2, species_list[,c('ASV_ID','SpeciesID')], by = c("ASV_ID" = "ASV_ID"), relationship = "many-to-many")


# we also need to add in the missing samples. these are samples that are in siteinfo but not in sd_long3.
# for every species, need to give a 0 value.
# most of missing id is richmond park V4 which we don't have results for yet.
missing_ids1 = missing_ids[grepl('CF', missing_ids)] # canada farm
# convert to include PCR reps
missing_ids2 = c(paste0(missing_ids1, "rep1"), paste0(missing_ids1, "rep2"), paste0(missing_ids1, "rep1"))
missing_sd
# remove extraction blanks
missing_sd = missing_sd[!grepl("CFEB", missing_sd)]

missing_sd = c(missing_ids2, missing_sd)
missing_sd = missing_sd[missing_sd!=""]

missing_sample = unique(sd_long2[,1:14])

names(sd_long3)
missing_sq = sq[1,]
missing_sq[,4:7] = 0


sd_long4 = sd_long3
for(i in 1:length(missing_sd)){
  # remove rep
  sampleid = gsub("rep\\d+", "", missing_sd[i])
  missing_sample_temp = missing_sample %>%
    mutate(Sample = missing_sd[i],
           read_count = 0,
           pcr_replicate = as.numeric(gsub(".*rep(\\d+).*", "\\1", missing_sd[i])),
           SampleID = sampleid,)
  temp = siteinfo %>%
    filter(SampleID == sampleid)

  temp2 = left_join(missing_sample_temp, temp, by = 'SampleID')

  # add sq data
  sq_temp = sq1 %>%
    filter(MiseqRunID == missing_ids[i])
  if (nrow(sq_temp) > 0) {
    temp3 = left_join(temp2, sq_temp, by = c("Sample" = "MiseqRunID"))
  } else {
    sq_temp = missing_sq
    sq_temp$MiseqRunID = missing_sd[i]
    sq_temp$Sample.ID = sampleid
    temp3 = left_join(temp2, sq_temp, by = c("Sample" = "MiseqRunID"))
  }
  temp3 = left_join(temp3, species_list[,c('ASV_ID','SpeciesID')], by = c("ASV_ID" = "ASV_ID"))
  sd_long4 = rbind(sd_long4, temp3)
}

siteinfo_ids = unique(siteinfo$SampleID)
sd_ids = unique(sd_long4$SampleID)
sq_ids = unique(sq$MiseqRunID)
sdsq_ids = unique(sd_long4$Sample)
# which sqsq are missing from sq_ids?
missing = sdsq_ids[!sdsq_ids %in% sq_ids]


write.csv(sd_long4, "output/spatial_data_longform_withquant_speciesID.csv", row.names = FALSE)
write.csv(siteinfo, "output/all_site_vars2.csv", row.names = FALSE)


