# SPATIAL DATA PREP
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# explanation of data files.
# as of March 2026, the most up to date file is NatureAirSpatialCombined. This
# includes data from all visits, including RP V4 which was run in a later sequencing batch
# so was not in the original data.
# it is also different to NatureAir_16SmamSpatial_ALL as it has the corrected sample
# that originally suggested there was a horseshoe bat in RP, which was due to a
# labelling mistake

# 'Spatial1' and 'Spatial2' represent data from the first sequencing batch and
# the second, respectively.
# Spatial2 only includes results from RP V4.
# SpatialCombined is the combination of these two, and should be the most up to
# date file.
# sd = read.csv("Data/NatureAir_16SmamSpatial_ALL.csv")
Sd = read.delim2('Data/NatureAirSpatialCombined.txt')
# sd1 = read.delim2('Data/NatureAirSpatial1.txt')
# sd2 = read.delim2('Data/NatureAirSpatial2.txt')

# difference in row numbers - there are 251 ASVs x 1506-15 samples in OG
# there are 255 ASVs x 1496-15 samples in new.
# need to work out which samples are different

siteinfo = read.csv("output/allsitevars.csv")

# fix typos
siteinfo$SampleID[siteinfo$SampleID == "CFANN3V1" & siteinfo$Visit == 2] = "CFANN3V2"
siteinfo$SampleID[siteinfo$SampleID == "CFAON3V1" & siteinfo$Visit == 2] = "CFAON3V2"
siteinfo$SampleID[siteinfo$SampleID == "CFASN3V1" & siteinfo$Visit == 2] = "CFASN3V2"
# there are two rows where sampleID = RPAS1N1V1
# one needs to be RPASN1V1, the other RPANN1V1
siteinfo[163,]$SampleID = "RPASN1V1"
siteinfo[164,]$SampleID = "RPANN1V1"
siteinfo[164,]$Point_ID = "AN"

# there are some weird points in siteinfo
# CFAN1V5, CFAN2V5, CFAN3V5
# CFAO is all accounted for, CFAS has nothing on visit 5, CFAN has no visit 5
# relabelling these as AN for now.
siteinfo$SampleID[siteinfo$SampleID == "CFAN1V5"] = "CFANN1V5"
siteinfo$SampleID[siteinfo$SampleID == "CFAN2V5"] = "CFANN2V5"
siteinfo$SampleID[siteinfo$SampleID == "CFAN3V5"] = "CFANN3V5"

# convert Sd to longform
sd_long = pivot_longer(Sd, cols = 15:ncol(Sd), names_to = 'Sample', values_to = 'read_count')

# remove string from Sample column containing 'rep' followed by a number. Add the number to a new column called 'replicate'
sd_long = sd_long %>%
  mutate(pcr_replicate = as.numeric(gsub(".*rep(\\d+).*", "\\1", Sample))) %>%
  mutate(SampleID = gsub("rep\\d+", "", Sample))
# warning that NAs are introduced is ok
# to check - there are some mock community results that to not have PCR rep numbers. ignoring for now, these are NA

# correct typos
sd_long$SampleID[grep("CFAAN1V2", sd_long$SampleID)] = "CFANN1V2"
sd_long$SampleID[grep("CFAN1V5", sd_long$SampleID)] = "CFANN1V5"
sd_long$SampleID[grep("CFAN2V5", sd_long$SampleID)] = "CFANN2V5"
sd_long$SampleID[grep("CFAN3V5", sd_long$SampleID)] = "CFANN3V5"
sd_long$SampleID[grep("RPAN1N1V1", sd_long$SampleID)] = "RPANN1V1"
sd_long$SampleID[grep("RPAN1N1V2", sd_long$SampleID)] = "RPANN1V2"
sd_long$SampleID[grep("RPAN1N1V3", sd_long$SampleID)] = "RPANN1V3"
sd_long$SampleID[grep("RPAS1N1V1", sd_long$SampleID)] = "RPASN1V1"
sd_long$SampleID[grep("RPAS1N1V2", sd_long$SampleID)] = "RPASN1V2"
sd_long$SampleID[grep("RPAS1N1V3", sd_long$SampleID)] = "RPASN1V3"

# when Sample = "RPIN1V3rep1.1", change PCR rep to 3 and sample id to "RPIN1V3"
sd_long$SampleID[sd_long$Sample == "RPIN1V3rep1.1"] = "RPIN1V3"
sd_long$pcr_replicate[sd_long$Sample == "RPIN1V3rep1.1"] = 3

# in the new dataset (03/2026), there is no CFLN1V2 rep 3
# sd_long$SampleID[sd_long$Sample == "CFLN1V2rep2.1"] = "CFLN1V2"
# sd_long$pcr_replicate[sd_long$Sample == "CFLN1V2rep2.1"] = 3

# join with metadata from siteinfo
sd_long1 = left_join(sd_long, siteinfo, by = c("SampleID" = "SampleID"), relationship = "many-to-many")

# check which samples did not line up.
mismatch = sd_long1 %>%
  filter(is.na(Point_ID)) %>%
  filter(!grepl("Mock", SampleID)) %>% #not interested in mock community
  filter(!grepl("EBA|EBT", SampleID)) # not interested in extraction blanks

unique(mismatch$SampleID)
# do all the rows in siteinfo match up with something?
siteinfo_ids = unique(siteinfo$SampleID)
sd_ids = unique(sd_long$SampleID)
missing_ids = siteinfo_ids[!siteinfo_ids %in% sd_ids]
unmatches = sd_ids[!sd_ids %in% siteinfo_ids]
# these are samples that we collected but did not work some reason

# missing from the metabarcoding data
# CFASN1V2
# CFAON1v2
# RPKN2V4
# RPANN3V4
# "RPASN3V4" "RPPN3V4"
# CFASN1-3V5 OR CFANN1-3V5 - labelled as CFA
# EB is extraction blank

# remove any rows that didn't match up to site info - these are all EBs  and mock samples
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

ggplot(sd_detections[sd_detections$Location=='Canada Farm',], aes(x = Point_ID, y = Family, fill = total_reads)) +
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
# remove Myotis as there are no reads here
bat_detections = bat_detections %>%
  filter(!is.na(Point_ID)) %>%
  filter(SpeciesID != "Myotis") %>%
  droplevels()

bat_detections <- bat_detections %>%
  group_by(Location) %>%
  mutate(SpeciesID = factor(SpeciesID, levels = unique(SpeciesID[total_reads > 0]))) %>%
  ungroup()

bat_detections = bat_detections %>%
 filter(!is.na(SpeciesID)) %>%
  droplevels()

# bin read numbers to five groups
max(bat_detections$total_reads) #133107
# 0, 1-1000, 1001 - 10,000, 10,001 - 20,000, 20,0001 - 40,000, 40,000-140,000
bat_detections = bat_detections %>%
  mutate(read_bin = case_when(
    total_reads == 0 ~ '0',
    total_reads <= 1000 & total_reads > 0 ~ '1 - 1,000',
    total_reads > 1000 & total_reads <= 10000 ~ '1,001 - 10,000',
    total_reads > 10000 & total_reads <= 20000 ~ '10,001 - 20,000',
    total_reads > 20000 & total_reads <= 40000 ~ '20,001 - 40,000',
    total_reads > 40000 ~ '40,001 - 140,000',
    .default = NA),
    month = case_when(
      Location == 'Canada Farm' & Visit == 1 ~ "March",
      Location == 'Canada Farm' & Visit == 2 ~ "May",
      Location == 'Canada Farm' & Visit == 3 ~ "June",
      Location == 'Canada Farm' & Visit == 4 ~ "July",
      Location == 'Canada Farm' & Visit == 5 ~ "October",
      Location == 'Richmond Park' & Visit == 1 ~ "July",
      Location == 'Richmond Park' & Visit == 2 ~ "August",
      Location == 'Richmond Park' & Visit == 3 ~ "September",
      Location == 'Richmond Park' & Visit == 4 ~ "October",
    ),
    Point_ID = factor(Point_ID, levels = c("AN", "AS", "AO","B",
                                           "C",  "D",  "E",  "F",  "G",  "H",
                                           "I",  "J",  "K",  "L",  "M",  "N",
                                           "O",  "P",  "Q"))
    )

colour_pal = c("0" = 'white',
               '1 - 1,000' = "#440154FF",
               '1,001 - 10,000' = "#3B528BFF",
               '10,001 - 20,000' = "#21908CFF",
               '20,001 - 40,000' = "#5DC863FF",
               '40,001 - 140,000' = "#FDE725FF")
bat_detections$read_bin = factor(bat_detections$read_bin, levels = c("0", '1 - 1,000', '1,001 - 10,000', '10,001 - 20,000', '20,001 - 40,000', '40,001 - 140,000'))
bat_detections$month = factor(bat_detections$month, levels = c("March", "May", "June", "July", "August", "September", "October"))


cf_bats = unique(bat_detections$SpeciesID[bat_detections$Location == "Canada Farm" & bat_detections$total_reads > 0])
rp_bats = unique(bat_detections$SpeciesID[bat_detections$Location == "Richmond Park" & bat_detections$total_reads > 0])
bat_detections = bat_detections %>%
  mutate(cf_include = case_when(total_reads > 0 ~ TRUE,
                                Location == "Canada Farm" & SpeciesID %in% cf_bats ~ TRUE,
                                .default = FALSE),
         rp_include = case_when(total_reads > 0 ~ TRUE,
                                Location == "Richmond Park" & Point_ID == "AN" & Visit %in% c(3,4) & SpeciesID %in% rp_bats ~ TRUE,
                                .default = FALSE))
# fix point name - places called A need to be sorted
cfds = bat_detections %>% filter(bat_detections$Location=='Canada Farm' & bat_detections$cf_include==TRUE) %>%
  ggplot(aes(x = Point_ID, y = SpeciesID, fill = read_bin)) +
  geom_tile(show.legend = TRUE) +
  scale_fill_manual(values = colour_pal, breaks = names(colour_pal[-1]), drop = FALSE) +
  facet_wrap(~month, nrow = 1) +
  theme_minimal() +
  labs(title = "a)", x = "Site ID", y = "Species ID", fill = "Total Reads") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        text = element_text(size = 16),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
        ) +
  geom_hline(yintercept = 1.5, color = "grey", size = 0.5) +
  geom_hline(yintercept = 2.5, color = "grey", size = 0.5) +
  geom_hline(yintercept = 3.5, color = "grey", size = 0.5) +
  geom_hline(yintercept = 4.5, color = "grey", size = 0.5) +
  geom_hline(yintercept = 5.5, color = "grey", size = 0.5) +
  geom_hline(yintercept = 6.5, color = "grey", size = 0.5) +
  geom_hline(yintercept = 7.5, color = "grey", size = 0.5) +
  geom_vline(xintercept = 3.5, color = "grey", size = 0.5, linetype = 'dashed') +
  geom_vline(xintercept = 11.5, color = "grey", size = 0.5, linetype ='dashed') +
  geom_vline(xintercept = 15.5, color = "grey", size = 0.5, linetype ='dashed')
  # annotate("text", x = 3, y = 8.4, label = "0m") +
  # annotate("text", x = 10.5, y = 8.4, label = "50m") +
  # annotate("text", x = 14.5, y = 8.4, label = '100m') +
  # annotate("text", x = 18.5, y = 8.4, label = '150m')

cfds
# we need to add an invisible value to a sample in richmond park v3 and v4 so these appear on the figure
# TO DO: change level order of POint_ID and workout out why N is not being included.

rpds = bat_detections %>% filter(bat_detections$Location=='Richmond Park' & bat_detections$rp_include==TRUE) %>%
  ggplot(aes(x = Point_ID, y = SpeciesID, fill = read_bin)) +
  geom_tile() +
  scale_fill_manual(values = colour_pal, breaks = names(colour_pal[-1])) +
  facet_wrap(~month, nrow = 1) +
  theme_minimal() +
  labs(title = "b)", x = "Site ID", y = "Species ID", fill = "Total Reads") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        text = element_text(size = 16),
        axis.text.y = element_text(face = 'italic'),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        # axis.text.x = element_blank(),
        axis.title.x = element_blank()
  )+
  geom_hline(yintercept = 1.5, color = "grey", size = 0.5) +
  geom_hline(yintercept = 2.5, color = "grey", size = 0.5) +
  geom_hline(yintercept = 3.5, color = "grey", size = 0.5) +
  geom_hline(yintercept = 4.5, color = "grey", size = 0.5) +
  geom_hline(yintercept = 5.5, color = "grey", size = 0.5) +
  geom_vline(xintercept = 3.5, color = "grey", size = 0.5, linetype = 'dashed') +
  geom_vline(xintercept = 10.5, color = "grey", size = 0.5, linetype ='dashed') +
  geom_vline(xintercept = 13.5, color = "grey", size = 0.5, linetype ='dashed')
  # annotate("text", x = 3, y = 6.4, label = "0m") +
  # annotate("text", x = 10, y = 6.4, label = "50m") +
  # annotate("text", x = 12.5, y = 6.4, label = '100m') +
  # annotate("text", x = 15.5, y = 6.4, label = '150m')

rpds
cfds + rpds + plot_layout(nrow = 2, guides = 'collect')

library(ggh4x)


ggplot(bat_detections[bat_detections$read_bin!='0',], aes(x = Point_ID, y = SpeciesID, fill = read_bin)) +
  geom_tile() +
  scale_fill_manual(values = colour_pal, breaks = names(colour_pal[-1])) +
  facet_nested(Location~month, scales = 'free_y', drop= TRUE, nest_line = TRUE) +
  theme_minimal() +
  labs(title = "Bat Detections at Richmond Park", x = "Distance from bat roost (m)", y = "Species ID", fill = "Total Reads") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major.y = element_line(color = 'grey', size = 0.5),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 16),
        axis.text.y = element_text(face = 'italic'),
        axis.text.x = element_blank())
  # geom_hline(yintercept = 1.5, color = "grey", size = 0.5) +
  # geom_hline(yintercept = 2.5, color = "grey", size = 0.5) +
  # geom_hline(yintercept = 3.5, color = "grey", size = 0.5) +
  # geom_hline(yintercept = 4.5, color = "grey", size = 0.5)
  # geom_hline(yintercept = 5.5, color = "grey", size = 0.5) +
  # geom_vline(xintercept = 3.5, color = "grey", size = 0.5, linetype = 'dashed') +
  # geom_vline(xintercept = 10.5, color = "grey", size = 0.5, linetype ='dashed') +
  # geom_vline(xintercept = 13.5, color = "grey", size = 0.5, linetype ='dashed') +
  # annotate("text", x = 3, y = 6.4, label = "0m") +
  # annotate("text", x = 10, y = 6.4, label = "50m") +
  # annotate("text", x = 12.5, y = 6.4, label = '100m') +
  # annotate("text", x = 15.5, y = 6.4, label = '150m')
  #


# Combine detections with quant data ----------------------------------------
# note on the files in the folder -
# spatialQuant is without RP V4
# SpatialqQuant (1) is with RPV4
sq = read.csv('Data/NatureAir_SpatialQuant (1).csv')
head(sq)
# at the moment ignoring mock, EBA, and EBT
sq1 = sq %>%
  filter(!grepl("Mock", MiseqRunID)) %>%
  filter(!grepl("EBA|EBT", MiseqRunID))

# fixing a typo - again not sure if this Site I is missing a rep
sd_long1$Sample[sd_long1$Sample == "RPIN1V3rep1.1"] = "RPIN1V3rep3"
sq1[991,]$MiseqRunID = "RPIN1V3rep3"

# in the new dataset, this doesn't exist
# sd_long1$Sample[sd_long1$Sample == "CFLN1V2rep2.1"] = "CFLN1V2rep3"

# we need the pre-index concentration
# start by seeing what matches up
sq_id = unique(sq1$MiseqRunID)
sd_id = unique(sd_long1$Sample)
# which sq ids are not in sd_long?
missing_sq = sd_id[!sd_id %in% sq_id] #all the ids in sd_long, match with something in sq
missing_sd = sq_id[!sq_id %in% sd_id] #there are 64 samples that were quanted but not metabarcoding results
missing_ids #there are 6 samples in siteinfo that are not in sd_long

# there are lots of missing samples - assuming this is because nothing is detected.
# will need to include these as blanks at some point.

# join sq1 with sd_long1 by Sample and MiseqRunID
sd_long2 = left_join(sd_long1, sq1[,-2], by = c("Sample" = "MiseqRunID"), relationship = "many-to-many")

# change !Value to 0 - this is where the concentration was too low to estimate
sd_long2$PreIndexConcentration.nM.[sd_long2$PreIndexConcentration.nM.=="#VALUE!"] = 0
sd_long2$PreIndexConcentration.ng.ul.[sd_long2$PreIndexConcentration.ng.ul.=="NaN"] = 0
sd_long2$PostIndexConcentration.nM.[sd_long2$PostIndexConcentration.nM.=="#VALUE!"] = 0
sd_long2$PostIndexConcentration.ng.ul.[sd_long2$PostIndexConcentration.ng.ul.=="NaN"] = 0

sd_long2$PreIndexConcentration.nM. = as.numeric(sd_long2$PreIndexConcentration.nM.)
sd_long2$PreIndexConcentration.ng.ul. = as.numeric(sd_long2$PreIndexConcentration.ng.ul.)
sd_long2$PostIndexConcentration.nM. = as.numeric(sd_long2$PostIndexConcentration.nM.)
sd_long2$PostIndexConcentration.ng.ul. = as.numeric(sd_long2$PostIndexConcentration.ng.ul.)

# save this
write.csv(sd_long2, "output/spatial_data_longform_withquant.csv", row.names = FALSE)


# Combine species ID ------------------------------------------------------

# the other thing to do is combine species.
# there are 255 species level assignments, but many are likely to be the same species.
species_list = unique(sd_long2[,c('ASV_ID',"Class", "Order",'Family','Genus','Species', "Status")])

# lets group some of these together.
species_list$SpeciesID = ifelse(species_list$Species != "", species_list$Species, species_list$Genus)
species_list$SpeciesID = ifelse(species_list$SpeciesID != "", species_list$SpeciesID, species_list$Family)


species_list$SpeciesID[species_list$Status == 'Contaminant'] = "Contaminant"
species_list$SpeciesID[species_list$Status== "Unassigned"] = "Unassigned"

# but, we want to keep cows separate because there are actually cows in the area
species_list$SpeciesID[species_list$Genus== "Bos"] = "Bos taurus" #cow
# later on we will remove BOS labelled as NUMT

# there are some genus IDs that can only be one species in our location
species_list$SpeciesID[species_list$Genus== "Meles"] = "Meles meles" #badger
species_list$SpeciesID[species_list$Order== "Primates" & species_list$Species==""] = "Contaminant" #humans
species_list$SpeciesID[species_list$Genus== "Microtus"] = "Microtus agrestis" #vole
species_list$SpeciesID[species_list$Genus== "Myodes"] = "Myodes glareolus" #vole
species_list$SpeciesID[species_list$Family== "Sciuridae" & species_list$Species==""] = "Sciurus"
species_list$SpeciesID[species_list$Class== "Actinopterygii"] = "Fish"
species_list$SpeciesID[species_list$Genus== "Pica"] = "Pica pica" #magpie
species_list$SpeciesID[species_list$Class== "Aves" & species_list$Family==""] = "Unassigned bird"
# bufo - toads there is only one species
species_list$SpeciesID[species_list$Genus== "Bufo"] = "Bufo bufo"

# all birds that could not be assigned below family can be unassigned bird
species_list$SpeciesID[species_list$Class== "Aves" & species_list$Genus==""] = "Unassigned bird"

# now lets group some species together to genus level
# Apodemus - mice
species_list$SpeciesID[species_list$Genus== "Apodemus"] = "Apodemus spp."
# Sciurus - squirrels
species_list$SpeciesID[species_list$Family== "Sciuridae"] = "Sciurus spp."

# there is only one sample where pipistrelles nathusii was found, so lets group all pipistelles together
species_list$SpeciesID[species_list$Genus== "Pipistrellus"] = "Pipistrellus spp."
# myotis can be grouped together.
# all the myotis species have quite differnt ecology so its a shame to group
# them together. but there just isn't enough data to look at them separately
species_list$SpeciesID[species_list$Genus== "Myotis"] = "Myotis spp."

# group all birds together
species_list$SpeciesID[species_list$Class== "Aves"] = "Bird"

# we could group by functional group - deer, small mammal,
species_list$SpeciesID[species_list$Family=='Cervidae'] = "Cervidae"
# small mammals - hedgehog, rabbit, rodents, Soricomorpha (shrews),
species_list$SpeciesID[species_list$Order %in% c('Erinaceomorpha', "Lagomorpha", "Rodentia", 'Soricomorpha')] = "Small Mammal"
# foxes we keep seperate, Mustelidae we combine - the majority of this is badger,
# there are only three samples that are assigned to family only.
species_list$SpeciesID[species_list$Family== 'Mustelidae'] = "Mustelidae"
# amphibians
species_list$SpeciesID[species_list$Class== "Amphibia"] = "Amphibian"

# assign NUMT as NUMT
species_list$SpeciesID[species_list$Status== "NUMT"] = "Contaminant"
species_list$SpeciesID[grepl("Synthetic", species_list$Status)] = "Synthetic"
species_list$SpeciesID[species_list$Status == "Sytnthetic3"] = "Synthetic"

unique(species_list$SpeciesID)

# tidying up
species_list$SpeciesID[species_list$SpeciesID == "Rhinolophus"] = "Rhinolophus spp."
species_list$SpeciesID[species_list$SpeciesID == "Plecotus"] = "Plecotus spp."

# there is an ASV read for Plecotus auritus but no reads in any samples.

# i've got this down to 18 species groups.

# for now, join species_ID back up with sd_long2
sd_long3 = left_join(sd_long2, species_list[,c('ASV_ID','SpeciesID')], by = c("ASV_ID" = "ASV_ID"), relationship = "many-to-many")

# remove plecotus auritus as there are no reads for this species
sd_long3 = sd_long3 %>%
  filter(SpeciesID != "Plecotus auritus")
# check total read numbers for each species to see if there are any other zero values
sd_long3 %>%
  group_by(SpeciesID) %>%
  summarise(total_reads = sum(read_count))

# we also need to add in the missing samples. these are samples that are in siteinfo but not in sd_long3.
# for every species, need to give a 0 value.
# convert to include PCR reps
# from missing ids, the CF samples arev not in sd or sq.
# RPKN2V4, RPANN3V4, RPASN3V4, are in SQ and have NA values for concentration
# RPPN3V4 has conc values in Sq but is not in sd.
missing_ids2 = c(paste0(missing_ids, "rep1"), paste0(missing_ids, "rep2"),
                 paste0(missing_ids, "rep3"))
missing_sd
# remove extraction blanks
missing_sd = missing_sd[!grepl("CFEB", missing_sd)]

missing_sd = c(missing_ids2, missing_sd)
# check for duplicates in missing_sd
duplicates = missing_sd[duplicated(missing_sd)]
# remove duplicates from missing_s
missing_sd = missing_sd[!duplicated(missing_sd)]
missing_sample = unique(sd_long3[,1:14])

names(sd_long3)
missing_sq = sq1[1,]
missing_sq[,4:7] = NA

# "CFAN1V5rep3"
# this is getting into missing_sd.
# missing_sd
# these have been coverted to CFANN1V5... rep3 is missing from the metabarcoding data.
# therefore, a row for CFANN1V5rep3 needs to be added.
# the best way to do this is changing the value in missing_sd.
missing_sd[grep("CFAN1V5rep3", missing_sd)] = "CFANN1V5rep3"

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
    filter(MiseqRunID == missing_sd[i])
  if (nrow(sq_temp) > 0) { #if there is data for this sample in sq..
    temp3 = left_join(temp2, sq_temp[,-2], by = c("Sample" = "MiseqRunID"))
    temp3$PreIndexConcentration.nM.[temp3$PreIndexConcentration.nM.=="#VALUE!"] = 0
    temp3$PreIndexConcentration.ng.ul.[temp3$PreIndexConcentration.ng.ul.=="NaN"] = 0
    temp3$PostIndexConcentration.nM.[temp3$PostIndexConcentration.nM.=="#VALUE!"] = 0
    temp3$PostIndexConcentration.ng.ul.[temp3$PostIndexConcentration.ng.ul.=="NaN"] = 0

    temp3$PreIndexConcentration.nM. = as.numeric(temp3$PreIndexConcentration.nM.)
    temp3$PreIndexConcentration.ng.ul. = as.numeric(temp3$PreIndexConcentration.ng.ul.)
    temp3$PostIndexConcentration.nM. = as.numeric(temp3$PostIndexConcentration.nM.)
    temp3$PostIndexConcentration.ng.ul. = as.numeric(temp3$PostIndexConcentration.ng.ul.)

  } else {
    sq_temp = missing_sq
    sq_temp$MiseqRunID = missing_sd[i]
    sq_temp$Sample.ID = sampleid
    temp3 = left_join(temp2, sq_temp[,-2], by = c("Sample" = "MiseqRunID"))
  }
  temp3 = left_join(temp3, species_list[,c('ASV_ID','SpeciesID')], by = c("ASV_ID" = "ASV_ID"))
  sd_long4 = rbind(sd_long4, temp3)
}

siteinfo_ids = unique(siteinfo$SampleID)
sd_ids = unique(sd_long4$SampleID)
sq_ids = unique(sq$MiseqRunID)
sdsq_ids = unique(sd_long4$Sample)
# which samples in sd_long are missing from sq_ids?
missing = sdsq_ids[!sdsq_ids %in% sq_ids]

names(sd_long4)
unique(sd_long4$DominantLC_10m)
# open or closed habirat types
sd_long4 = sd_long4 %>%
  mutate(Habitat_type = case_when(DominantLC_10m %in% c('Broadleaf woodland') ~ 'Closed',
                                  DominantLC_10m %in% c('Improve grassland', 'Neutral grassland', 'Fen') ~ 'Open',
                                  TRUE ~ 'Other'))
table(sd_long4$Habitat_type)
table(sd_long4$SampleID)

write.csv(sd_long4, "output/spatial_data_longform_withquant_speciesID.csv", row.names = FALSE)
write.csv(siteinfo, "output/all_site_vars2.csv", row.names = FALSE)


#
sd_long4$primer = 'Mam16S'
# Bat specific primer -----------------------------------------------------

batsd = read.delim2('Data/NatureAirBats.txt')
# there are 510 ASVs!
bats_long = pivot_longer(batsd, cols = 15:ncol(batsd), names_to = 'Sample',
                         values_to = 'read_count')


# seperate by _, choose second option
bats_long$SampleID = sapply(strsplit(bats_long$Sample, "_"), function(x) x[2])
# sometimes the reps have been named weirdly i.e rep1_2. we need to parse these
# out
bats_long$pcr_replicate = gsub("_RD(\\d+).*", "", bats_long$Sample)
bats_long$pcr_replicate = gsub("ZG0006223_", "", bats_long$pcr_replicate)
# bats_long$pcr_replicate = gsub("_", ".", bats_long$pcr_replicate)
bats_long$pcr_replicate = gsub('r', "R", bats_long$pcr_replicate)
bats_long$pcr_replicate = gsub(".*Rep", "", bats_long$pcr_replicate)


# filter all samples that have a . in the pcr_replicate column -
# these are the ones that have been named weirdly and need to be fixed
probreps = bats_long[grepl("_", bats_long$pcr_replicate),]
unique(probreps$SampleID)
# kate confirmed that these samples had 4 PCR replicates.
unique(probreps$pcr_replicate)

teste = bats_long[grepl("CFNN1V4", bats_long$SampleID),]
unique(teste$pcr_replicate)
# for this case, '3_1 should be '3' and '3_2' should be 4
teste = bats_long[grepl("CFON1V4", bats_long$SampleID),]
unique(teste$pcr_replicate)
table(teste$pcr_replicate)
# there are equal row numbers for each replicate, so even though there are 6
# this would suggest that there are actually 6 replicates.
# so '1_1' is 1, '1_2' is 2, '2_1' is 3, '2_2' is 4, '3_1' is 5, and '3_2' is 6.
teste = bats_long[grepl("CFPN1V4", bats_long$SampleID),]
unique(teste$pcr_replicate)
# here there are five.
#  1 is 1, 2_1 is 2, 2_2 is 3, 3_1 is 4, and 3_2 is 5.

bats_long = bats_long %>%
  mutate(pcr_replicate = case_when(
    grepl("CFNN1V4", SampleID) & pcr_replicate == "3_1" ~ "3",
    grepl("CFNN1V4", SampleID) & pcr_replicate == "3_2" ~ "4",
    grepl("CFON1V4", SampleID) & pcr_replicate == "1_1" ~ "1",
    grepl("CFON1V4", SampleID) & pcr_replicate == "1_2" ~ "2",
    grepl("CFON1V4", SampleID) & pcr_replicate == "2_1" ~ "3",
    grepl("CFON1V4", SampleID) & pcr_replicate == "2_2" ~ "4",
    grepl("CFON1V4", SampleID) & pcr_replicate == "3_1" ~ "5",
    grepl("CFON1V4", SampleID) & pcr_replicate == "3_2" ~ "6",
    grepl("CFPN1V4", SampleID) & pcr_replicate == "2_1" ~ "2",
    grepl("CFPN1V4", SampleID) & pcr_replicate == "2_2" ~ "3",
    grepl("CFPN1V4", SampleID) & pcr_replicate == "3_1" ~ "4",
    grepl("CFPN1V4", SampleID) & pcr_replicate == "3_2" ~ "5",
    TRUE ~ pcr_replicate
  ))

unique(bats_long$pcr_replicate)

# remove string from Sample column containing 'rep' followed by a number.
# Add the number to a new column called 'replicate'
bats_long$SampleID = gsub('r', "R", bats_long$SampleID)
bats_long = bats_long %>%
  mutate(SampleID = gsub("Rep\\d+", "", SampleID))


# CFAN1V5 CFAN2V5 CFAN3V5
# there is no V5 for AS or AN, Kate confirmed only one 'A' sample was taken for this visit
# relabelling as AN so it is easier to link to site info, as i did above.
bats_long$SampleID[grep("CFAN1V5", bats_long$SampleID)] = "CFANN1V5"
bats_long$SampleID[grep("CFAN2V5", bats_long$SampleID)] = "CFANN2V5"
bats_long$SampleID[grep("CFAN3V5", bats_long$SampleID)] = "CFANN3V5"

# "CFBNN2V1" - this could be CFBN2V1 rep3 or CFNNN2V1 rep3.. both seem to be missing.
# relabelling as CFB for now.
bats_long$SampleID[grep("CFBNN2V1", bats_long$SampleID)] = "CFBN2V1"
bats_long$SampleID[grep("CFFEN2V1", bats_long$SampleID)] = "CFEN2V1"
bats_long$SampleID[grep('CFO3V1', bats_long$SampleID)] = "CFON3V1"

# left join with site_info
# join with metadata from siteinfo
bats_long1 = left_join(bats_long, siteinfo, by = c("SampleID" = "SampleID"), relationship = "many-to-many")

# check which samples did not line up.
mismatch = bats_long1 %>%
  filter(is.na(Point_ID)) %>%
  filter(!grepl("Mock", SampleID)) %>% #not interested in mock community
  filter(!grepl("EBA|EBT", SampleID)) # not interested in extraction blanks

unique(mismatch$SampleID)

# do all the rows in siteinfo match up with something?
siteinfo_ids = unique(siteinfo$SampleID[siteinfo$Location=="Canada Farm"])
sd_ids = unique(bats_long$SampleID)

missing_ids = siteinfo_ids[!siteinfo_ids %in% sd_ids] #these are the samples
# that are in siteinfo but not in the bats_long dataset.
unmatches = sd_ids[!sd_ids %in% siteinfo_ids]
# these are samples that we collected but are not in siteinfo - all mocks

# sort out species ID. These will need to match up with categories from
# the main dataset.
unique(bats_long1$Species)
bats_sp_list = unique(bats_long1[,c('NMSeqID',"Class", "Order",'Family','Genus','Species', "Status")])

# lets group some of these together.
bats_sp_list$SpeciesID = ifelse(bats_sp_list$Species != "", bats_sp_list$Species, bats_sp_list$Genus)
bats_sp_list$SpeciesID = ifelse(bats_sp_list$SpeciesID != "", bats_sp_list$SpeciesID, bats_sp_list$Family)

bats_sp_list$SpeciesID[bats_sp_list$Status == 'Contaminant'] = "Contaminant"
bats_sp_list$SpeciesID[bats_sp_list$Status== "Unassigned"] = "Unassigned"

# but, we want to keep cows separate because there are actually cows in the area
bats_sp_list$SpeciesID[bats_sp_list$Genus== "Bos"] = "Bos taurus" #cow

unique(bats_sp_list$SpeciesID)
unique(species_list$SpeciesID)
# in the previous dataset, all bird species are grouped together

bats_sp_list$SpeciesID[bats_sp_list$Status== "NUMT"] = "Contaminant"
bats_sp_list$SpeciesID[bats_sp_list$Class == 'Aves'] = "Bird"
bats_sp_list$SpeciesID[bats_sp_list$Genus == 'Pipistrellus'] = "Pipistrellus spp."
bats_sp_list$SpeciesID[bats_sp_list$Genus == 'Myotis'] = "Myotis spp."
bats_sp_list$SpeciesID[bats_sp_list$SpeciesID == "Scardinius erythrophthalmus"] = "Fish"
bats_sp_list$SpeciesID[bats_sp_list$SpeciesID == "Cyprinidae"] = "Fish"
bats_sp_list$SpeciesID[bats_sp_list$SpeciesID == "Salmo salar"] = "Fish"
bats_sp_list$SpeciesID[bats_sp_list$Class == "Amphibia"] = "Amphibian"
# Mustelidae
bats_sp_list$SpeciesID[bats_sp_list$Family == "Mustelidae"] = "Mustelidae"
bats_sp_list$SpeciesID[bats_sp_list$Family=='Cervidae'] = "Cervidae"
# small mammals - hedgehog, rabbit, rodents, Soricomorpha (shrews),
bats_sp_list$SpeciesID[bats_sp_list$Order %in% c('Erinaceomorpha', "Lagomorpha", "Rodentia", 'Soricomorpha')] = "Small Mammal"
# "Aporrectodea"  #Invertebrate - earth worm
bats_sp_list$SpeciesID[bats_sp_list$SpeciesID == "Aporrectodea"] = "Earthworm"
# "Sabellaria"  marine polychaete worm  #contaminant
bats_sp_list$SpeciesID[bats_sp_list$SpeciesID == "Sabellaria"] = "Contaminant"
# "Helophilus pendulus" hover fly

bats_sp_list[bats_sp_list$SpeciesID == "",] #synthetic
bats_sp_list$SpeciesID[bats_sp_list$SpeciesID == ""] = "Synthetic"

bats_sp_list[bats_sp_list$SpeciesID == "Vulpes",] #this can only be red fox in UK
bats_sp_list$SpeciesID[bats_sp_list$SpeciesID == "Vulpes"] = "Vulpes vulpes"

# combine speices ID with sample data
# also combine quant data
# join sq1 with sd_long1 by Sample and MiseqRunID

# we cannot connect quant data directly with reps as they come from different samples.
# so instead we will just put averages from each sample and connect them

# do all sq1 samples match up with bats_long2?
sq1samples = unique(sq1$Sample.ID[grepl("CF", sq1$Sample.ID)])
bats_long1samples = unique(bats_long1$SampleID)

# what samples are in sq1 but not bats_long?
# this means that there are samples that were quanted but not metabarcoded.
missing_sq1 = sq1samples[!sq1samples %in% bats_long1samples]
# what samples are in bats_long but not sq?
# this means there is probably a miss-labelled sample.
mssing_bats = bats_long1samples[!bats_long1samples %in% sq1samples] #just mocks
# just the ANs V5
sq1$Sample.ID[grepl("CFAN1", sq1$Sample.ID)] = 'CFANN1V5'
sq1$Sample.ID[grepl("CFAN2", sq1$Sample.ID)] = 'CFANN2V5'
sq1$Sample.ID[grepl("CFAN3", sq1$Sample.ID)] = 'CFANN3V5'
sq1$Sample.ID = gsub('CFAA', 'CFAN', sq1$Sample.ID)

# combine quant data with bats_long1
# get average quant data for each sample
sq1_avg = sq1 %>%
  group_by(Sample.ID) %>%
  summarise(PreIndexConcentration.ng.ul. = mean(PreIndexConcentration.ng.ul., na.rm = TRUE))

bats_long2 = left_join(bats_long1, sq1_avg, by = c("SampleID" = "Sample.ID"), relationship = "many-to-many")

# change !Value to 0 - this is where the concentration was too low to estimate
bats_long2$PreIndexConcentration.ng.ul.[bats_long2$PreIndexConcentration.ng.ul.=="NaN"] = 0
bats_long2$PreIndexConcentration.ng.ul. = as.numeric(bats_long2$PreIndexConcentration.ng.ul.)

names(sd_long4)
names(bats_long2)
sd_long4$MiseqNo.[2]
bats_long2$MiseqNo = NA

# add SpeciesID
bats_long3 = left_join(bats_long2, bats_sp_list[,c('NMSeqID','SpeciesID')], by = c("NMSeqID" = "NMSeqID"), relationship = "many-to-many")
names(sd_long4)
names(bats_long3)

bats_long3$primer = 'Bat'

# need to add information from missing samples - those that were analyses but got 0 DNA

# we also need to add in the missing samples. these are samples that are in siteinfo but not in bats_long3.
# for every species, need to give a 0 value.
# convert to include PCR reps
# from missing ids, the CF samples are not in sd or sq.

# those that are in site info but not in the metabarcoding data are:
missing_ids2 = c(paste0(missing_ids, "rep1"), paste0(missing_ids, "rep2"),
                 paste0(missing_ids, "rep3"))

# those that are in sq but not in metabarcoding are:
missing_sq1
# remove extraction blanks
missing_sq1 = missing_sq1[!grepl("CFEB", missing_sq1)]
# add rep numbers
missing_sq2 = c(paste0(missing_sq1, "rep1"), paste0(missing_sq1, "rep2"),
                 paste0(missing_sq1, "rep3"))


missing_sd = c(missing_ids2, missing_sq2)

missing_sample = unique(bats_long3[,1:14])

names(bats_long3)
missing_sq = sq_avg[1,] #we only care about preIndexConcentration.ng.ul.
missing_sq[,4] = NA

bats_long4 = bats_long3
bats_long4$pcr_replicate = as.numeric(bats_long4$pcr_replicate)
# remove miseqno columns
bats_long4 = bats_long4 %>%
  select(-MiseqNo) %>%
  rename(biomassInSample = PreIndexConcentration.ng.ul.)

# check for duplicates in missing_sd
missing_sd[duplicated(missing_sd)]
missing_sd = missing_sd[!duplicated(missing_sd)]


for(i in 1:length(missing_sd)){
  # remove rep
  sampleid = gsub("rep\\d+", "", missing_sd[i])
  missing_sample_temp = missing_sample %>%
    mutate(Sample = missing_sd[i],
           read_count = 0,
           SampleID = sampleid,
           pcr_replicate = as.numeric(gsub(".*rep(\\d+).*", "\\1", missing_sd[i]))
           )
  temp = siteinfo %>%
    filter(SampleID == sampleid)

  temp2 = left_join(missing_sample_temp, temp, by = 'SampleID')

  # add sq data - again we just do the average for the sample, not the rep.
  sq_temp = sq_avg %>%
    filter(SampleID == sampleid)
  if (nrow(sq_temp) > 0) { #if there is data for this sample in sq..
    temp2$biomassInSample = sq_temp$biomassInSample
#     if its NA make it 0 instead
    temp2$biomassInSample[is.na(temp2$biomassInSample)] = 0
    temp2$biomassInSample = as.numeric(temp2$biomassInSample)

  } else {
    sq_temp = missing_sq
    temp2$biomassInSample = sq_temp$biomassInSample
  }
  temp3 = left_join(temp2, bats_sp_list[,c('NMSeqID','SpeciesID')], by = c("NMSeqID" = "NMSeqID"))
  temp3$primer = 'Bat'
  bats_long4 = rbind(bats_long4, temp3)
}

siteinfo_ids = unique(siteinfo$SampleID)
sd_ids = unique(bats_long4$SampleID)
sq_ids = unique(sq1$Sample.ID)
sdsq_ids = unique(bats_long4$SampleID)
# which samples in sd_long are missing from sq_ids?
missing = sdsq_ids[!sdsq_ids %in% sq_ids]

# remove mocks
bats_long4 = bats_long4 %>%
  filter(!grepl("Mock", SampleID))

table(bats_long4$Sample)
write.csv(bats_long4, "output/bat_primer_longform_withquant_speciesID.csv", row.names = FALSE)
