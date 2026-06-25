
library(sf)
library(rstan)
library(bayesplot)
library(tidyr)
library(dplyr)
library(patchwork)
# model plotting
#
stan = readRDS("output/model_output/results_stan_batspecies_Canada_Farm_mcmc_20260326.rds")


# check model convergence
rstan::traceplot(stan, pars = c("tau", "sigma", "mu0"))
rstan::check_hmc_diagnostics(stan)
print(stan, pars = c("beta_z", "beta_theta", "tau", "sigma"))

posterior = as.array(stan)
dim(posterior)
dimnames(posterior)
mcmc_intervals(posterior, pars = c("beta_z[1,1]"))
# the chains are not mixing well, except for mu0
range(rhat(stan))
# data
sd = read.csv("output/spatial_data_longform_withquant_speciesID.csv")
siteinfo = read.csv("output/all_site_vars2.csv")
siteinfo = siteinfo %>%
  mutate(Habitat_type = case_when(DominantLC_10m %in% c('Broadleaf woodland') ~ 'Closed',
                                  DominantLC_10m %in% c('Improve grassland', 'Neutral grassland', 'Fen') ~ 'Open',
                                  TRUE ~ 'Other'))

# remove thresholded read counts
# sum lib size of each sample
libsizes = sd %>%
  group_by(Sample) %>%
  summarise(total_libsize = sum(read_count))
libsizes$thresh = libsizes$total_libsize * 0.0005
sd = sd %>%
  left_join(libsizes, by = c("Sample")) %>%
  filter(read_count >= thresh)

# new sampleID
sd$SampleIDrep = paste0(sd$SampleID, "_", sd$pcr_replicate)


# for now, just do canadafarm
sd = sd %>% filter(Location == "Canada Farm")
siteinfo = siteinfo %>% filter(Location == "Canada Farm") %>% droplevels()

# at the moment, the model is only looking at bats
# sd = sd %>% filter(Order == "Chiroptera") %>% droplevels()

# get species order
sd_summ = sd %>%
  group_by(Point_ID, Visit, Night, pcr_replicate, SpeciesID) %>%
  summarise(reads = sum(read_count)) %>%
  arrange(Point_ID, Visit, Night, pcr_replicate, SpeciesID) %>%
  ungroup()
sd_summ2 = pivot_wider(sd_summ, names_from = SpeciesID, values_from = reads, values_fill = 0)
logy1 = as.matrix(sd_summ2[, -(1:4)])
# remove columns with no detections (colSUm = 0)
logy1 = logy1[, colSums(logy1) > 0]
species_names = colnames(logy1)

sites_ecol = sites_ecol %>%
  arrange(Point_ID, Visit)
sites_ecol = siteinfo %>%
  select(Point_ID, Visit, dist_m, Habitat_type) %>%
  distinct()

# prep model results
matrix_results <- as.matrix(stan) %>%
  as.data.frame

texts <- colnames(matrix_results)
t = 5

# i also want total number of detections for each species
sd_bats = sd %>% filter(Order == "Chiroptera") %>% droplevels()  %>% filter(read_count > 0)

detection_summ = sd_bats %>%
  group_by(SpeciesID) %>%
  filter(read_count > 0) %>%
  summarise(detections = n_distinct(SampleIDrep)) %>%
  arrange(SpeciesID) %>%
  filter(detections > 0)
detection_visits = sd_bats %>%
  group_by(SpeciesID, Visit) %>%
  summarise(detections = n_distinct(SampleIDrep[read_count > 0])) %>%
  arrange(SpeciesID, Visit)


# effect of distancef

{
  # select only columns using distance coefficient. this is the same as t, as there are t-1 columns for each time point and
  # then one column for distance.
  distance_effects <- matrix_results[,grepl(paste0("beta_z\\[",t), texts) ]

  # on each column, calculate the 95% credible intervals & mean
  # t here means transpose to give one row per species
  beta_z_results <- distance_effects %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame
  beta_z_results$mean = colMeans(distance_effects)
  beta_z_results$species = species_names
  beta_z_results = arrange(beta_z_results, mean)
  # add detection numbers

  # connect species names to columns
  # only plot bats
  bats = c("Eptesicus serotinus", "Myotis spp.", "Nyctalus","Pipistrellus spp.", "Plecotus", "Rhinolophus")
  betaz_bats = beta_z_results %>%
    filter(species %in% bats) %>%
    arrange(mean)
  # add detection results
  betaz_bats = betaz_bats %>%
    left_join(detection_summ, by = c("species" = "SpeciesID"))
  betaz_bats$species = factor(betaz_bats$species,
                                  levels = c("Myotis spp.","Pipistrellus spp.",
                                             "Eptesicus serotinus","Rhinolophus"))

  distanceef = ggplot(betaz_bats, aes(x = species,
                                          y = mean,
                                          ymin = `2.5%`,
                                          ymax = `97.5%`)) +
    geom_errorbar(linewidth = 1) +
    geom_point(size = 4) +
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of distance from roost on \n DNA biomass') +
    theme_classic() +
    ylim(-6.5,4) +
    geom_text(aes(label = paste0("n=", detections)),
              y = 3.5,
              size = 5) +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(angle = 20, hjust = 0.6, vjust = 0.75)
          # axis.title =element_text(size = 18),
          # legend.text = element_text(size = 14),
          ) +
    ggtitle("Canada Farm")


  distanceef
  }
# now look at effect of time.
# effect of distance
t-1
visit_effects <- matrix_results[,grepl("beta_z\\[[1-4],", texts) ]
month_cols = c(
  "March" = "#440154",
  "May" = "#31688e",
  "June" = "#21908CFF",
  "July" = "#35b779",
  "August" = "#fde725",
  "September" = "#fbb03b",
  "October" = "#ff7f0e"
)


{
  # pattern <- sprintf("\\[%d,", 1)
  # pattern = paste(sprintf("\\[%d,", 2), collapse = "|")
  # selected <- grepl(pattern, texts) & grepl("beta_z", texts)
  # on each column, calculate the 95% credible intervals
  # beta_z_results <- matrix_results[,selected] %>%
  beta_z_results <- visit_effects %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame
  beta_z_results$mean = colMeans(visit_effects)

  # not including march as this is reference value
  Month = c("May","June", "July", "October")
  beta_z_results$Species <- rep(species_names, each = t - 1)
  beta_z_results$Month <- factor(rep(Month, times = n_distinct(species_names)), levels = Month)

  # just bats
  beta_z_results = beta_z_results %>%
    filter(Species %in% bats)
  beta_z_results$Species = factor(beta_z_results$Species,
                                  levels = c("Myotis spp.","Pipistrellus spp.",
                                             "Eptesicus serotinus","Rhinolophus"))


  seasonef = ggplot(beta_z_results, aes(x = Species,
                                        y = mean,
                                        ymin = `2.5%`,
                                        ymax = `97.5%`,
                                        colour = Month)) +
    geom_errorbar(position = position_dodge(width = 0.7)) +
    geom_point(size = 3, position = position_dodge(width = 0.7)) +
    # scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen")) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # scale_x_continuous(breaks = 1:t, labels = Visits) +
    xlab('Species') +
    ylab('Effect of season on DNA biomass') +
    # facet_wrap(~Species, scales = 'free') +
    ylim(-6.5,4) +
    theme_classic() +
    scale_color_manual(values = month_cols, breaks = levels(beta_z_results$Month)) +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(angle = 20, hjust = 0.6, vjust = 0.75),
          # axis.title =element_text(size = 12),
          # legend.text = element_text(size = 12),
          strip.background = element_blank())
          # strip.text = element_text(size = 12)) +
    # ylim(-1,8) +

  seasonef
}

{

  # select only columns using distance coefficient. this is the same as t, as there are t-1 columns for each time point and
  # then one column for distance.
  hab_effects <- matrix_results[,grepl(paste0("beta_z\\[",t+1), texts) ]

  # on each column, calculate the 95% credible intervals & mean
  # t here means transpose to give one row per species
  beta_z_results <- hab_effects %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame
  beta_z_results$mean = colMeans(hab_effects)
  beta_z_results$species = species_names
  beta_z_results = arrange(beta_z_results, mean)
  # add detection numbers

  # connect species names to columns
  # only plot bats
  betaz_bats = beta_z_results %>%
    filter(species %in% bats) %>%
    arrange(mean)
  # add detection results
  betaz_bats = betaz_bats %>%
    left_join(detection_summ, by = c("species" = "SpeciesID"))
  betaz_bats$species = factor(betaz_bats$species, levels = betaz_bats$species)

  habef = ggplot(betaz_bats, aes(x = species,
                                 y = mean,
                                 ymin = `2.5%`,
                                 ymax = `97.5%`)) +
    geom_errorbar(linewidth = 1) +
    geom_point(size = 4) +
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of open habitat on DNA biomass') +
    theme_classic() +
    geom_text(aes(label = paste0("n=", detections)),
              y = 3.5,
              size = 5) +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(angle = 20, hjust = 0.6, vjust = 0.75),
          axis.title =element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.position = "top") +
    ylim(-6.5,4)


  habef

}
distanceef + habef + seasonef + plot_layout()
# actual DNA biomass
logl_results = matrix_results %>%
  as.data.frame()

texts <- colnames(logl_results)
loglsel <- grepl('logl', texts)
# for each species, there are 94 estimates of logl
# this will be in the order of 1-19 for visit 1, visit 2 etc, or visit 1-4 for site 1, site 2 etc.
# it is the same format as the original data - sites_ecol, repeated for each species.

logl_results = logl_results[,loglsel] %>%
  apply(., 2, function(x) c(
    mean = mean(x),
    quantile(x, probs = c(0.025, 0.975)))) %>% t %>%
  as.data.frame()

logl_results$species = rep(species_names, each = nrow(sites_ecol))
#merge cfvarssub with loglresults (each needs to repeat)
logl_results$Visit = rep(sites_ecol$Visit, times = length(species_names))
logl_results$Point_ID = rep(sites_ecol$Point_ID, times = length(species_names))
library(sf)
# get unique site coordinates
sites_unique = siteinfo %>%
  group_by(Point_ID, Longitude, Latitude, dist_m) %>%
  filter(row_number() == 1)
# join longitude and latitude columns from sites ecol using point_id
logl_results2 = logl_results %>%
  left_join(sites_unique[,c("Point_ID", "Longitude", "Latitude", "dist_m")],
            by = 'Point_ID')
logl_sf = st_as_sf(logl_results2, coords = c("Longitude", "Latitude"), crs = 4326)

ggplot(data = logl_sf) +
  geom_sf(aes(size = mean), alpha = 0.6, colour = 'grey') +
  theme_minimal() + facet_grid(species ~ Visit)

# relabel visit
logl_sf$Visit <- factor(
  logl_sf$Visit,
  levels = 1:5,
  labels = c("March", "May", "June", "July", "October")
)

library(ggh4x)
logl_sf_bat = filter(logl_sf, species %in% bats)
logl_sf_bat$species = factor(logl_sf_bat$species, levels = c("Myotis spp.","Pipistrellus spp.",
           "Eptesicus serotinus","Rhinolophus"))
ggplot(data = logl_sf_bat, aes(x = dist_m, y = mean, ymax = `97.5%`, ymin = `2.5%`,
                           fill = Visit,
                           colour = Visit)) +
  geom_smooth() +
  facet_wrap(~species) +
  # ggh4x::facet_grid2(species~Visit, labeller = labeller(species = species_names), scales = 'free_x', independent = "x") +
  labs(fill = 'Visit', colour = 'Visit', x = 'Distance from roost (m)',
       y = 'DNA Biomass') +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_classic() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
  )

# logl is a proxy for how much biomass is at the site. but we don't know amplification rates, we don't
# know shedding rate. so its on an abstract scale because we cannot actually measure how much DNA is in the area.
ggplot(data = logl_sf[logl_sf$species %in% bats,], aes(x = dist_m, y = exp(mean), ymax = `97.5%`, ymin = `2.5%`)) +
  geom_smooth() +
  facet_wrap(~species, scales = 'free') +
  labs(x = 'Distance from roost (m)',
       y = 'DNA Biomass') +
  ylim(0,0.65) +
  theme_classic() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  ggtitle('Model: Approximate Sampling, all species, Canda Farm')


# plot effect of open habitat




# Bats as one group - Canada Farm -----------------------------------------
# Not up to date as of 30th Jan

# model plotting
stan_grp = readRDS("output/model_output/results_stan_batgrpd_Canada_Farm_mcmc_20260401.rds")
# get species order

# prep model results
matrix_grp <- as.matrix(stan_grp) %>%
  as.data.frame

texts <- colnames(matrix_grp)
t = 5

sd$SpeciesID[sd$Order == "Chiroptera"] = "Chiroptera"

# get species order
sd_summ = sd %>%
  group_by(Point_ID, Visit, Night, pcr_replicate, SpeciesID) %>%
  summarise(reads = sum(read_count)) %>%
  arrange(Point_ID, Visit, Night, pcr_replicate, SpeciesID) %>%
  ungroup()
sd_summ2 = pivot_wider(sd_summ, names_from = SpeciesID, values_from = reads, values_fill = 0)
logy1 = as.matrix(sd_summ2[, -(1:4)])
# remove columns with no detections (colSUm = 0)
logy1 = logy1[, colSums(logy1) > 0]
species_names = colnames(logy1)

sites_ecol = sites_ecol %>%
  arrange(Point_ID, Visit)
sites_ecol = siteinfo %>%
  select(Point_ID, Visit, dist_m, Habitat_type) %>%
  distinct()

# i also want total number of detections for each species
sd_bats = sd %>% filter(Order == "Chiroptera") %>% droplevels()
detection_sumb = sd_bats %>%
  summarise(detections = n_distinct(SampleIDrep[read_count > 0])) %>%
  filter(detections > 0)
detection_visitsb = sd_bats %>%
  group_by(Visit) %>%
  summarise(detections = n_distinct(SampleIDrep[read_count > 0])) %>%
  arrange(Visit)

# effect of distance
{

  # select only columns using distance coefficient. this is the same as t, as there are t-1 columns for each time point and
  # then one column for distance.
  distance_effects <- matrix_grp[,grepl(paste0("beta_z\\[",t), texts) ]

  # on each column, calculate the 95% credible intervals & mean
  # t here means transpose to give one row per species
  beta_z_results <- distance_effects %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame
  beta_z_results$mean = colMeans(distance_effects)
  beta_z_results$species = species_names
  beta_z_results = arrange(beta_z_results, mean)
  # add detection numbers

  # connect species names to columns
  # only plot bats
  betaz_bats = beta_z_results %>%
    filter(species == 'Chiroptera') %>%
    arrange(mean)
  # add detection results
  betaz_bats$detections = detection_sumb$detections
  betaz_bats$species = factor(betaz_bats$species, levels = betaz_bats$species)

  distanceef = ggplot(betaz_bats, aes(x = species,
                                      y = mean,
                                      ymin = `2.5%`,
                                      ymax = `97.5%`)) +
    geom_errorbar(linewidth = 1) +
    geom_point(size = 4) +
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of distance from roost on \n DNA biomass') +
    theme_classic() +
    geom_text(aes(label = paste0("n=", detections)),
              y = 3.5,
              size = 5) +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(angle = 20, hjust = 0.6, vjust = 0.75),
          axis.title =element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.position = "top") +
    ylim(-6.5,4)

  distanceef

  }

# now look at effect of time.
# effect of distance
t-1
visit_effects <- matrix_grp[,grepl("beta_z\\[[1-4],", texts) ]

{
  # pattern <- sprintf("\\[%d,", 1)
  # pattern = paste(sprintf("\\[%d,", 2), collapse = "|")
  # selected <- grepl(pattern, texts) & grepl("beta_z", texts)
  # on each column, calculate the 95% credible intervals
  # beta_z_results <- matrix_results[,selected] %>%
  beta_z_results <- visit_effects %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame
  beta_z_results$mean = colMeans(visit_effects)

  # not including march as this is reference value
  Months = c("May","June", "July", "October")
  beta_z_results$Species <- rep(species_names, each = t - 1)
  beta_z_results$Month <- factor(rep(Months, times = n_distinct(species_names)), levels = Months)

  # just bats
  beta_z_results = beta_z_results %>%
    filter(Species == 'Chiroptera')


  seasonef = ggplot(beta_z_results, aes(x = Species,
                                        y = mean,
                                        ymin = `2.5%`,
                                        ymax = `97.5%`,
                                        colour = Month)) +
    geom_errorbar(position = position_dodge(width = 0.7)) +
    geom_point(size = 3, position = position_dodge(width = 0.7)) +
    # scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen")) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # scale_x_continuous(breaks = 1:t, labels = Visits) +
    xlab('Species') +
    ylab('Effect of season on DNA biomass') +
    # facet_wrap(~Species, scales = 'free') +
    ylim(-6.5,4.15) +
    theme_classic() +
    scale_color_manual(values = month_cols, breaks = levels(beta_z_results$Month)) +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(angle = 20, hjust = 0.6, vjust = 0.75),
          # axis.title =element_text(size = 12),
          # legend.text = element_text(size = 12),
          strip.background = element_blank())
  # strip.text = element_text(size = 12)) +
  # ylim(-1,8) +

  seasonef
}

{

  # select only columns using distance coefficient. this is the same as t, as there are t-1 columns for each time point and
  # then one column for distance.
  hab_effects <- matrix_grp[,grepl(paste0("beta_z\\[",t+1), texts) ]

  # on each column, calculate the 95% credible intervals & mean
  # t here means transpose to give one row per species
  beta_z_results <- hab_effects %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame
  beta_z_results$mean = colMeans(hab_effects)
  beta_z_results$species = species_names
  beta_z_results = arrange(beta_z_results, mean)
  # add detection numbers

  # connect species names to columns
  # only plot bats
  betaz_bats = beta_z_results %>%
    filter(species == "Chiroptera") %>%
    arrange(mean)
  # add detection results

  betaz_bats$detections = detection_sumb$detections
  betaz_bats$species = factor(betaz_bats$species, levels = betaz_bats$species)


  habef = ggplot(betaz_bats, aes(x = species,
                                 y = mean,
                                 ymin = `2.5%`,
                                 ymax = `97.5%`)) +
    geom_errorbar(linewidth = 1) +
    geom_point(size = 4) +
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of open habitat on DNA biomass') +
    theme_classic() +
    geom_text(aes(label = paste0("n=", detections)),
              y = 3.5,
              size = 5) +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(angle = 20, hjust = 0.6, vjust = 0.75),
          axis.title =element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.position = "top") +
    ylim(-6.5,4)


  habef

}
distanceef + habef + seasonef + plot_layout()

# plot effect of distance, with real data underneath

logl_results = matrix_grp %>%
  as.data.frame()

texts <- colnames(logl_results)
loglsel <- grepl('logl', texts)
# for each species, there are 94 estimates of logl
# this will be in the order of 1-19 for visit 1, visit 2 etc, or visit 1-4 for site 1, site 2 etc.
# it is the same format as the original data - sites_ecol, repeated for each species.

logl_results = logl_results[,loglsel] %>%
  apply(., 2, function(x) c(
    mean = mean(x),
    quantile(x, probs = c(0.025, 0.975)))) %>% t %>%
  as.data.frame()

logl_results$species = rep(species_names, each = nrow(sites_ecol))
#merge cfvarssub with loglresults (each needs to repeat)
logl_results$Month = rep(sites_ecol$Visit, times = length(species_names))
logl_results$Point_ID = rep(sites_ecol$Point_ID, times = length(species_names))
library(sf)

# get unique site coordinates
sites_unique = siteinfo %>%
  group_by(Point_ID, Longitude, Latitude, dist_m) %>%
  filter(row_number() == 1)
# join longitude and latitude columns from sites ecol using point_id
logl_results2 = logl_results %>%
  left_join(sites_unique[,c("Point_ID", "Longitude", "Latitude", "dist_m")],
            by = 'Point_ID')

# relabel visit
logl_results2$Month <- factor(
  logl_results2$Month,
  levels = 1:5,
  labels = c("March", "May", "June", "August", "October")
)

logl_results2 = filter(logl_results2, species == 'Chiroptera')

# logl is a proxy for how much biomass is at the site. but we don't know amplification rates, we don't
# know shedding rate. so its on an abstract scale because we cannot actually measure how much DNA is in the area.
ggplot(data = logl_results2, aes(x = dist_m, y = mean, ymax = `97.5%`, ymin = `2.5%`,
                                                       fill = Month,
                                                       colour = Month)) +
  geom_smooth() +
  labs(fill = 'Visit', colour = 'Visit', x = 'Distance from roost (m)',
       y = 'DNA Biomass') +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_classic() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))


# Richmond park - all bat species -----------------------------------------
# Plotting richmond park results

# all bat species
# model plotting
stanrp = readRDS("output/model_output/results_stan_batspecies_Richmond_Park_mcmc_20260326.rds")
# data
sd = read.csv("output/spatial_data_longform_withquant_speciesID.csv")
siteinfo = read.csv("output/all_site_vars2.csv")

# new sampleID
sd$SampleIDrep = paste0(sd$SampleID, "_", sd$pcr_replicate)

sd = sd %>% filter(Location == "Richmond Park")
siteinfo = siteinfo %>% filter(Location == "Richmond Park") %>% droplevels()
siteinfo = siteinfo %>%
  mutate(Habitat_type = case_when(DominantLC_10m %in% c('Broadleaf woodland') ~ 'Closed',
                                  DominantLC_10m %in% c('Improve grassland', 'Neutral grassland', 'Fen') ~ 'Open',
                                  TRUE ~ 'Other'))

# remove thresholded read counts
# sum lib size of each sample
libsizes = sd %>%
  group_by(Sample) %>%
  summarise(total_libsize = sum(read_count))
libsizes$thresh = libsizes$total_libsize * 0.0005
sd = sd %>%
  left_join(libsizes, by = c("Sample")) %>%
  filter(read_count >= thresh)

# get species order
sd_summ = sd %>%
  group_by(Point_ID, Visit, Night, pcr_replicate, SpeciesID) %>%
  summarise(reads = sum(read_count)) %>%
  arrange(Point_ID, Visit, Night, pcr_replicate, SpeciesID) %>%
  ungroup()
sd_summ2 = pivot_wider(sd_summ, names_from = SpeciesID, values_from = reads, values_fill = 0)
logy1 = as.matrix(sd_summ2[, -(1:4)])
# remove columns with no detections (colSUm = 0)
logy1 = logy1[, colSums(logy1) > 0]
species_names = colnames(logy1)

sites_ecol = siteinfo %>%
  select(Point_ID, Visit, dist_m, Habitat_type) %>%
  distinct()
sites_ecol = sites_ecol %>%
  arrange(Point_ID, Visit)

# prep model results
matrix_results <- as.matrix(stanrp) %>%
  as.data.frame

texts <- colnames(matrix_results)
t = 4

# i also want total number of detections for each species
sd_bats = sd %>% filter(Order == "Chiroptera") %>% droplevels()
detection_summ = sd_bats %>%
  group_by(SpeciesID) %>%
  summarise(detections = n_distinct(SampleIDrep[read_count > 0])) %>%
  arrange(SpeciesID) %>%
  filter(detections > 0)
detection_visits = sd_bats %>%
  group_by(SpeciesID, Visit) %>%
  summarise(detections = n_distinct(SampleIDrep[read_count > 0])) %>%
  arrange(SpeciesID, Visit)

# check model performance
rstan::traceplot(stanrp, pars = c("tau", "sigma", "mu0"))
# the chains are not moxing well
rstan::check_hmc_diagnostics(stanrp)

print(stanrp, pars = c("tau", "sigma", "mu0", "phi", "p", "q"))

posterior_array <- as.array(stanrp)

# plot marginal posteriors
mcmc_areas(posterior_array,
           pars = c("mu0", "sigma0", "phi[1]", "phi[2]", "p[1]", "q[1]"),
           prob = 0.8)  # shades the 80% credible interval

# check for any parameters that look identical to their prior
# (suggests the data isn't informing that parameter)

# pair plots
mcmc_pairs(posterior_array,
           pars = c("tau[1]", "sigma[1]", "phi[1]"),
           off_diag_args = list(size = 0.5, alpha = 0.5))

# effect of distance

{
  # select only columns using distance coefficient. this is the same as t, as there are t-1 columns for each time point and
  # then one column for distance.
  distance_effects <- matrix_results[,grepl(paste0("beta_z\\[",t), texts) ]

  # on each column, calculate the 95% credible intervals & mean
  # t here means transpose to give one row per species
  beta_z_results <- distance_effects %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame
  beta_z_results$mean = colMeans(distance_effects)
  beta_z_results$species = species_names
  beta_z_results = arrange(beta_z_results, mean)
  # add detection numbers

  # connect species names to columns
  # only plot bats
  bats = c("Nyctalus","Pipistrellus spp.", "Plecotus spp.")
  betaz_bats = beta_z_results %>%
    filter(species %in% bats) %>%
    arrange(mean)
  # add detection results
  betaz_bats = betaz_bats %>%
    left_join(detection_summ, by = c("species" = "SpeciesID"))
  betaz_bats$species = factor(betaz_bats$species, levels = betaz_bats$species)

  #   change to Nyctalus spp.
  betaz_bats$species = recode(betaz_bats$species, "Nyctalus" = "Nyctalus spp.")

  distanceefrp = ggplot(betaz_bats, aes(x = species,
                                      y = mean,
                                      ymin = `2.5%`,
                                      ymax = `97.5%`)) +
    geom_errorbar(linewidth = 1) +
    geom_point(size = 4) +
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of distance from roost on \n DNA biomass') +
    theme_classic() +
    ylim(-6.1,4.25) +
    geom_text(aes(label = paste0("n=", detections)),
              y = 4.1,
              size = 5) +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(angle = 20, hjust = 0.6, vjust = 0.75)
          # axis.title =element_text(size = 18),
          # legend.text = element_text(size = 14),
    ) +
    ggtitle("Richmond Park")


  distanceefrp
  }

# now look at effect of time.
# effect of distance
t-1
visit_effects <- matrix_results[,grepl("beta_z\\[[1-3],", texts) ]
month_cols = c(
  "March" = "#440154",
  "May" = "#31688e",
  "June" = "#21908CFF",
  "July" = "#35b779",
  "August" = "#fde725",
  "September" = "#fbb03b",
  "October" = "#ff7f0e"
)


{
  # pattern <- sprintf("\\[%d,", 1)
  # pattern = paste(sprintf("\\[%d,", 2), collapse = "|")
  # selected <- grepl(pattern, texts) & grepl("beta_z", texts)
  # on each column, calculate the 95% credible intervals
  # beta_z_results <- matrix_results[,selected] %>%
  beta_z_results <- visit_effects %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame
  beta_z_results$mean = colMeans(visit_effects)

  # not including march as this is reference value
  Month = c("August", "September", "October")
  beta_z_results$Species <- rep(species_names, each = t - 1)
  beta_z_results$Month <- factor(rep(Month, times = n_distinct(species_names)), levels = Month)

  # just bats
  beta_z_results = beta_z_results %>%
    filter(Species %in% bats)
  beta_z_results$Species = factor(beta_z_results$Species,
                                  levels = c("Pipistrellus spp.", "Plecotus spp.", "Nyctalus"))

#   change to Nyctalus spp.
  beta_z_results$Species = recode(beta_z_results$Species, "Nyctalus" = "Nyctalus spp.")
  seasonefrp = ggplot(beta_z_results, aes(x = Species,
                                        y = mean,
                                        ymin = `2.5%`,
                                        ymax = `97.5%`,
                                        colour = Month)) +
    geom_errorbar(position = position_dodge(width = 0.7)) +
    geom_point(size = 3, position = position_dodge(width = 0.7)) +
    # scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen")) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # scale_x_continuous(breaks = 1:t, labels = Visits) +
    xlab('Species') +
    ylab('Effect of season on DNA biomass') +
    # facet_wrap(~Species, scales = 'free') +
    ylim(-6,4.2) +
    theme_classic() +
    scale_color_manual(values = month_cols, breaks = levels(beta_z_results$Month)) +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(angle = 20, hjust = 0.6, vjust = 0.75),
          # axis.title =element_text(size = 12),
          # legend.text = element_text(size = 12),
          strip.background = element_blank())
  # strip.text = element_text(size = 12)) +
  # ylim(-1,8) +

  seasonefrp
}

{

  # select only columns using distance coefficient. this is the same as t, as there are t-1 columns for each time point and
  # then one column for distance.
  hab_effects <- matrix_results[,grepl(paste0("beta_z\\[",t+1), texts) ]

  # on each column, calculate the 95% credible intervals & mean
  # t here means transpose to give one row per species
  beta_z_results <- hab_effects %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame
  beta_z_results$mean = colMeans(hab_effects)
  beta_z_results$species = species_names
  beta_z_results = arrange(beta_z_results, mean)
  # add detection numbers

  # connect species names to columns
  # only plot bats
  betaz_bats = beta_z_results %>%
    filter(species %in% bats) %>%
    arrange(mean)
  # add detection results
  betaz_bats = betaz_bats %>%
    left_join(detection_summ, by = c("species" = "SpeciesID"))
  betaz_bats$species = factor(betaz_bats$species, levels = c("Pipistrellus spp.", "Plecotus spp.", "Nyctalus"))

# change to nyctalus spp.
  betaz_bats$species = recode(betaz_bats$species, "Nyctalus" = "Nyctalus spp.")
  habefrp = ggplot(betaz_bats, aes(x = species,
                                 y = mean,
                                 ymin = `2.5%`,
                                 ymax = `97.5%`)) +
    geom_errorbar(linewidth = 1) +
    geom_point(size = 4) +
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of open habitat on DNA biomass') +
    ylim(-6.1,4.25) +
    theme_classic() +
    geom_text(aes(label = paste0("n=", detections)),
              y = 4.5,
              size = 5) +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(angle = 20, hjust = 0.6, vjust = 0.75),
          axis.title =element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.position = "top")


  habefrp

}
distanceefrp +habefrp + seasonefrp

# actual DNA biomass

logl_results = matrix_results %>%
  as.data.frame()

texts <- colnames(logl_results)
loglsel <- grepl('logl', texts)
# for each species, there are 94 estimates of logl
# this will be in the order of 1-19 for visit 1, visit 2 etc, or visit 1-4 for site 1, site 2 etc.
# it is the same format as the original data - sites_ecol, repeated for each species.

logl_results = logl_results[,loglsel] %>%
  apply(., 2, function(x) c(
    mean = mean(x),
    quantile(x, probs = c(0.025, 0.975)))) %>% t %>%
  as.data.frame()

logl_results$species = rep(species_names, each = nrow(sites_ecol))
#merge cfvarssub with loglresults (each needs to repeat)
logl_results$Visit = rep(sites_ecol$Visit, times = length(species_names))
logl_results$Point_ID = rep(sites_ecol$Point_ID, times = length(species_names))

# get unique site coordinates
sites_unique = siteinfo %>%
  group_by(Point_ID, Longitude, Latitude, dist_m) %>%
  filter(row_number() == 1)
# join longitude and latitude columns from sites ecol using point_id
logl_results2 = logl_results %>%
  left_join(sites_unique[,c("Point_ID", "Longitude", "Latitude", "dist_m")],
            by = 'Point_ID')
logl_sf = st_as_sf(logl_results2, coords = c("Longitude", "Latitude"), crs = 4326)

ggplot(data = logl_sf[logl_sf$species %in% bats,]) +
  geom_sf(aes(size = mean), alpha = 0.6, colour = 'grey') +
  theme_minimal() + facet_grid(species ~ Visit)

# relabel visit
logl_sf$Visit <- factor(
  logl_sf$Visit,
  levels = 1:4,
  labels = c("June", "July", "August", "September")
)

library(ggh4x)
logl_sf_b = logl_sf %>% filter(species %in% bats)
ggplot(data = logl_sf_b, aes(x = dist_m, y = mean, ymax = `97.5%`, ymin = `2.5%`,
                           fill = Visit,
                           colour = Visit)) +
  geom_smooth() +
  facet_wrap(~species) +
  # ggh4x::facet_grid2(species~Visit, labeller = labeller(species = species_names), scales = 'free_x', independent = "x") +
  labs(fill = 'Visit', colour = 'Visit', x = 'Distance from roost (m)',
       y = 'DNA Biomass') +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_classic() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
  )

# logl is a proxy for how much biomass is at the site. but we don't know amplification rates, we don't
# know shedding rate. so its on an abstract scale because we cannot actually measure how much DNA is in the area.
ggplot(data = logl_sf[logl_sf$species %in% bats,], aes(x = dist_m, y = mean, ymax = `97.5%`, ymin = `2.5%`,
                                                       fill = Visit,
                                                       colour = Visit)) +
  geom_smooth() +
  facet_wrap(~species) +
  labs(fill = 'Visit', colour = 'Visit', x = 'Distance from roost (m)',
       y = 'DNA Biomass') +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_classic() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  ggtitle('Model: Approximate Sampling, all species, Richmond Park')

# effect of open habitat
{

  # select only columns using distance coefficient. this is the same as t, as there are t-1 columns for each time point and
  # then one column for distance.
  hab_effects <- matrix_results[,grepl(paste0("beta_z\\[",t+1), texts) ]

  # on each column, calculate the 95% credible intervals & mean
  # t here means transpose to give one row per species
  beta_z_results <- hab_effects %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame
  beta_z_results$mean = colMeans(hab_effects)
  beta_z_results$species = species_names
  beta_z_results = arrange(beta_z_results, mean)
  # add detection numbers

  # connect species names to columns
  # only plot bats
  betaz_bats = beta_z_results %>%
    filter(species %in% bats) %>%
    arrange(mean)
  # add detection results
  betaz_bats = betaz_bats %>%
    left_join(detection_summ, by = c("species" = "SpeciesID"))
  betaz_bats$species = factor(betaz_bats$species, levels = betaz_bats$species)

  habef = ggplot(betaz_bats, aes(x = species,
                                 y = mean,
                                 ymin = `2.5%`,
                                 ymax = `97.5%`)) +
    geom_errorbar(linewidth = 1) +
    geom_point(size = 4) +
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of open habitat on DNA biomass') +
    theme_classic() +
    geom_text(aes(label = paste0("n=", detections)),
              y = 3,
              size = 5) +
    theme(axis.text = element_text(size = 18),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title =element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.position = "top") +
    ggtitle("Model: Approximate Sampling, all species, Richmond Park")


  habef

  }

