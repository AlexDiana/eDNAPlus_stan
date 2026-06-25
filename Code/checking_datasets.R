# quickly checking the difference between two files
library(dplyr); library(tidyr)
Sdnew = read.csv("Data/NatureAir_16SmamSpatial_ALL.csv")
siteinfo = read.csv("output/allsitevars.csv")

sdold = read.csv("Data/NatureAir_SpatialDATA.csv")



# convert Sd to longformß
sdnew_long = pivot_longer(Sdnew, cols = 15:ncol(Sdnew), names_to = 'Sample', values_to = 'read_count')
sdold_long = pivot_longer(sdold, cols = 15:ncol(sdold), names_to = 'Sample', values_to = 'read_count')


# lets look at rhinolophus detections only
new_rhin = sdnew_long %>% filter(grepl("Rhinolophus", Genus)) %>%
  filter(read_count > 0) %>%
  select(ASV_ID, Genus, Species, Sample, read_count)
old_rhin = sdold_long %>% filter(grepl("Rhinolophus", Genus)) %>%
  filter(read_count > 0) %>%
  select(ASV_IS, Genus, Species, Sample, read_count)


new_bat = sdnew_long %>% filter(grepl("Chiroptera", Order)) %>%
  filter(read_count > 0) %>%
  select(ASV_ID, Order, Genus, Species, Sample, read_count)
old_bat = sdold_long %>% filter(grepl("Chiroptera", Order)) %>%
  filter(read_count > 0) %>%
  select(ASV_ID, Order, Genus, Species, Sample, read_count)

bats = full_join(new_bat, old_bat, by = c("ASV_ID", "Order", "Genus", "Species", "Sample"), suffix = c("_new", "_old")) %>%
  mutate(read_count_new = ifelse(is.na(read_count_new), 0, read_count_new),
         read_count_old = ifelse(is.na(read_count_old), 0, read_count_old),
         same = ifelse(read_count_new == read_count_old, 1,0))

write.csv(bats, "output/bats_comparison.csv", row.names = FALSE)
# fix typo
# whenname
