#### COMBINING ALL CARRIAGE DATA - SEROTYPE LEVEL ####

# read in data
nepal <- read.csv("nepal.csv")
nepal.gps <- read.csv("nepal_gps.csv")

# make list of serotypes included in each vaccine
pcv10_serotypes <- c('1', '4', '5', '6B', '7F', '9V', '14', '18C', '19F', '23F')
pcv13_serotypes <- c('3', '6A', '19A')
pcv15_serotypes <- c('22F', '33F')
pcv20_serotypes <- c('8', '10A', '11A', '12F', '15B')
all_pcv_serotypes <- c(pcv10_serotypes, pcv13_serotypes, pcv15_serotypes, pcv20_serotypes)

# putting carriage samples from gps dataset into non-gps dataset form
nepal.gps.carriage <- nepal.gps %>%
  # take only carriage data from 2009 on
  filter(Clinical_manifestation == "CARRIAGE") %>%
  filter(Year >= 2009) %>%
  # filter out irrelevant serotypes
  filter(!(In_silico_serotype %in% c("COVERAGE TOO LOW", "UNTYPABLE", "SWISS_NT"))) %>%
  # combine non-PCV serotypes into "other" category
  mutate(In_silico_serotype = if_else(In_silico_serotype %in% all_pcv_serotypes, 
                                      In_silico_serotype, "Other")) %>%
  # group by year, serotype, and count
  group_by(Year, In_silico_serotype) %>%
  summarize(count = n(), .groups = "drop") %>%
  rename(Serotype = In_silico_serotype) %>%
  # pivot into same form as nepal dataset: years as cols, values as counts, all NAs are 0
  pivot_wider(names_from = Year, values_from = count, values_fill = 0)

# select only counts from the nepal dataset (not proportions), fix name to be just an int
nepal <- nepal %>%
  dplyr::select(Serotype, X2014n, X2015n, X2017n, X2018n, X2019n, X2021n) %>%
  rename_with(~ gsub("X|n", "", .x), starts_with("X"))

# join datasets together using full_join by serotype to get full carriage dataset
nepal_carriage <- full_join(nepal, nepal.gps.carriage, by = "Serotype") %>%
  mutate(across(-Serotype, ~ replace_na(.x, 0))) %>%
  # get rid of Total row: it is no longer correct
  filter(Serotype != "Total") %>%
  # get rid of .x and .y from years that only had data from one of the two datasets
  mutate(across(ends_with(".x"), ~ coalesce(.x, get(sub(".x$", ".y", cur_column()))),
                .names = "{sub('.x$', '', .col)}")) %>%
  dplyr::select(-matches("\\.x$|\\.y$"))

#### GETTING JUST DISEASE DATA FROM GPS DATASET ####
nepal_disease <- nepal.gps %>%
  filter(Clinical_manifestation != "CARRIAGE") %>%
  mutate(Clinical_manifestation = "DISEASE") %>% # 
  filter(Year >= 2009) %>%
  # filter out irrelevant serotypes
  filter(!(In_silico_serotype %in% c("COVERAGE TOO LOW", "UNTYPABLE", "SWISS_NT"))) %>%
  # combine non-PCV serotypes into "other" category
  mutate(In_silico_serotype = if_else(In_silico_serotype %in% all_pcv_serotypes, 
                                      In_silico_serotype, "Other")) %>%
  rename(Serotype = In_silico_serotype) %>%
  group_by(Year, Serotype) %>%
  summarize(count = n(), .groups = "drop")
