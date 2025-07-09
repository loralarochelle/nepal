## fixing

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
  filter(!(In_silico_serotype %in% c("COVERAGE TOO LOW", "UNTYPABLE", "SWISS_NT", "ALTERNATIVE_ALIB_NT"))) %>%
  # combine non-PCV serotypes into "other" category
  mutate(In_silico_serotype = if_else(In_silico_serotype %in% all_pcv_serotypes, 
                                      In_silico_serotype, "Other")) %>%
  # group by year, serotype, and count
  group_by(Year, In_silico_serotype) %>%
  dplyr::summarize(count = n(), .groups = "drop") %>%
  rename(Serotype = In_silico_serotype) %>%
  # pivot into same form as nepal dataset: years as cols, values as counts, all NAs are 0
  pivot_wider(names_from = Year, values_from = count, values_fill = 0)

# select only counts from the nepal dataset (not proportions), fix name to be just an int
nepal_new <- nepal %>%
  dplyr::select(Serotype, X2014n, X2015n, X2017n, X2018n, X2019n, X2021n) %>%
  rename_with(~ gsub("X|n", "", .x), starts_with("X")) %>%
  filter(!(Serotype == "PCV10serotypes"))

# join datasets together using full_join by serotype to get full carriage dataset
nepal_carriage <- full_join(nepal_new, nepal.gps.carriage, by = "Serotype") %>%
  mutate(across(-Serotype, ~ replace_na(.x, 0))) %>%
  # get rid of Total row: it is no longer correct
  filter(Serotype != "Total") %>%
  # get rid of .x and .y from years that only had data from one of the two datasets
  mutate(across(ends_with(".x"), ~ coalesce(.x, get(sub(".x$", ".y", cur_column()))),
                .names = "{sub('.x$', '', .col)}")) %>%
  dplyr::select(-matches("\\.x$|\\.y$"))

carriage_serotype_year_counts <- nepal_carriage %>%
  pivot_longer(cols = -Serotype, names_to = "Year", values_to = "cases")

df_sums <- carriage_serotype_year_counts %>%
  mutate(group = case_when(Serotype %in% pcv10_serotypes ~ "PCV10serotypes",
                           Serotype %in% pcv13_serotypes ~ "PCV13serotypes",
                           Serotype %in% pcv15_serotypes ~ "PCV15serotypes",
                           Serotype %in% pcv20_serotypes ~ "PCV20serotypes")) %>%
  filter(!is.na(group)) %>%
  group_by(Serotype = group, Year) %>%
  summarise(cases = sum(cases), .groups = "drop")
carriage_groups_df <- bind_rows(carriage_serotype_year_counts, df_sums)

carriage_groups <- carriage_groups_df %>%
  filter(Serotype %in% c("Other", "PCV10serotypes", "PCV13serotypes", "PCV15serotypes", "PCV20serotypes"))

View(carriage_groups)
carriage_groups

# get total cases per serotype per year
carriage_serotype_year_totals <- carriage_groups %>%
  group_by(Year, Serotype) %>%
  summarise(total_cases = sum(cases), .groups = "drop")
# get total cases per year
carriage_year_totals <- carriage_serotype_year_totals %>%
  group_by(Year) %>%
  summarise(year_total = sum(total_cases), .groups = "drop")
# join & compute proportions, then multiply by 100 and round down
carriage_serotype_proportions <- carriage_serotype_year_totals %>%
  left_join(carriage_year_totals, by = "Year") %>%
  mutate(proportion = total_cases / year_total,
         percent = floor(proportion * 100))
