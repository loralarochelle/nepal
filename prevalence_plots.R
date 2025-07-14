library(tidyverse)
nepal <- read.csv("nepal.csv")
nepal_gps <- read.csv("nepal_gps.csv")
pcv10_serotypes <- c('1', '4', '5', '6B', '7F', '9V', '14', '18C', '19F', '23F')
pcv13_serotypes <- c('3', '6A', '19A')
pcv15_serotypes <- c('22F', '33F')
pcv20_serotypes <- c('8', '10A', '11A', '12F', '15B')

serotypes_all <- unique(nepal_gps$In_silico_serotype)
non_vaccine_serotypes <- setdiff(serotypes_all, pcv10_serotypes)
non_vaccine_serotypes <- non_vaccine_serotypes[non_vaccine_serotypes != "ALTERNATIVE_ALIB_NT"
                                               & non_vaccine_serotypes != "COVERAGE TOO LOW" 
                                               & non_vaccine_serotypes != "SWISS_NT"]

nepal_shortened <- nepal_gps[, c("Year", "Clinical_manifestation", "In_silico_serotype")]

nepal_shortened$Serotype_group <- ifelse(
  nepal_shortened$In_silico_serotype %in% pcv10_serotypes, "PCV10",
  ifelse(nepal_shortened$In_silico_serotype %in% pcv13_serotypes, "PCV13",
  ifelse(nepal_shortened$In_silico_serotype %in% pcv15_serotypes, "PCV15",
  ifelse(nepal_shortened$In_silico_serotype %in% pcv20_serotypes, "PCV20",
  "OTHER"))))
View(nepal_shortened)

# breaking up data by year, serotype, and clinical manifestation
# excludes sparse data: clinical_manifestation
nepal_serotype_counts <- nepal_shortened %>%
  filter(Year >= 2009) %>%
  mutate(Disease_Status = ifelse(Clinical_manifestation == "CARRIAGE", "Carriage", "Disease")) %>%
  group_by(Year, Disease_Status, Serotype_group) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Disease_Status)

# carriage, meningitis, pneumonia plot over time - most prevalent. 
# no meningitis obs after 2013
ggplot(data=nepal_serotype_counts, aes(x=Year, y=count, color=Serotype_group)) +
  geom_line()+
  geom_point()+
  facet_grid(Disease_Status ~ .)+
  labs(title="Pneumococcal Disease Prevalence by Serotype Group over Time")+
  theme_bw()


ggplot(nepal_serotype_props, aes(x = Year, y = prop, color = Serotype_group)) +
  geom_line() +
  geom_point() +
  facet_grid(Disease_Status ~ .) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Proportion of Pneumococcal Serotype Groups by Year",
    y = "Proportion of Cases",
    x = "Year",
    color = "Serotype Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

