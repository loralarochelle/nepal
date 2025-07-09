### Tree metatdata
library(tidyverse)
sample_isolates <- read.csv("diverse_isolates_sample.csv")

sample_isolates$macrolide_resistant <- ifelse(sample_isolates$ermB == "NEG" & sample_isolates$mefA == "NEG", "macrolide", "macrolide_resistant")
pcv10_serotypes <- c('1', '4', '5', '6B', '7F', '9V', '14', '18C', '19F', '23F') 
pcv13_serotypes <- c('3', '6A', '19A')
pcv15_serotypes <- c('22F', '33F')
pcv20_serotypes <- c('8', '10A', '11A', '12F', '15B')
untypable = c("ALTERNATIVE_ALIB_NT","COVERAGE TOO LOW", "SWISS_NT", "UNTYPABLE")
sample_isolates$Serotype_Group <- ifelse(sample_isolates$In_silico_serotype %in% pcv10_serotypes, "PCV10",
                                   ifelse(sample_isolates$In_silico_serotype %in% pcv13_serotypes, "PCV13",
                                          ifelse(sample_isolates$In_silico_serotype %in% pcv15_serotypes, "PCV15",
                                                 ifelse(sample_isolates$In_silico_serotype %in% pcv20_serotypes, "PCV20", 
                                                        ifelse(sample_isolates$In_silico_serotype %in% untypable, "UNTYPABLE",
                                                               "Other")))))
sample_isolates$PCV_Period <- ifelse(sample_isolates$Vaccine_period == "PREPCV", "PREPCV", "POSTPCV")
sample_isolates <- sample_isolates %>%
  select(Public_name, Age_years, Clinical_manifestation, PCV_Period, GPSC, WGS_PEN_SIR_Meningitis, Cot, folA_I100L, cat, macrolide_resistant, Serotype_Group)
cols_to_change1 <- c("Cot", "folA_I100L", "cat")
for (col in cols_to_change1) {
  sample_isolates[[col]] <- ifelse(sample_isolates[[col]] == "POS",
                             paste0("POS", col),
                             paste0("NEG", col))}

write.csv(sample_isolates, "tree_metadata.csv", row.names = FALSE)
