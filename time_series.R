### TIME SERIES ###

# loading packages
suppressPackageStartupMessages({
  library("readr")      # For reading in the data
  library("tibble")     # For handling tidyverse tibble data classes
  library("tidyr")      # For tidying data 
  library("dplyr")      # For data manipulation 
  library("stringr")    # For string manipulation
  library("MASS")       # Functions/datasets for statistical analysis
  library("lubridate")  # For date manipulation
  library("ggplot2")    # For creating static visualizations
  library("scales")     # For formatting plots axis
  library("gridExtra")  # Creates multiple grid-based plots
})
'%!in%' <- function(x,y)!('%in%'(x,y))

df <- read.csv("nepal.csv")
df.gps <- read.csv("nepal_gps.csv")
df.gps <- df.gps %>%
  filter(Year >= 2009) %>%
  filter(In_silico_serotype != "ALTERNATIVE_ALIB_NT" & 
           In_silico_serotype != "SWISS_NT" &
           In_silico_serotype != "COVERAGE TOO LOW" &
           In_silico_serotype != "UNTYPABLE") %>%
  mutate(Disease_Status = ifelse(Clinical_manifestation == "CARRIAGE", 
                                 "Carriage", "Disease"))
pcv10_serotypes <- c('1', '4', '5', '6B', '7F', '9V', '14', '18C', '19F', '23F')
pcv13_serotypes <- c('3', '6A', '19A')
pcv15_serotypes <- c('22F', '33F')
pcv20_serotypes <- c('8', '10A', '11A', '12F', '15B')
df.gps$Serotype_group <- ifelse(
  df.gps$In_silico_serotype %in% pcv10_serotypes, "PCV10",
  ifelse(df.gps$In_silico_serotype %in% pcv13_serotypes, "PCV13",
  ifelse(df.gps$In_silico_serotype %in% pcv15_serotypes, "PCV15",
  ifelse(df.gps$In_silico_serotype %in% pcv20_serotypes, "PCV20",
  "OTHER"))))

df.gps2 <- df.gps %>%
  filter(Year >= 2009) %>%
  mutate(Disease_Status = ifelse(Clinical_manifestation == "CARRIAGE", "Carriage", "Disease")) %>%
  group_by(Year, Disease_Status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Disease_Status)
View(df.gps2)

p1 <- ggplot(df.gps2, aes(x = Year, y = count, color=Disease_Status)) +
      geom_line() +
      labs(title = "Cases of Pneumococcal Disease in Nepal Over time",
        x = "Year", y = "Count") +
      geom_vline(xintercept = 2015, col = 'red', lty = 2) +
      theme_linedraw()
p1
