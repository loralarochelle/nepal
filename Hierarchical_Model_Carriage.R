library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(lubridate)
# library(rgdal)
# library(maptools)
library(MASS)
library(rjags)
library(HDInterval)
library(ggplot2)
library(sf)
#library(spData)
library(areal)
library(ggthemes)
library(gridExtra)
library(RColorBrewer)
library(readr)
library(wesanderson)
library(viridis)
library(ggspatial)
library(janitor)

carriage_serotype_year_counts <- nepal_carriage %>%
  pivot_longer(cols = -Serotype, names_to = "Year", values_to = "cases")

carriage_serotype_year_counts <- carriage_serotype_year_counts %>% 
  filter(Year >= 2009 & Year != 2014)

# getting sums of counts for PCV13, PCV15, and PCV20 serotypes as rows as well 
df_sums <- carriage_serotype_year_counts %>%
  mutate(group = case_when(Serotype %in% pcv13_serotypes ~ "PCV13serotypes",
                           Serotype %in% pcv15_serotypes ~ "PCV15serotypes",
                           Serotype %in% pcv20_serotypes ~ "PCV20serotypes")) %>%
  filter(!is.na(group)) %>%
  group_by(Serotype = group, Year) %>%
  summarise(cases = sum(cases), .groups = "drop")
carriage_serotype_year_counts <- bind_rows(carriage_serotype_year_counts, df_sums)
View(carriage_serotype_year_counts)
# Assign numeric IDs for JAGS
carriage_serotype_year_counts <- carriage_serotype_year_counts %>%
  mutate(sero_id = as.integer(factor(Serotype)), # assigning IDs for serotypes
         year_id = as.integer(factor(Year))) # assigning IDs for the year

# Create JAGS data list
jdat <- list(
  N = nrow(carriage_serotype_year_counts), # number of rows = number of cases
  cases = carriage_serotype_year_counts$cases, # cases = cases
  sero_id = carriage_serotype_year_counts$sero_id, # serotype IDs
  year_id = carriage_serotype_year_counts$year_id, # year IDs
  n_sero = length(unique(carriage_serotype_year_counts$sero_id)), # number of unique serotypes
  n_year = length(unique(carriage_serotype_year_counts$year_id)) # number of years
)

library(rjags)
source("Model2.R")
model_path <- "Model2.R"

# creates jags model based on the textConnection to the jcode file (Model2), using jdat
mod <- jags.model(textConnection(jcode), data = jdat, n.chains = 2)
update(mod, 1000)  # burn-in
# generates posterior samples based on mu and beta, iterates 5000 times
samp <- coda.samples(mod, variable.names = c("mu", "beta"), n.iter = 5000) 
# summarizes samples: gives quantile for each variable
summary(samp)

library(coda)
par(mar = c(2, 2, 2, 2))  # smaller margins
# plot(samp)
# Density plots for parameters: these are what is random?
# densplot(samp, main = "Posterior Density for Parameters")

posterior_summary <- summary(samp)
# rounds estimates for each posterior value to 3 places
round(posterior_summary$statistics, 3)  # means, SDs
round(posterior_summary$quantiles, 3)   # 2.5%, 50%, 97.5%

# creates matrix of all betas across all simulation runs
beta_samples <- as.matrix(samp)[, grep("beta", colnames(as.matrix(samp)))]
# gives mean of each beta across all simulation runs
beta_means <- apply(beta_samples, 2, mean)
# makes CI for each beta
beta_ci <- apply(beta_samples, 2, quantile, probs = c(0.025, 0.975))

# gives mean and CI
df <- data.frame(
  param = colnames(beta_samples),
  mean = beta_means,
  lower = beta_ci[1, ],
  upper = beta_ci[2, ]
)

library(ggplot2)
ggplot(df, aes(x = param, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() +
  labs(title = "Posterior Estimates for Beta", y = "Effect Size", x = "") +
  theme_minimal()

#heatmap
# Your original year and serotype info
year_seq <- sort(unique(jdat$year_id))  # e.g., 2000:2006
sero_seq <- 1:jdat$n_sero
mean_year <- (length(year_seq) + 1) / 2

# Create a dataframe of all serotype-year combinations
pred_grid <- expand.grid(
  sero = sero_seq,
  year = year_seq
)

samp_mat <- as.matrix(samp)
mu_means <- colMeans(samp_mat[, grep("^mu\\[", colnames(samp_mat))])


# Predict expected log incidence
pred_grid$log_lambda <- mu_means[pred_grid$sero] +
  beta_means[pred_grid$sero] * (pred_grid$year - mean_year)

# Back-transform to incidence
pred_grid$lambda <- exp(pred_grid$log_lambda)

ggplot(pred_grid, aes(x = year, y = factor(sero), fill = lambda)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "Expected\nCases", option = "C") +
  labs(x = "Year", y = "Serotype", title = "Expected Cases by Year and Serotype") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

library(ggplot2)

ggplot(carriage_serotype_year_counts, aes(x = Year, y = Serotype, fill = cases)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "C", name = "Cases") +
  labs(title = "IPD Cases by Serotype and Year",
       x = "Year",
       y = "Serotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#model version
mean_year <- (max(carriage_serotype_year_counts$year_id) + 1) / 2

# Create a data frame with expected log lambda for each serotype-year
carriage_serotype_year_counts$expected_log_lambda <- mu_means[carriage_serotype_year_counts$sero_id] + 
  beta_means[carriage_serotype_year_counts$sero_id] * (carriage_serotype_year_counts$year_id - mean_year)

# Convert to expected cases (lambda)
carriage_serotype_year_counts$expected_cases <- exp(carriage_serotype_year_counts$expected_log_lambda)

#plot
ggplot(carriage_serotype_year_counts, aes(x = Year)) +
  geom_point(aes(y = cases), color = "blue", alpha = 0.5) +
  geom_line(aes(y = expected_cases, group = Serotype), color = "red") +
  facet_wrap(~ Serotype, scales = "free_y") +
  labs(y = "Cases", x = "Year", title = "Observed (points) vs Expected (lines) Cases by Serotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
