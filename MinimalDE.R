library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(lubridate)
#library(rgdal)
#library(maptools)
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
#install.packages("rgdal")
d1.adults<-read.csv("./testDE2.csv")


# Aggregate data
serotype_year_counts <- d1.adults %>%
  count(SeroType, epiyr, name = "cases") %>%
  complete(SeroType, epiyr, fill = list(cases = 0))

# Assign numeric IDs for JAGS
serotype_year_counts <- serotype_year_counts %>%
  mutate(sero_id = as.integer(factor(SeroType)),
         year_id = as.integer(factor(epiyr)))

# Create JAGS data list
jdat <- list(
  N = nrow(serotype_year_counts),
  cases = serotype_year_counts$cases,
  sero_id = serotype_year_counts$sero_id,
  year_id = serotype_year_counts$year_id,
  n_sero = length(unique(serotype_year_counts$sero_id)),
  n_year = length(unique(serotype_year_counts$year_id))
)

library(rjags)
source("Model2.R")
model_path <- "Model2.R"

mod <- jags.model(textConnection(jcode), data = jdat, n.chains = 2)
update(mod, 1000)  # burn-in
samp <- coda.samples(mod, variable.names = c("mu", "beta"), n.iter = 5000)

summary(samp)

library(coda)
par(mar = c(2, 2, 2, 2))  # smaller margins
plot(samp)
# Density plots
densplot(samp, main = "Posterior Density for Parameters")

posterior_summary <- summary(samp)
round(posterior_summary$statistics, 3)  # means, SDs
round(posterior_summary$quantiles, 3)   # 2.5%, 50%, 97.5%

beta_samples <- as.matrix(samp)[, grep("beta", colnames(as.matrix(samp)))]
beta_means <- apply(beta_samples, 2, mean)
beta_ci <- apply(beta_samples, 2, quantile, probs = c(0.025, 0.975))

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
mu_means<- colMeans(samp_mat[, grep("^mu\\[", colnames(samp_mat))])

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

ggplot(serotype_year_counts, aes(x = epiyr, y = SeroType, fill = cases)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "C", name = "Cases") +
  labs(title = "IPD Cases by Serotype and Year",
       x = "Year",
       y = "Serotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#model version
mean_year <- (max(serotype_year_counts$year_id) + 1) / 2

# Create a data frame with expected log lambda for each serotype-year
serotype_year_counts$expected_log_lambda <- mu_means[serotype_year_counts$sero_id] + 
  beta_means[serotype_year_counts$sero_id] * (serotype_year_counts$year_id - mean_year)

# Convert to expected cases (lambda)
serotype_year_counts$expected_cases <- exp(serotype_year_counts$expected_log_lambda)

#plot
ggplot(serotype_year_counts, aes(x = epiyr)) +
  geom_point(aes(y = cases), color = "blue", alpha = 0.5) +
  geom_line(aes(y = expected_cases, group = SeroType), color = "red") +
  facet_wrap(~ SeroType, scales = "free_y") +
  labs(y = "Cases", x = "Year", title = "Observed (points) vs Expected (lines) Cases by Serotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))