# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(conflicted)
library(rjags)
library(coda)

disease_groups2 <- disease_groups %>%
  mutate(
    Serotype = case_when(
      Serotype %in% c("PCV13serotypes", "PCV15serotypes", "PCV20serotypes") ~ "Additional PCV",
      TRUE ~ Serotype
    )
  ) %>%
  group_by(Year, Serotype) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop")

# Step 1: Pivot data wide
df_wide <- disease_groups2 %>%
  pivot_wider(
    names_from = Serotype,
    values_from = cases,
    values_fill = 0
  )

df_wide$Year <- as.numeric(df_wide$Year)

# Step 2: Calculate total cases per year
df_wide <- df_wide %>%
  mutate(total_cases = rowSums(dplyr::select(., -Year))) %>%
  dplyr::select(Year, total_cases, everything())

# Step 3: Identify serotype columns
serotype_cols <- colnames(df_wide)[!(colnames(df_wide) %in% c("Year", "total_cases"))]
sero_id_map <- setNames(seq_along(serotype_cols), serotype_cols)

# Step 4: Reshape to long format for JAGS
df_long <- df_wide %>%
  pivot_longer(
    cols = all_of(serotype_cols),
    names_to = "Serotype",
    values_to = "cases"
  ) %>%
  mutate(
    sero_id = as.integer(factor(Serotype, levels = serotype_cols)),
    year_id = as.integer(factor(Year))
  )

jdat <- list(
  N = length(unique(df_long$Year)),        # one row per year
  S = length(unique(df_long$sero_id)),     # number of serotype groups
  counts = matrix(df_long$cases,           # reshape to matrix: rows = years, cols = serotypes
                  nrow = length(unique(df_long$Year)),
                  ncol = length(unique(df_long$sero_id)),
                  byrow = TRUE),
  total_cases = df_long$total_cases[seq(1, nrow(df_long), by = length(unique(df_long$sero_id)))],
  year_id = df_long$year_id[seq(1, nrow(df_long), by = length(unique(df_long$sero_id)))],
  mean_year = (length(unique(df_long$year_id)) + 1) / 2,
  sero_id = df_long$sero_id
)

jcode_og <- "model {
  for (i in 1:N) {
    counts[i, 1:S] ~ dmulti(p[i, 1:S], total_cases[i])
    
    for (j in 1:S) {
      logit(p[i, j]) <- alpha[sero_id[j]] + beta[sero_id[j]] * (year_id[i] - mean_year)
    }
  }
  
  for (s in 1:S) {
    alpha[s] ~ dnorm(0, 0.001)
    beta[s] ~ dnorm(0, 1)
  }
}" 

jcode <- "model {
  for (i in 1:N) {
    cases[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha[sero_id[i]] + beta[sero_id[i]] * (year_id[i] - mean_year) + log(total_cases[i])
  }
  
  for (s in 1:S) {
    alpha[s] ~ dnorm(0, 0.001)       # Intercept per serotype
    beta[s] ~ dnorm(0, 0.001)        # Slope per serotype
  }
}"


# Step 7: Run JAGS model
mod <- jags.model(textConnection(jcode), data = jdat, n.chains = 2)
update(mod, 1000)  # Burn-in
samp <- coda.samples(mod, variable.names = c("alpha", "beta"), n.iter = 5000)

# Step 8: Summarize and plot
posterior_summary <- summary(samp)
round(posterior_summary$statistics, 3)
round(posterior_summary$quantiles, 3)

# Extract beta samples
samp_mat <- as.matrix(samp)
beta_samples <- samp_mat[, grep("^beta", colnames(samp_mat))]
beta_means <- colMeans(beta_samples)
beta_ci <- apply(beta_samples, 2, quantile, probs = c(0.025, 0.975))

# Step 9: Plot beta estimates
df_beta <- data.frame(
  param = names(beta_means),
  mean = beta_means,
  lower = beta_ci[1, ],
  upper = beta_ci[2, ]
)

ggplot(df_beta, aes(x = param, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() +
  labs(title = "Posterior Estimates for Beta", y = "Effect Size", x = "") +
  theme_minimal()

# Step 10: Heatmap of predicted incidence
mu_means <- colMeans(samp_mat[, grep("^alpha\\[", colnames(samp_mat))])

year_seq <- sort(unique(jdat$year_id))
sero_seq <- 1:jdat$S

pred_grid <- expand.grid(
  sero = sero_seq,
  year = year_seq
)

pred_grid$log_lambda <- mu_means[pred_grid$sero] +
  beta_means[pred_grid$sero] * (pred_grid$year - jdat$mean_year)
# pred_grid$lambda <- exp(pred_grid$log_lambda)

total_cases_by_year <- disease_groups %>%
  group_by(Year) %>%
  summarise(total_cases = sum(cases))

year_totals <- total_cases_by_year$total_cases
names(year_totals) <- total_cases_by_year$Year
pred_grid$Year <- as.numeric(levels(factor(df_long$Year)))[pred_grid$year]  # map year_id back to Year
pred_grid$total_cases <- year_totals[as.character(pred_grid$Year)]
pred_grid$log_lambda <- mu_means[pred_grid$sero] +
  beta_means[pred_grid$sero] * (pred_grid$year - jdat$mean_year)
pred_grid$lambda <- exp(pred_grid$log_lambda + log(pred_grid$total_cases))

# Optional: Heatmap visualization (serotype vs year)
library(ggplot2)
ggplot(pred_grid, aes(x = year, y = sero, fill = lambda)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Expected\nIncidence") +
  labs(title = "Predicted Incidence by Serotype and Year", x = "Year ID", y = "Serotype ID") +
  theme_minimal()

# Step 1: Recreate or confirm serotype group names
serotype_labels <- names(sero_id_map)  # original names like PCV10serotypes, etc.

# Add serotype group names back to df_long
df_long <- df_long %>%
  mutate(Serotype = factor(Serotype, levels = serotype_labels))

# Step 2: Compute total cases per year to get proportions
df_long <- df_long %>%
  group_by(Year) %>%
  mutate(year_total = sum(cases)) %>%
  ungroup() %>%
  mutate(Percent = (cases / year_total) * 100)

# Step 3: Get expected percent from model output
# Get mu and beta means
mu_means <- colMeans(samp_mat[, grep("^alpha\\[", colnames(samp_mat))])
beta_means <- colMeans(samp_mat[, grep("^beta\\[", colnames(samp_mat))])

# Add sero_id and year_id back (if not already in df_long)
df_long <- df_long %>%
  mutate(
    sero_id = as.integer(factor(Serotype, levels = serotype_labels)),
    year_id = as.integer(factor(Year)),
    mean_year = (length(unique(year_id)) + 1) / 2
  )

# Expected log incidence and proportion (rescaled to %)
df_long <- df_long %>%
  mutate(
    log_lambda = mu_means[sero_id] + beta_means[sero_id] * (year_id - mean_year),
    lambda = exp(log_lambda)
  ) %>%
  group_by(Year) %>%
  mutate(expected_cases = lambda / sum(lambda) * 100) %>%
  ungroup()

df_long2 <- df_long %>%
  mutate(Serotype = dplyr::recode(Serotype,
                                  "PCV10serotypes" = "PCV10 Serotypes",
                                  "PCV13serotypes" = "Additional PCV13",
                                  "PCV15serotypes" = "Additional PCV15",
                                  "PCV20serotypes" = "Additional PCV20",
                                  "Other" = "Non-PCV Serotypes"),
  Serotype = factor(Serotype, levels = c(
    "PCV10 Serotypes",
    "Additional PCV13",
    "Additional PCV15",
    "Additional PCV20",
    "Non-PCV Serotypes"
  )))

# Step 5: Plot observed vs expected
custom_colors <- c(
  "PCV10 Serotypes"    = "#00356B",
  "Additional PCV13"   = "#6A0DAD",
  "Additional PCV15"   = "#008080",
  "Additional PCV20"   = "#E69F00",
  "Non-PCV Serotypes"  = "#999999"
)

disease_plot <- ggplot(df_long, aes(x = Year, group = Serotype, color = Serotype)) +
  geom_point(aes(y = Percent), alpha = 0.6) +
  geom_line(aes(y = expected_cases), size = 1) +
  scale_x_continuous(
    breaks = seq(floor(min(df_long2$Year)), ceiling(max(df_long2$Year)), by = 2)) +
  # scale_color_manual(values=custom_colors)+
  labs(
    x = "Year",
    y = "Percent of Cases Per Year",
    title = "Observed vs Expected Serotype Group Trends",
    color = "Serotype Group"
  ) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 15, hjust = 0.5),
  )
disease_plot
