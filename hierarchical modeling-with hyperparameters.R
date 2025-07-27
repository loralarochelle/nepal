# hierarchical w hyperparameters

library(dplyr)
library(tidyr)
library(ggplot2)
library(rjags)
library(coda)

# 1. Prepare data wide
df_wide <- carriage_groups %>%
  pivot_wider(
    names_from = Serotype,
    values_from = cases,
    values_fill = 0
  )
df_wide$Year <- as.numeric(df_wide$Year)

# 2. Calculate total cases per year
df_wide <- df_wide %>%
  mutate(total_cases = rowSums(dplyr::select(., -Year))) %>%
  select(Year, total_cases, everything())

# 3. Identify serotype columns and map
serotype_cols <- colnames(df_wide)[!(colnames(df_wide) %in% c("Year", "total_cases"))]
sero_id_map <- setNames(seq_along(serotype_cols), serotype_cols)

# 4. Reshape to long format for JAGS
df_long3 <- df_wide %>%
  pivot_longer(
    cols = all_of(serotype_cols),
    names_to = "Serotype",
    values_to = "cases"
  ) %>%
  mutate(
    sero_id = as.integer(factor(Serotype, levels = serotype_cols)),
    year_id = as.integer(factor(Year))
  )

# 5. Prepare JAGS data list
jdat <- list(
  N = length(unique(df_long3$Year)),           # number of years
  S = length(unique(df_long3$sero_id)),        # number of serotypes
  counts = matrix(df_long3$cases,
                  nrow = length(unique(df_long3$Year)),
                  ncol = length(unique(df_long3$sero_id)),
                  byrow = TRUE),
  total_cases = df_long3$total_cases[seq(1, nrow(df_long3), by = length(unique(df_long3$sero_id)))],
  year_id = df_long3$year_id[seq(1, nrow(df_long3), by = length(unique(df_long3$sero_id)))],
  mean_year = (length(unique(df_long3$year_id)) + 1) / 2
)

# 6. JAGS model with hierarchy
jcode4 <- "model {
  for (i in 1:N) {
    counts[i, 1:S] ~ dmulti(p[i, 1:S], total_cases[i])
    
    for (j in 1:S) {
      eta[i, j] <- alpha[j] + beta[j] * (year_id[i] - mean_year)
    }
    
    for (j in 1:S) {
      exp_eta[i, j] <- exp(eta[i, j])
    }
    
    for (j in 1:S) {
      p[i, j] <- exp_eta[i, j] / sum(exp_eta[i, 1:S])
    }
  }
  
  for (j in 1:S) {
    alpha[j] ~ dnorm(mu_alpha, 0.001)
    beta[j]  ~ dnorm(mu_beta, 0.001)
  }
  mu_alpha ~ dnorm(0, 0.001)
  mu_beta  ~ dnorm(0, 0.001)
}
"


# 7. Run model
mod <- jags.model(textConnection(jcode4), data = jdat, n.chains = 2)
update(mod, 1000)  # Burn-in
samp <- coda.samples(mod, variable.names = c("alpha", "beta"), n.iter = 5000)

# 8. Posterior summaries
posterior_summary <- summary(samp)
round(posterior_summary$statistics, 3)
round(posterior_summary$quantiles, 3)

samp_mat <- as.matrix(samp)
beta_samples <- samp_mat[, grep("^beta", colnames(samp_mat))]
beta_means <- colMeans(beta_samples)
beta_ci <- apply(beta_samples, 2, quantile, probs = c(0.025, 0.975))

alpha_samples <- samp_mat[, grep("^alpha", colnames(samp_mat))]
alpha_means <- colMeans(alpha_samples)

# 9. Plot beta estimates
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

# 10. Calculate predicted probabilities with softmax for plotting

year_seq <- sort(unique(jdat$year_id))
sero_seq <- 1:jdat$S

pred_grid <- expand.grid(
  year = year_seq,
  sero = sero_seq
)

pred_grid <- pred_grid %>%
  group_by(year) %>%
  mutate(
    eta = alpha_means[sero] + beta_means[sero] * (year - jdat$mean_year),
    exp_eta = exp(eta),
    sum_exp_eta = sum(exp_eta),
    p = exp_eta / sum_exp_eta
  ) %>%
  ungroup()

# Optional: convert year_id back to actual Year for plotting (if desired)
year_map <- data.frame(year_id = year_seq, Year = sort(unique(df_long3$Year)))
pred_grid <- left_join(pred_grid, year_map, by = c("year" = "year_id"))

# 11. Add serotype labels back to df_long3 for plotting
serotype_labels <- names(sero_id_map)

df_long3 <- df_long3 %>%
  mutate(Serotype = factor(Serotype, levels = serotype_labels)) %>%
  group_by(Year) %>%
  mutate(
    eta = alpha_means[sero_id] + beta_means[sero_id] * (year_id - jdat$mean_year),
    exp_eta = exp(eta),
    sum_exp_eta = sum(exp_eta),
    expected_cases = (exp_eta / sum_exp_eta) * 100,  # percent expected cases
    Percent = (cases / sum(cases)) * 100            # observed percent cases
  ) %>%
  ungroup()

# 12. Rename serotype groups for plotting
df_long3 <- df_long3 %>%
  mutate(Serotype = recode(Serotype,
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

# 13. Plot observed vs expected
custom_colors <- c(
  "PCV10 Serotypes"    = "#00356B",
  "Additional PCV13"   = "#6A0DAD",
  "Additional PCV15"   = "#008080",
  "Additional PCV20"   = "#E69F00",
  "Non-PCV Serotypes"  = "#999999"
)

carriage_plot <- ggplot(df_long3, aes(x = Year, group = Serotype, color = Serotype)) +
  geom_point(aes(y = Percent), alpha = 0.6) +
  geom_line(aes(y = expected_cases), size = 1) +
  scale_x_continuous(
    breaks = seq(floor(min(df_long3$Year)), ceiling(max(df_long3$Year)), by = 2),
    limits = c(2009, 2021)
  ) +
  scale_color_manual(values=custom_colors) +
  ylim(0, 80) +
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
    plot.title = element_text(size = 15, hjust = 0.5)
  )

carriage_plot

carriage_plot <- ggplot(df_long3, aes(x = Year, group = Serotype, color = Serotype)) +
  geom_point(aes(y = Percent), alpha = 0.6) +
  geom_line(aes(y = expected_cases), size = 1) +
  scale_x_continuous(
    breaks = seq(floor(min(df_long3$Year)), ceiling(max(df_long3$Year)), by = 2),
    limits = c(2009, 2021)
  ) +
  scale_color_manual(values = custom_colors) +
  ylim(0, 80) +
  labs(
    x = "Year",
    y = "Percent of Cases Per Year",
    title = "Hierarchical Modeling: Percent of Carriage Cases",
    subtitle = expression(
      "For year " * t * " and serotype group " * s * 
        ":  " * 
        p[ts] == frac(exp(alpha[s] + beta[s] * (Year[t] - bar(Year))), 
                      sum(exp(alpha[k] + beta[k] * (Year[t] - bar(Year)))))
    )
    ,
    color = "Serotype Group"
  ) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, color = "purple", hjust = 0.5)
  )

ggsave("Carriage_Plot_HM.jpeg", plot=carriage_plot, width=9, height=5)
