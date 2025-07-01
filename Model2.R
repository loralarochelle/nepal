# weakly informative priors
# made using knowledge of the space
# dpois is frequency of each serotype going in: basing this on count of each serotypes
# assigning prior from normal to make a model for IPD from each serotype using predictors of time
# and the matrix of all the serotypes, cases per year of each, etc.
# 83 outputs from each model (one for each serotype), averages across mu1 for serotype 1 from all 
# simulation runs. adds extra variation to give a more sound estimate
# to ask How many cases from serotype 1? You're asking the whole vector of 83 which are kinda influenced
# by the rest of it
# running w fewer serotypes: probably worse fits? added stability w more serotypes
# looking only at vaccine serotypes: possibly an issue, pulls towards overall average if there are fewer
# point estimates, w sparse counts you don't want the estimate at the bounds of its CI, you want to pull
# it closer to a conservative guess (create more narrow, plausible range)

jcode <- "
model {
  for (i in 1:N) {
    cases[i] ~ dpois(lambda[i])
    log(lambda[i]) <- mu[sero_id[i]] + beta[sero_id[i]] * (year_id[i] - mean_year)
  }

  for (s in 1:n_sero) {
    mu[s] ~ dnorm(0, 0.001)
    beta[s] ~ dnorm(0, 0.001)
  }

  mean_year <- (n_year + 1) / 2
}
"
