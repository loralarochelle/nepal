jcode <- " model {
  for (i in 1:N) {
    counts[i, 1:S] ~ dmulti(p[i, 1:S], total_cases[i])
    
    for (j in 1:S) {
      logit(p[i, j]) <- alpha[sero[j]] + beta[sero[j]] * (year_id[i] - mean_year)
    }
  }
  
  for (s in 1:S) {
    alpha[s] ~ dnorm(0, 0.001)
    beta[s] ~ dnorm(0, 0.001)
  }
  
  mean_year <- (n_year + 1) / 2
}"