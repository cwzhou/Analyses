# Set parameters
n <- 3 # Sample size
a <- 0  # Lower limit of censoring time
b <- 5  # Upper limit of censoring time
# Define separate beta values for cause 1 and cause 2
beta_cause1 <- c(0.1, 0.2)  # Example beta values for cause 1
beta_cause2 <- c(0.321, 4.2)  # Example beta values for cause 2

# Generate covariates Zi = (Zi1, Zi2)
Zi1 <- rnorm(n, mean = 0, sd = 1)
Zi2 <- rnorm(n, mean = 0, sd = 1)
Z = cbind(Zi1,Zi2)

# Generate subdistribution for type 1 failures
p <- 0.2  # Parameter p
T1_subdist <- numeric(n)
for (i in 1:n) {
  # t <- runif(1, 0, 1)
  if (Zi1[i] == 0 && Zi2[i] == 0) {
    T1_subdist[i] <- rexp(1, 1)
  } else {
    # cump1 <- 1 - (1 - p) * exp(-t)
    hazard_ratio1 <- as.numeric(exp(Z %*% beta_cause1))
    # T1_subdist[i] <- -log(1 - runif(1) * (1 - cump1)^hazard_ratio) / hazard_ratio
    num = log(1-runif(1))
    denom = hazard_ratio1
    above0 = exp(num/denom)
    above = 1-above0
    T1_subdist[i] <- -log(1-above/p)
    print(T1_subdist)
  }
}

# Generate subdistribution for type 2 failures
T2_subdist <- numeric(n)
for (i in 1:n) {
  hazard_ratio1 <- as.numeric(exp(Z %*% beta_cause1))
  hazard_ratio2 = as.numeric(exp(Z %*% beta_cause2))
  # Calculate the conditional probability for type 2 failures
  # Pr_E2 <- 1 - (1 - p)^exp(Zi1[i] * 0.1 + Zi2[i] * 0.2)
  Pr_E2 <- (1-p)^(hazard_ratio1)
  # Generate random values following the conditional distribution
  T2_subdist[i] <- -log(1 - runif(1) / Pr_E2)/ hazard_ratio2
  print(T2_subdist)
}

# Generate event times (minimum of T1 and T2 subdistributions)
event_time <- pmin(T1_subdist, T2_subdist)

# Generate censoring times from the uniform [a, b] distribution
censor_time <- runif(n, a, b)

# Merge event and censoring times
time_final <- pmin(event_time, censor_time)

# Generate event/censoring indicators
event_censor <- ifelse(event_time < censor_time, 1, 0)

# Create a survival object for analysis
survival_obj <- Surv(time_final, event_censor)

# Further analysis or modeling using the simulated data
print(event_time)
