# Set the parameters
rho_values <- c(0, 0.0625, 0.125, 0.25); rho = rho_values[1]
alpha_values = rho_values*4; alpha = alpha[1]
lambda_0D <- 0.25
lambda_0R <- 1
n_values <- c(50, 100); n = n_values[1]
num_samples <- 10000; num_samples = 1
tau = 10

# Parameters and Covariates
beta_D1 = 1
beta_D2 = 2
beta_D = c(beta_D1, beta_D2)
beta_R1 = 1.1
beta_R2 = 2.1
beta_R = c(beta_R1, beta_R2)
# Sample covariate cov1 from N(0,1)
cov1 = rnorm(n)
# Sample covariate cov2 from a Bernoulli distribution
cov2 = rbinom(n, 1, 0.5)
cov = cbind(cov1,cov2)

# Create marginal exponential rates
lambda_D = lambda_0D*exp(cov%*%beta_D)
lambda = lambda_0R*exp(cov%*%beta_R)

# Simulate the scenario
results <- list()

for (rho in rho_values) {
  print(rho)
  for (n in n_values) {
    print(n)
    for (i in 1:num_samples) {
      if (i %% 25 == 0) {
        cat("Iteration:", i, "\n")
      }
      # Generate data
      censoring_time <- runif(n, min = 0, max = tau)
      
      survival_time <- rexp(n, rate = lambda_D)
      gap_time1 <- rexp(n, rate = lambda)
      recurrent_times <- cumsum(c(0, gap_time1))
      R1 = recurrent_gaps[1]
      # Create dependence among variables using Gumbel's bivariate exponential distribution
      # Correlation between the first recurrent time and survival time
      u1 <- runif(n)
      
      u2 <- runif(n)
      x1 <- -log(-log(u1))
      x2 <- -log(-log(u2))
      v1 <- (x1 - log(1 - rho)) / sqrt(-2 * log(1 - rho))
      v2 <- (x2 - log(1 - rho)) / sqrt(-2 * log(1 - rho))
      y1 <- exp(-exp(-v1))
      y2 <- exp(-exp(-v2))
      recurrent_times[1] <- y1[1] * survival_time[1]
      survival_time <- y1 * survival_time
      recurrent_times[-1] <- recurrent_times[-1] + y2[-1] * survival_time[-1]
      

      observed_events <- pmin(recurrent_times, censoring_time)
    }
  }
}

                        