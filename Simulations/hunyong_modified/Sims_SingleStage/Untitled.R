backsolve_t1 <- function(u1i, p, Zi, beta1) {
  # Function to solve for t
  solve_t <- function(t, u1i, p, Zi, beta1) {
    (1 - (1 - p * (1 - exp(-t))) ^ exp((Zi %*% beta1))) / (1 - (1 - p) ^ exp(Zi %*% beta1)) - u1i
  }
  # Use uniroot to solve for t
  result <- uniroot(solve_t, interval = c(0, 100), u1i = u1i, p = p, Zi = Zi, beta1 = beta1)
  return(result$root)
}

# Define input values
N = 2000
ncov = 5
# seed1 = (N+ncov)*0.2
# set.seed(seed1)
u1 <- runif(N)  # Observed value
p <- 0.4  # Probability
Z <- matrix(rnorm(N * ncov), nrow = N, ncol = ncov)  # Covariate vector
beta1 <- c(-0.1,-0.1,0.1,0.1,-0.4)  # Beta coefficients
beta2 <- c(0.1,0.5,-0.1,0.1,0.1)  # Beta coefficients
#censoring
a = 0
b = 1
failure_t1 = numeric(0)
for (i in 1:N){
  u1i = u1[i]
  Zi = Z[i,]
  # Use backsolve_t function
  failure_t1[i] <- backsolve_t1(u1i, p, Zi, beta1)
  # message("id:",i,"\nfailure_t1:", failure_t1)
}
failure_t2_rate = exp(Z %*% beta2)
failure_t2 = rexp(N, rate = failure_t2_rate)
# message("id:",i,"\nfailure_t2:",failure_t2)

# censoring
C = runif(N, min = a, max = b)

data = data.frame(ID = 1:N, 
                failure1_time = failure_t1,
                failure2_time = failure_t2,
                obs_failureCR_time = pmin(failure_t1, failure_t2),
                censoring_time = C
                )

# Calculate obs_time as the minimum of failure times and censoring time
data$obs_time <- do.call(pmin, c(data[c("failure1_time", 
                                        "failure2_time", 
                                        "censoring_time")], na.rm = TRUE))

# Create status variable based on obs_time and different failure times
data$status <- ifelse(data$obs_time == data$failure1_time, 1, 
                      ifelse(data$obs_time == data$failure2_time, 2, 
                             0))

type1FailureRate = sum(data$status == 1)/N
type2FailureRate = sum(data$status == 2)/N
censoringRate = sum(data$status == 0)/N
message("type 1 failure rate: ", type1FailureRate)
message("type 2 failure rate: ", type2FailureRate)
message("censoring rate: ", censoringRate)
View(data)
data = data %>% dplyr::select(ID, obs_time, status)
