# Load required packages
library(randomForestSRC)

# Generate sample data for recurrent events without a terminating event
n <- 100  # Number of individuals

# Simulate recurrent event times and counts
recurrent_event_times <- rexp(n, rate = 0.1)
recurrent_event_counts <- rpois(n, lambda = 3)

# Simulate covariates (just an example, replace with your own)
covariate_1 <- rnorm(n)
covariate_2 <- rnorm(n)

# Create the data frame
data <- data.frame(
  Individual = rep(1:n, times = recurrent_event_counts),
  Recurrent_Event_Time = unlist(lapply(recurrent_event_counts, function(count) rexp(count, rate = 0.1))),
  Recurrent_Event_Count = rep(recurrent_event_counts),
  Covariate_1 = rep(covariate_1, times = recurrent_event_counts),
  Covariate_2 = rep(covariate_2, times = recurrent_event_counts)
)

# Fit the RSF model
rsf_model <- rfsrc(
  Surv(Recurrent_Event_Time, Recurrent_Event_Count) ~ Covariate_1 + Covariate_2,
  data = data,
  nsplit = 5,  # Number of splits for RSF
  ntree = 500  # Number of trees
)

# Display the RSF model summary
summary(rsf_model)
