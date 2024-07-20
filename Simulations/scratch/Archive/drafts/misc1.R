# Load required packages
library(randomForestSRC)

# Set random seed for reproducibility
set.seed(123)

# Generate sample data in the format with recurrent and terminal events
n <- 100  # Number of individuals

# Simulate recurrent event times and counts
recurrent_event_times <- rexp(n, rate = 0.1)
recurrent_event_counts <- rpois(n, lambda = 3)

# Simulate terminal event times and indicators (0: censored, 1: event occurred)
terminal_event_times <- rexp(n, rate = 0.05)
terminal_event_indicators <- rbinom(n, size = 1, prob = 0.8)

# Simulate covariates (just an example, replace with your own)
covariate_1 <- rnorm(n)
covariate_2 <- rnorm(n)

# Create the data frame
data <- data.frame(
  Individual = rep(1:n, times = recurrent_event_counts),
  Recurrent_Event_Time = unlist(lapply(recurrent_event_counts, function(count) rexp(count, rate = 0.1))),
  Recurrent_Event_Count = rep(recurrent_event_counts),
  Terminal_Event_Time = rep(rexp(n * max(recurrent_event_counts), rate = 0.05)),
  Terminal_Event_Indicator = rep(rbinom(n * max(recurrent_event_counts), size = 1, prob = 0.8)),
  Covariate_1 = rep(covariate_1, times = recurrent_event_counts),
  Covariate_2 = rep(covariate_2, times = recurrent_event_counts)
)

# Create the event indicator variable
data$Event_Indicator <- with(data, ifelse(
  Recurrent_Event_Count > 0 | Terminal_Event_Indicator == 1, 1, 0
))

# Create the survival object
surv_object <- Surv(
  time = pmin(data$Recurrent_Event_Time, data$Terminal_Event_Time),
  event = data$Event_Indicator
)

# Fit the RSF model
rsf_model <- rfsrc(
  surv_object ~ .,
  data = data,
  nsplit = 5,  # Number of splits for RSF
  ntree = 500  # Number of trees
)

# Display the RSF model summary
summary(rsf_model)
